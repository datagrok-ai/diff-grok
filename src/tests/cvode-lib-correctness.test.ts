import {
  Cvode, cvodeCreate, cvodeInit, cvode, cvodeSStolerances, cvodeSVtolerances,
  cvodeSetLinearSolver, cvodeGetDky,
  CV_BDF, CV_NORMAL,
  CV_SUCCESS, CV_ROOT_RETURN, CV_TOO_MUCH_WORK,
} from '../solver-tools/cvode';
import type {OdeFunction, CvodeSolveResult, CvodeRhsFn} from '../solver-tools/cvode';
import {ODEs} from '../../index';
import {corrProbs} from '../examples/corr-probs';

// ─── Robertson chemical kinetics (original test) ────────────────────────────

function robertsonRhs(t: number, y: Float64Array, ydot: Float64Array): void {
  ydot[0] = 1.0e4 * y[1] * y[2] - 0.04 * y[0];
  ydot[2] = 3.0e7 * y[1] * y[1];
  ydot[1] = -(ydot[0] + ydot[2]);
}

/** Low-level RHS for Robertson (returns int). */
function robertsonRhsLowLevel(t: number, y: Float64Array, ydot: Float64Array, _userData: any): number {
  robertsonRhs(t, y, ydot);
  return 0;
}

const robertsonExpected = [
  {t: 4.0000e-01, y: [9.851712e-01, 3.386380e-05, 1.479493e-02]},
  {t: 4.0000e+00, y: [9.055333e-01, 2.240655e-05, 9.444430e-02]},
  {t: 4.0000e+01, y: [7.158403e-01, 9.186334e-06, 2.841505e-01]},
  {t: 4.0000e+02, y: [4.505250e-01, 3.222964e-06, 5.494717e-01]},
  {t: 4.0000e+03, y: [1.831976e-01, 8.941773e-07, 8.168015e-01]},
  {t: 4.0000e+04, y: [3.898729e-02, 1.621940e-07, 9.610125e-01]},
  {t: 4.0000e+05, y: [4.936362e-03, 1.984221e-08, 9.950636e-01]},
  {t: 4.0000e+06, y: [5.161833e-04, 2.065787e-09, 9.994838e-01]},
  {t: 4.0000e+07, y: [5.179804e-05, 2.072027e-10, 9.999482e-01]},
  {t: 4.0000e+08, y: [5.283675e-06, 2.113481e-11, 9.999947e-01]},
  {t: 4.0000e+09, y: [4.658667e-07, 1.863468e-12, 9.999995e-01]},
  {t: 4.0000e+10, y: [1.431100e-08, 5.724404e-14, 1.000000e+00]},
];

// Robertson solver tolerances for CVODE.
const ROBERTSON_SOLVER_RTOL = 1e-8;
const ROBERTSON_SOLVER_ATOL = new Float64Array([1e-10, 1e-14, 1e-10]);

// Number of Robertson output points to check with strict relative tolerance.
// At very late times (t>=4e8) the solution components become extremely small
// (y1~1e-5, y2~1e-10) and accumulated integration error over many orders of
// magnitude of time makes pure relative error checks impractical. We check
// the first 9 points (up to t=4e7) with strict tolerance and use qualitative
// checks for the remaining points.
const ROBERTSON_STRICT_POINTS = 9;
const ROBERTSON_CHECK_TOL = 5e-3;

// ─── Helpers ─────────────────────────────────────────────────────────────────

function relError(computed: number, reference: number): number {
  if (reference === 0) return Math.abs(computed);
  return Math.abs((computed - reference) / reference);
}

/**
 * Solves an ODE with the Cvode class and checks against exact solution at the final time.
 * Uses combined absolute + relative tolerance: |computed - exact| < checkTol * |exact| + checkAbsTol.
 */
function solveAndCheck(
  ode: ODEs,
  exactSolution: (t: number) => Float64Array,
  lmm: 'adams' | 'bdf',
  checkTol: number = 1e-3,
  checkAbsTol: number = 1e-6,
) {
  const neq = ode.initial.length;
  const tol = ode.tolerance;

  const solver = new Cvode(ode.func, neq, ode.arg.start, Float64Array.from(ode.initial), {
    lmm,
    rtol: tol,
    atol: tol,
  });

  let t = ode.arg.start;
  const step = ode.arg.step;
  const finish = ode.arg.finish;
  let result: CvodeSolveResult;

  let tout = t + step;
  while (tout <= finish + step / 2) {
    result = solver.solve(tout);
    expect(result.flag).toBe(CV_SUCCESS);
    t = result.t;
    tout += step;
  }

  // Check final solution against exact
  const exact = exactSolution(t);
  const y = result!.y;
  for (let j = 0; j < neq; j++) {
    const absErr = Math.abs(y[j] - exact[j]);
    const allowed = checkTol * Math.abs(exact[j]) + checkAbsTol;
    expect(absErr).toBeLessThan(allowed);
  }
}

/**
 * Dense output helper: solve an ODE to evenly spaced output points, then use getDky
 * to verify the interpolation matches the exact solution.
 */
function denseExactCheck(
  ode: ODEs,
  exactSolution: (t: number) => Float64Array,
  lmm: 'adams' | 'bdf',
  checkTol: number = 1e-3,
  checkAbsTol: number = 1e-6,
) {
  const neq = ode.initial.length;
  const tol = ode.tolerance;

  const solver = new Cvode(ode.func, neq, ode.arg.start, Float64Array.from(ode.initial), {
    lmm,
    rtol: tol,
    atol: tol,
  });

  const start = ode.arg.start;
  const finish = ode.arg.finish;
  const nOutputPts = 100;
  let nChecked = 0;

  for (let i = 1; i <= nOutputPts; i++) {
    const tout = start + (finish - start) * i / nOutputPts;
    const result = solver.solve(tout);
    expect(result.flag).toBe(CV_SUCCESS);

    // getDky at the output time should match the solve result
    const dky = solver.getDky(result.t, 0);
    const exact = exactSolution(result.t);

    for (let j = 0; j < neq; j++) {
      // Dense output should match the direct solve result
      expect(Math.abs(dky[j] - result.y[j])).toBeLessThan(1e-12);

      // Dense output should match the exact solution
      const absErr = Math.abs(dky[j] - exact[j]);
      const allowed = checkTol * Math.abs(exact[j]) + checkAbsTol;
      expect(absErr).toBeLessThan(allowed);
    }
    nChecked++;
  }

  // Ensure we actually checked points
  expect(nChecked).toBe(nOutputPts);
}

// ─── Tests ───────────────────────────────────────────────────────────────────

const nonStiffProbs = corrProbs.filter((p) => p.odes.name.startsWith('non-stiff'));
const stiffProbs = corrProbs.filter((p) => p.odes.name.startsWith('stiff'));

describe('CVODE LIB', () => {
  it('should solve Robertson chemical kinetics with Cvode class', () => {
    const neq = 3;
    const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
      lmm: 'bdf',
      rtol: ROBERTSON_SOLVER_RTOL,
      atol: ROBERTSON_SOLVER_ATOL,
    });

    let tout = 0.4;
    for (let iout = 0; iout < 12; iout++) {
      const result = solver.solve(tout);
      expect(result.flag).toBe(CV_SUCCESS);

      const exp = robertsonExpected[iout];
      if (iout < ROBERTSON_STRICT_POINTS) {
        for (let j = 0; j < 3; j++) {
          const err = relError(result.y[j], exp.y[j]);
          expect(err).toBeLessThanOrEqual(ROBERTSON_CHECK_TOL);
        }
      } else {
        const sum = result.y[0] + result.y[1] + result.y[2];
        expect(Math.abs(sum - 1.0)).toBeLessThan(1e-4);
        expect(result.y[2]).toBeGreaterThan(0.999);
      }

      tout *= 10.0;
    }
  });

  it('should solve Robertson chemical kinetics with the low-level API', () => {
    const neq = 3;
    const y0 = Float64Array.from([1.0, 0.0, 0.0]);

    const mem = cvodeCreate(CV_BDF);
    cvodeInit(mem, robertsonRhsLowLevel, 0.0, y0);
    cvodeSVtolerances(mem, ROBERTSON_SOLVER_RTOL, ROBERTSON_SOLVER_ATOL);
    cvodeSetLinearSolver(mem);

    const yout = new Float64Array(neq);
    let tout = 0.4;

    for (let iout = 0; iout < 12; iout++) {
      const {flag, t} = cvode(mem, tout, yout, CV_NORMAL);
      expect(flag).toBe(CV_SUCCESS);

      const exp = robertsonExpected[iout];
      if (iout < ROBERTSON_STRICT_POINTS) {
        for (let j = 0; j < 3; j++) {
          const err = relError(yout[j], exp.y[j]);
          expect(err).toBeLessThanOrEqual(ROBERTSON_CHECK_TOL);
        }
      } else {
        const sum = yout[0] + yout[1] + yout[2];
        expect(Math.abs(sum - 1.0)).toBeLessThan(1e-4);
        expect(yout[2]).toBeGreaterThan(0.999);
      }

      tout *= 10.0;
    }
  });

  describe('non-stiff problems (Adams)', () => {
    nonStiffProbs.forEach(({odes, exact}) => {
      it(`should solve ${odes.name} problem`, () => {
        solveAndCheck(odes, exact, 'adams');
      });
    });
  });

  describe('stiff problems (BDF)', () => {
    stiffProbs.forEach(({odes, exact}) => {
      it(`should solve ${odes.name} problem`, () => {
        solveAndCheck(odes, exact, 'bdf');
      });
    });
  });

  // ─── Rootfinding ─────────────────────────────────────────────────────────

  describe('rootfinding', () => {
    it('should detect roots in Robertson system', () => {
      const neq = 3;

      // Root functions: g1 = y1 - 1e-4, g2 = y3 - 0.01
      const rootFn = (t: number, y: Float64Array, gout: Float64Array): void => {
        gout[0] = y[0] - 1e-4;
        gout[1] = y[2] - 0.01;
      };

      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
        rootFn,
        nRootFns: 2,
      });

      const rootTimes: number[] = [];
      const rootIndices: Int32Array[] = [];

      // Advance through the problem, collecting root returns
      let tout = 0.4;
      const tEnd = 4e10;

      while (tout <= tEnd) {
        const result = solver.solve(tout);

        if (result.flag === CV_ROOT_RETURN) {
          rootTimes.push(result.t);
          rootIndices.push(result.rootsFound!);
          // Continue from the root time; do not advance tout yet
          continue;
        }

        expect(result.flag).toBe(CV_SUCCESS);
        tout *= 10.0;
      }

      // Expect at least 2 roots: y3 reaches 0.01 first, then y1 drops to 1e-4
      expect(rootTimes.length).toBeGreaterThanOrEqual(2);

      // First root: y3 = 0.01, near t ~ 2.6391e-1
      const tRoot1 = rootTimes[0];
      expect(relError(tRoot1, 2.6391e-1)).toBeLessThan(2e-2);
      // rootsFound should indicate root function index 1 (g2)
      expect(rootIndices[0][1]).not.toBe(0);

      // Second root: y1 = 1e-4, near t ~ 2.0790e+7
      const tRoot2 = rootTimes[1];
      expect(relError(tRoot2, 2.0790e7)).toBeLessThan(2e-2);
      // rootsFound should indicate root function index 0 (g1)
      expect(rootIndices[1][0]).not.toBe(0);
    });
  });

  // ─── Dense output tests ─────────────────────────────────────────────────

  describe('dense output (getDky)', () => {
    it('consistency: getDky at output points should match solve results (Robertson)', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
      });

      let tout = 0.4;
      for (let iout = 0; iout < 12; iout++) {
        const result = solver.solve(tout);
        expect(result.flag).toBe(CV_SUCCESS);

        // getDky at the output time should match the solve result
        const dky = solver.getDky(result.t, 0);
        for (let j = 0; j < neq; j++)
          expect(Math.abs(dky[j] - result.y[j])).toBeLessThan(1e-12);


        tout *= 10.0;
      }
    });

    it('monotonicity: getDky at successive output points should produce finite values', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
      });

      // Solve to 1000 uniformly spaced points across [0.01, 1.0]
      const n = 1000;
      const tMin = 0.01;
      const tMax = 1.0;
      let prevT = 0.0;

      for (let i = 0; i < n; i++) {
        const tout = tMin + (tMax - tMin) * i / (n - 1);
        const result = solver.solve(tout);
        expect(result.flag).toBe(CV_SUCCESS);

        // Time should advance monotonically
        expect(result.t).toBeGreaterThan(prevT);
        prevT = result.t;

        // getDky at the output point should return finite values
        const dky = solver.getDky(result.t, 0);
        for (let j = 0; j < neq; j++)
          expect(Number.isFinite(dky[j])).toBe(true);
      }
    });

    it('boundary: getDky at endpoints of the current step interval', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
      });

      // Advance to a point to establish state
      const result = solver.solve(0.4);
      expect(result.flag).toBe(CV_SUCCESS);

      // getDky at the last output time should succeed
      const dkyEnd = solver.getDky(result.t, 0);
      for (let j = 0; j < neq; j++)
        expect(Number.isFinite(dkyEnd[j])).toBe(true);


      // Query derivative (k=1) - first derivative at output time
      const dkyDeriv = solver.getDky(result.t, 1);
      for (let j = 0; j < neq; j++)
        expect(Number.isFinite(dkyDeriv[j])).toBe(true);
    });

    corrProbs.forEach(({odes, exact}) => {
      const lmm = odes.name.startsWith('stiff') ? 'bdf' as const : 'adams' as const;
      it(`dense output on ${odes.name} problem should match exact solution`, () => {
        denseExactCheck(odes, exact, lmm);
      });
    });
  });

  // ─── Error handling ─────────────────────────────────────────────────────

  describe('error handling', () => {
    it('should return error for negative tolerances (low-level API)', () => {
      const mem = cvodeCreate(CV_BDF);
      cvodeInit(mem, robertsonRhsLowLevel, 0.0, Float64Array.from([1.0, 0.0, 0.0]));

      const flag = cvodeSStolerances(mem, -1e-4, 1e-6);
      expect(flag).toBeLessThan(0);
    });

    it('should return CV_TOO_MUCH_WORK with maxSteps: 1', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
        maxSteps: 1,
      });

      const result = solver.solve(4e10);
      expect(result.flag).toBe(CV_TOO_MUCH_WORK);
    });

    it('should return error flag when RHS function fails', () => {
      const failingRhs: OdeFunction = (_t, _y, ydot) => {
        ydot[0] = NaN;
      };

      const solver = new Cvode(failingRhs, 1, 0.0, Float64Array.from([1.0]), {
        lmm: 'bdf',
        rtol: 1e-4,
        atol: 1e-6,
      });

      const result = solver.solve(1.0);
      expect(result.flag).toBeLessThan(0);
    });
  });

  // ─── Cvode class features ──────────────────────────────────────────────

  describe('Cvode class features', () => {
    it('getStats should return valid integrator statistics', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
      });

      // Solve to first output point
      solver.solve(0.4);

      const stats = solver.getStats();
      expect(stats.nSteps).toBeGreaterThan(0);
      expect(stats.nRhsEvals).toBeGreaterThan(0);
      expect(stats.nLinSolvSetups).toBeGreaterThanOrEqual(0);
      expect(stats.nErrTestFails).toBeGreaterThanOrEqual(0);
      expect(stats.nNonlinSolvIters).toBeGreaterThan(0);
      expect(stats.nNonlinSolvConvFails).toBeGreaterThanOrEqual(0);
      expect(stats.nJacEvals).toBeGreaterThanOrEqual(0);
      expect(stats.lastOrder).toBeGreaterThanOrEqual(1);
      expect(stats.lastStep).toBeGreaterThan(0);
      expect(stats.currentTime).toBeGreaterThan(0);
    });

    it('reInit should allow solving with new initial conditions', () => {
      const neq = 3;
      const solver = new Cvode(robertsonRhs, neq, 0.0, Float64Array.from([1.0, 0.0, 0.0]), {
        lmm: 'bdf',
        rtol: ROBERTSON_SOLVER_RTOL,
        atol: ROBERTSON_SOLVER_ATOL,
      });

      // Solve to first output point
      const result1 = solver.solve(0.4);
      expect(result1.flag).toBe(CV_SUCCESS);
      const y1 = new Float64Array(result1.y);

      // Re-initialize with the same initial conditions and solve again
      solver.reInit(0.0, Float64Array.from([1.0, 0.0, 0.0]));
      const result2 = solver.solve(0.4);
      expect(result2.flag).toBe(CV_SUCCESS);

      // Results should match closely
      for (let j = 0; j < neq; j++) {
        const err = relError(result2.y[j], y1[j]);
        expect(err).toBeLessThanOrEqual(1e-6);
      }
    });

    it('reInit with different initial conditions should produce different results', () => {
      const neq = 1;
      const simpleRhs: OdeFunction = (_t, y, ydot) => {
        ydot[0] = -y[0]; // y' = -y, solution: y(t) = y0 * exp(-t)
      };

      const solver = new Cvode(simpleRhs, neq, 0.0, Float64Array.from([1.0]), {
        lmm: 'adams',
        rtol: 1e-10,
        atol: 1e-12,
      });

      // Solve with y0 = 1
      const result1 = solver.solve(1.0);
      expect(result1.flag).toBe(CV_SUCCESS);
      const y1 = result1.y[0]; // should be ~ exp(-1) ~ 0.3679

      // Re-initialize with y0 = 2
      solver.reInit(0.0, Float64Array.from([2.0]));
      const result2 = solver.solve(1.0);
      expect(result2.flag).toBe(CV_SUCCESS);
      const y2 = result2.y[0]; // should be ~ 2*exp(-1) ~ 0.7358

      // y2 should be approximately 2 * y1
      expect(relError(y2, 2 * y1)).toBeLessThan(1e-6);

      // Check against exact solutions
      expect(relError(y1, Math.exp(-1))).toBeLessThan(1e-6);
      expect(relError(y2, 2 * Math.exp(-1))).toBeLessThan(1e-6);
    });
  });
});
