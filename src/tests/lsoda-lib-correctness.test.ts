import {Lsoda} from '../solver-tools/lsoda';
import {ODEs} from '../../index';
import {corrProbs} from '../examples/corr-probs';

// ─── Robertson chemical kinetics (original test) ────────────────────────────

function fex(t: number, y: Float64Array, ydot: Float64Array): number {
  ydot[0] = 1.0e4 * y[1] * y[2] - 0.04 * y[0];
  ydot[2] = 3.0e7 * y[1] * y[1];
  ydot[1] = -(ydot[0] + ydot[2]);
  return 0;
}

const expected = [
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

// ─── Helpers ─────────────────────────────────────────────────────────────────

function relError(computed: number, reference: number): number {
  if (reference === 0) return Math.abs(computed);
  return Math.abs((computed - reference) / reference);
}

/** Wraps an ODEs func into an OdeFunction compatible with Lsoda (returns 0). */
function wrapFunc(ode: ODEs): (t: number, y: Float64Array, ydot: Float64Array) => number {
  return (t, y, ydot) => {
    ode.func(t, y, ydot);
    return 0;
  };
}

/**
 * Solves an ODE problem and checks against the exact solution at the final time point.
 * Uses combined absolute + relative tolerance: |computed - exact| < checkTol * |exact| + checkAbsTol.
 */
function solveAndCheck(
  ode: ODEs,
  exactSolution: (t: number) => Float64Array,
  checkTol: number = 1e-3,
  checkAbsTol: number = 1e-6,
) {
  const neq = ode.initial.length;
  const tol = ode.tolerance;
  const rtol = new Float64Array(neq).fill(tol);
  const atol = new Float64Array(neq).fill(tol);

  const solver = new Lsoda(wrapFunc(ode), neq, {rtol, atol, itask: 1});

  let y: ArrayLike<number> = [...ode.initial];
  let t = ode.arg.start;
  const step = ode.arg.step;
  const finish = ode.arg.finish;

  let tout = t + step;
  while (tout <= finish + step / 2) {
    const result = solver.solve(y, t, tout);
    y = result.y;
    t = result.t;

    expect(solver.state).toBeGreaterThan(0);
    tout += step;
  }

  // Check final solution against exact using combined abs + rel tolerance
  const exact = exactSolution(t);
  for (let j = 0; j < neq; j++) {
    const absErr = Math.abs(y[j] - exact[j]);
    const allowed = checkTol * Math.abs(exact[j]) + checkAbsTol;
    expect(absErr).toBeLessThan(allowed);
  }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

const nonStiffProbs = corrProbs.filter((p) => p.odes.name.startsWith('non-stiff'));
const stiffProbs = corrProbs.filter((p) => p.odes.name.startsWith('stiff'));

describe('LSODA solver', () => {
  it('should solve Robertson chemical kinetics with expected accuracy', () => {
    const neq = 3;
    const rtol = new Float64Array([1.0e-4, 1.0e-4, 1.0e-4]);
    const atol = new Float64Array([1.0e-6, 1.0e-10, 1.0e-6]);

    const solver = new Lsoda(fex, neq, {
      rtol,
      atol,
      itask: 1,
    });

    let y: ArrayLike<number> = [1.0, 0.0, 0.0];
    let t = 0.0;
    let tout = 0.4;

    for (let iout = 0; iout < 12; iout++) {
      const result = solver.solve(y, t, tout);
      y = result.y;
      t = result.t;

      expect(solver.state).toBeGreaterThan(0);

      const exp = expected[iout];
      for (let j = 0; j < 3; j++) {
        const err = relError(y[j], exp.y[j]);
        expect(err).toBeLessThanOrEqual(5e-3);
      }

      tout *= 10.0;
    }
  });

  describe('non-stiff problems', () => {
    nonStiffProbs.forEach(({odes, exact}) => {
      it(`should solve ${odes.name} problem`, () => {
        solveAndCheck(odes, exact);
      });
    });
  });

  describe('stiff problems', () => {
    stiffProbs.forEach(({odes, exact}) => {
      it(`should solve ${odes.name} problem`, () => {
        solveAndCheck(odes, exact);
      });
    });
  });

  // ─── Dense output tests ───────────────────────────────────────────────────

  describe('dense output', () => {
    // Helper: run Robertson to completion with dense output enabled
    function solveRobertsonDense() {
      const neq = 3;
      const rtol = new Float64Array([1.0e-4, 1.0e-4, 1.0e-4]);
      const atol = new Float64Array([1.0e-6, 1.0e-10, 1.0e-6]);

      const solver = new Lsoda(fex, neq, {rtol, atol, itask: 1, dense: true});

      let y: ArrayLike<number> = [1.0, 0.0, 0.0];
      let t = 0.0;
      let tout = 0.4;

      for (let iout = 0; iout < 12; iout++) {
        const result = solver.solve(y, t, tout);
        y = result.y;
        t = result.t;
        tout *= 10.0;
      }

      return solver.getDenseOutput();
    }

    it('consistency: dense output should match direct solver results (Robertson)', () => {
      const dense = solveRobertsonDense();

      // Query dense output at the same tout points used in the reference test
      const tPoints = Float64Array.from(expected.map((e) => e.t));
      const results = dense.solveAtTimes(tPoints);

      for (let i = 0; i < expected.length; i++) {
        for (let j = 0; j < 3; j++) {
          const err = relError(results[j][i], expected[i].y[j]);
          expect(err).toBeLessThanOrEqual(5e-3);
        }
      }
    });

    it('monotonicity: uniform query should have increasing t and finite y', () => {
      const dense = solveRobertsonDense();

      // Query at 1000 uniform points across the full range
      const tMin = dense.tMin;
      const tMax = dense.tMax;
      const n = 1000;
      const tArray = new Float64Array(n);
      for (let i = 0; i < n; i++)
        tArray[i] = tMin + (tMax - tMin) * i / (n - 1);

      const results = dense.solveAtTimes(tArray);
      expect(results.length).toBe(3);
      for (let j = 0; j < 3; j++)
        expect(results[j].length).toBe(n);

      for (let i = 1; i < n; i++)
        expect(tArray[i]).toBeGreaterThan(tArray[i - 1]);

      for (let i = 0; i < n; i++) {
        for (let j = 0; j < 3; j++)
          expect(Number.isFinite(results[j][i])).toBe(true);
      }
    });

    it('boundary: query at tMin and tMax should succeed', () => {
      const dense = solveRobertsonDense();

      const yMin = dense.evaluateAt(dense.tMin);
      const yMax = dense.evaluateAt(dense.tMax);

      expect(yMin.length).toBe(3);
      expect(yMax.length).toBe(3);
      for (let j = 0; j < 3; j++) {
        expect(Number.isFinite(yMin[j])).toBe(true);
        expect(Number.isFinite(yMax[j])).toBe(true);
      }
    });

    it('solveOnGrid: should return correct grid dimensions and quality', () => {
      const dense = solveRobertsonDense();

      const {t, y} = dense.solveOnGrid(0, 4e10, 1e9);

      expect(t.length).toBe(41);
      expect(y.length).toBe(3);
      for (let j = 0; j < 3; j++)
        expect(y[j].length).toBe(41);
      expect(t[0]).toBeCloseTo(0, 10);
      expect(t[t.length - 1]).toBeCloseTo(4e10, 0);

      // Check solution quality at a sampled point (t=4e10)
      const refLast = expected[expected.length - 1].y;
      for (let j = 0; j < 3; j++) {
        const err = relError(y[j][t.length - 1], refLast[j]);
        expect(err).toBeLessThanOrEqual(5e-3);
      }
    });

    function denseExactCheck(
      ode: ODEs,
      exactSolution: (t: number) => Float64Array,
      checkTol: number = 1e-3,
      checkAbsTol: number = 1e-6,
    ) {
      const neq = ode.initial.length;
      const tol = ode.tolerance;
      const rtol = new Float64Array(neq).fill(tol);
      const atol = new Float64Array(neq).fill(tol);

      const solver = new Lsoda(wrapFunc(ode), neq, {rtol, atol, itask: 1, dense: true});

      let y: ArrayLike<number> = [...ode.initial];
      let t = ode.arg.start;
      const finish = ode.arg.finish;
      const result = solver.solve(y, t, finish);
      y = result.y;
      t = result.t;

      const dense = solver.getDenseOutput();

      // Query at 100 points and compare against exact solution
      const nPts = 100;
      const tArr = new Float64Array(nPts + 1);
      for (let i = 0; i <= nPts; i++)
        tArr[i] = ode.arg.start + (finish - ode.arg.start) * i / nPts;

      const results = dense.solveAtTimes(tArr);
      for (let i = 0; i < tArr.length; i++) {
        const exact = exactSolution(tArr[i]);
        for (let j = 0; j < neq; j++) {
          const absErr = Math.abs(results[j][i] - exact[j]);
          const allowed = checkTol * Math.abs(exact[j]) + checkAbsTol;
          expect(absErr).toBeLessThan(allowed);
        }
      }
    }

    corrProbs.forEach(({odes, exact}) => {
      it(`dense output on ${odes.name} problem should match exact solution`, () => {
        denseExactCheck(odes, exact);
      });
    });
  });
});
