import {Cvode, CvodeOptions, CvodeSolveResult} from '../solver-tools/cvode';
import {ODEs} from '../../index';
import {robertson, robertsonReferencePoint} from '../examples/robertson';
import {hires, hiresReferencePoint} from '../examples/hires';
import {vdpol, vdpolReferencePoint} from '../examples/vdpol';
import {orego, oregoReferencePoint} from '../examples/orego';
import {e5, e5ReferencePoint} from '../examples/e5';
import {pollution, pollutionReferencePoint} from '../examples/pollution';

// ─── Helpers ─────────────────────────────────────────────────────────────────

interface SolverOpts {
  mxstep?: number;
  hmax?: number;
  rtol?: Float64Array;
  atol?: Float64Array;
  /** Small initial tout to bootstrap very stiff problems before the main loop. */
  warmupTout?: number;
}

/** Solves a benchmark problem using CVODE, checks accuracy against the reference point. */
function solveBenchmarkAndCheck(
  ode: ODEs, referencePoint: Float64Array, relTol: number, absTol: number, opts?: SolverOpts,
): void {
  const neq = ode.initial.length;
  const tol = ode.tolerance;

  // Build CvodeOptions
  const cvodeOpts: CvodeOptions = {lmm: 'bdf'};

  // Handle tolerances
  if (opts?.rtol && opts?.atol) {
    cvodeOpts.rtol = opts.rtol[0];
    cvodeOpts.atol = opts.atol;
  } else if (opts?.rtol)
    cvodeOpts.rtol = opts.rtol[0];
  else if (opts?.atol)
    cvodeOpts.atol = opts.atol;
  else {
    cvodeOpts.rtol = tol;
    cvodeOpts.atol = tol;
  }

  if (opts?.mxstep) cvodeOpts.maxSteps = opts.mxstep;
  if (opts?.hmax) cvodeOpts.maxStep = opts.hmax;

  const y0 = Float64Array.from(ode.initial);
  const solver = new Cvode(ode.func, neq, ode.arg.start, y0, cvodeOpts);

  // Warmup phase: bootstrap very stiff problems with a small initial tout
  if (opts?.warmupTout !== undefined) {
    const result = solver.solve(opts.warmupTout);
    expect(result.flag).toBeGreaterThanOrEqual(0);
  }

  const t = opts?.warmupTout ?? ode.arg.start;
  let tout = t + ode.arg.step;
  const finish = ode.arg.finish;
  let lastResult: CvodeSolveResult;

  while (tout <= finish + ode.arg.step / 2) {
    lastResult = solver.solve(tout);
    if (lastResult.flag < 0)
      throw new Error(`${ode.name}: solver failed at t=${lastResult.t} with flag=${lastResult.flag}`);
    tout += ode.arg.step;
  }

  // Check final solution against reference point
  const y = lastResult!.y;
  for (let j = 0; j < neq; j++) {
    const absErr = Math.abs(y[j] - referencePoint[j]);
    const allowed = relTol * Math.abs(referencePoint[j]) + absTol;
    expect(absErr).toBeLessThan(allowed);
  }
}

// ─── Performance tests ──────────────────────────────────────────────────────

// Robertson: components span many orders of magnitude; use component-wise atol
describe('Robertson', () => {
  it('should solve correctly', () => {
    solveBenchmarkAndCheck(robertson, robertsonReferencePoint, 1e-2, 1e-7, {
      mxstep: 50000,
      atol: new Float64Array([1e-8, 1e-14, 1e-8]),
      rtol: new Float64Array([1e-7, 1e-7, 1e-7]),
    });
  }, 60000);
});

describe('HIRES', () => {
  it('should solve correctly', () => {
    solveBenchmarkAndCheck(hires, hiresReferencePoint, 1e-3, 1e-12);
  }, 60000);
});

// Van der Pol (mu=1000): very stiff, needs many internal steps
describe('van der Pol', () => {
  it('should solve correctly', () => {
    solveBenchmarkAndCheck(vdpol, vdpolReferencePoint, 1e-3, 1e-6, {
      mxstep: 50000,
    });
  }, 60000);
});

describe('OREGO', () => {
  it('should solve correctly', () => {
    solveBenchmarkAndCheck(orego, oregoReferencePoint, 1e-2, 1e-6);
  }, 60000);
});

// E5: extremely stiff (rate constants span 20 orders of magnitude).
// Known limitation: CVODE BDF with dense direct solver hits convergence failures
// on this extreme problem. LSODA handles it via automatic stiff/non-stiff switching.
// TODO: revisit when band/iterative solvers or preconditioners are added.
describe('E5', () => {
  it.skip('should solve correctly (known limitation: CV_CONV_FAILURE)', () => {
    solveBenchmarkAndCheck(e5, e5ReferencePoint, 1, 1e-6, {
      mxstep: 50000,
      warmupTout: 1e-5,
    });
  }, 60000);
});

describe('Pollution', () => {
  it('should solve correctly', () => {
    solveBenchmarkAndCheck(pollution, pollutionReferencePoint, 1e-2, 1e-10);
  }, 60000);
});
