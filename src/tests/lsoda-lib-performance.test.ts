import {Lsoda} from '../solver-tools/lsoda';
import {ODEs} from '../../index';
import {robertson, robertsonReferencePoint} from '../examples/robertson';
import {hires, hiresReferencePoint} from '../examples/hires';
import {vdpol, vdpolReferencePoint} from '../examples/vdpol';
import {orego, oregoReferencePoint} from '../examples/orego';
import {e5, e5ReferencePoint} from '../examples/e5';
import {pollution, pollutionReferencePoint} from '../examples/pollution';

// ─── Helpers ─────────────────────────────────────────────────────────────────

/** Wraps an ODEs func into an OdeFunction compatible with Lsoda (returns 0). */
function wrapFunc(ode: ODEs): (t: number, y: Float64Array, ydot: Float64Array) => number {
  return (t, y, ydot) => {
    ode.func(t, y, ydot);
    return 0;
  };
}

interface SolverOpts {
  mxstep?: number;
  hmax?: number;
  rtol?: Float64Array;
  atol?: Float64Array;
  /** Small initial tout to bootstrap very stiff problems before the main loop. */
  warmupTout?: number;
}

/** Solves a benchmark problem, checks accuracy against the reference point, and verifies solver health. */
function solveBenchmarkAndCheck(
  ode: ODEs, referencePoint: Float64Array, relTol: number, absTol: number, opts?: SolverOpts,
): void {
  const neq = ode.initial.length;
  const tol = ode.tolerance;
  const rtol = opts?.rtol ?? new Float64Array(neq).fill(tol);
  const atol = opts?.atol ?? new Float64Array(neq).fill(tol);

  const solver = new Lsoda(wrapFunc(ode), neq, {
    rtol, atol, itask: 1,
    ...(opts?.mxstep ? {mxstep: opts.mxstep} : {}),
    ...(opts?.hmax ? {hmax: opts.hmax} : {}),
  });

  let y: ArrayLike<number> = [...ode.initial];
  let t = ode.arg.start;
  const step = ode.arg.step;
  const finish = ode.arg.finish;

  // Warmup phase: bootstrap very stiff problems with a small initial tout
  if (opts?.warmupTout !== undefined) {
    const result = solver.solve(y, t, opts.warmupTout);
    y = result.y;
    t = result.t;
    expect(solver.state).toBeGreaterThan(0);
  }

  let tout = t + step;
  while (tout <= finish + step / 2) {
    const result = solver.solve(y, t, tout);
    y = result.y;
    t = result.t;

    if (solver.state <= 0)
      throw new Error(`${ode.name}: solver failed at t=${t} with state=${solver.state}, error: ${solver.error}`);

    tout += step;
  }

  // Check final solution against reference point
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

// E5: extremely stiff (rate constants span 20 orders), needs warmup and many steps.
// Reference values are near-zero (y1~1e-290, y2~1e-22), so absTol must be generous.
describe('E5', () => {
  it('should solve correctly', () => {
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
