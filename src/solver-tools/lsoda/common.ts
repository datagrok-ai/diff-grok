// Constants from common.h
export const ETA = 2.2204460492503131e-16;
export const SQRTETA = 1.4901161193847656e-08;
export const CCMAX = 0.3;
export const MAXCOR = 3;
export const MSBP = 20;
export const MXNCF = 10;
export const RATIO = 5.0;

// Stability region limits for Adams method (from common.c)
export const sm1: number[] = [
  0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025,
];

// Method switch coefficients (from methodswitch.c)
// cm1 and cm2 are calculated from tesco/elco 1 and 2 formula in cfode.c
export const cm1: number[] = [
  0, 2.0, 5.999999999999998, 4.0,
  1.5780346820809881, 0.44444444444444464, 0.09712509625876042, 0.017636684303350972,
  0.002666977809498494, 0.00033760591160674403, 0.000035714285714285924, 3.1001984126984256e-06,
  2.1543369643350994e-07,
];

export const cm2: number[] = [
  0, 2.0, 1.5, 0.6666666666666667,
  0.20833333333333348, 0.04999999999999998, 0.09712509625876042, 0.017636684303350972,
  0.002666977809498494, 0.00033760591160674403, 0.000035714285714285924, 3.1001984126984256e-06,
  2.1543369643350994e-07,
];

// User-provided ODE function type
// t: current time, y: state vector (0-based), ydot: output derivatives (0-based), data: user data
export type LsodaFunc = (t: number, y: Float64Array, ydot: Float64Array, data: any) => number;

// Solver options
export interface LsodaOpt {
  ixpr: number; // print flag
  mxstep: number; // max internal steps per call
  mxhnil: number; // max step size warnings
  mxordn: number; // max Adams order (1-12)
  mxords: number; // max BDF order (1-5)
  tcrit: number; // critical time for itask 4,5
  h0: number; // initial step size (0 = auto)
  hmax: number; // max absolute step size
  hmin: number; // min absolute step size
  hmxi: number; // 1/hmax (computed)
  itask: number; // integration task (1-5)
  rtol: Float64Array; // relative tolerance (1-indexed)
  atol: Float64Array; // absolute tolerance (1-indexed)
}

// Internal solver state (mirrors lsoda_common_t)
export class LsodaCommon {
  // Solution history (2D: yh[j][i], 1-indexed)
  yh: Float64Array[] = [];
  // Working matrix for Jacobian/LU (2D: wm[i][j], 1-indexed)
  wm: Float64Array[] = [];
  // Error weights (1-indexed)
  ewt: Float64Array = new Float64Array(0);
  // Saved function values (1-indexed)
  savf: Float64Array = new Float64Array(0);
  // Accumulated corrections (1-indexed)
  acor: Float64Array = new Float64Array(0);
  // Pivot indices for LU (1-indexed)
  ipvt: Int32Array = new Int32Array(0);

  // Step size
  h: number = 0;
  hu: number = 0;
  rc: number = 0;
  tn: number = 0;
  tsw: number = 0;
  pdnorm: number = 0;

  // Convergence
  crate: number = 0;
  el: Float64Array = new Float64Array(14);
  hold: number = 0;
  rmax: number = 0;
  pdest: number = 0;
  pdlast: number = 0;

  // Method coefficients (2D arrays, 1-indexed)
  elco: Float64Array[] = [];
  tesco: Float64Array[] = [];

  // Counters and flags
  ialth: number = 0;
  ipup: number = 0;
  nslp: number = 0;
  icount: number = 0;
  irflag: number = 0;
  imxer: number = 0;
  illin: number = 0;
  nhnil: number = 0;
  nslast: number = 0;
  jcur: number = 0;
  meth: number = 0;
  mused: number = 0;
  nq: number = 0;
  nst: number = 0;
  ncf: number = 0;
  nfe: number = 0;
  nje: number = 0;
  nqu: number = 0;
  miter: number = 0;
}

// Snapshot of the Nordsieck array at a successful step (0-based arrays)
export interface NordsieckSnapshot {
  tn: number; // right edge of the step
  h: number; // step size used
  hu: number; // last successful step size (for range check)
  nq: number; // method order
  yh: Float64Array[]; // yh[0..nq], each of length neq (0-based copies)
}

// Solver context
export class LsodaContext {
  func: LsodaFunc;
  data: any;
  neq: number;
  state: number;
  error: string | null = null;
  common: LsodaCommon | null = null;
  opt: LsodaOpt | null = null;
  snapshots: NordsieckSnapshot[] | null = null;

  constructor(func: LsodaFunc, neq: number, data?: any) {
    this.func = func;
    this.neq = neq;
    this.data = data ?? null;
    this.state = 1;
  }
}
