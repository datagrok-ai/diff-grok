/**
 * TypeScript port of SUNDIALS CVODE v7.5.0 (https://github.com/LLNL/sundials)
 *
 * Copyright (c) 2025, Lawrence Livermore National Security, University of
 *   Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security and
 *   Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * Copyright (c) 2026 Datagrok.
 *
 * Licensed under the BSD-3-Clause License. See THIRD_PARTY_LICENSES
 * for the full license text and provenance chain of the original code.
 */

// CVODE TypeScript port - Foundation types, constants, and data structures
// Ported from SUNDIALS CVODE v7.5.0 (cvode_impl.h, cvode_ls_impl.h)
// All arrays are 0-based.

// ============================================================================
// Constants from cvode_impl.h
// ============================================================================

// --- Method order limits ---
export const ADAMS_Q_MAX = 12;
export const BDF_Q_MAX = 5;
export const Q_MAX = 12;
export const L_MAX = 13;       // Q_MAX + 1
export const NUM_TESTS = 5;

// --- Linear multistep method types ---
export const CV_ADAMS = 1;
export const CV_BDF = 2;

// --- Task types ---
export const CV_NORMAL = 1;    // integrate to tout, interpolate
export const CV_ONE_STEP = 2;  // take one internal step

// --- Return values (success) ---
export const CV_SUCCESS = 0;
export const CV_TSTOP_RETURN = 1;
export const CV_ROOT_RETURN = 2;

// --- Return values (failure) ---
export const CV_TOO_MUCH_WORK = -1;
export const CV_TOO_MUCH_ACC = -2;
export const CV_ERR_FAILURE = -3;
export const CV_CONV_FAILURE = -4;
export const CV_RHSFUNC_FAIL = -8;
export const CV_FIRST_RHSFUNC_ERR = -9;
export const CV_RTFUNC_FAIL = -12;
export const CV_TOO_CLOSE = -27;

// --- Default values ---
export const HMIN_DEFAULT = 0.0;
export const HMAX_INV_DEFAULT = 0.0;
export const MXHNIL_DEFAULT = 10;
export const MXSTEP_DEFAULT = 500;
export const MSBP_DEFAULT = 20;
export const DGMAX_LSETUP_DEFAULT = 0.3;

// --- Eta (step size ratio) defaults ---
export const ETA_MIN_FX_DEFAULT = 0.0;
export const ETA_MAX_FX_DEFAULT = 1.5;
export const ETA_MAX_FS_DEFAULT = 10000.0;
export const ETA_MAX_ES_DEFAULT = 10.0;
export const ETA_MAX_GS_DEFAULT = 10.0;
export const ETA_MIN_DEFAULT = 0.1;
export const ETA_MAX_EF_DEFAULT = 0.2;
export const ETA_MIN_EF_DEFAULT = 0.1;
export const ETA_CF_DEFAULT = 0.25;
export const SMALL_NST_DEFAULT = 10;
export const SMALL_NEF_DEFAULT = 2;
export const ONEPSM = 1.000001;

// --- Step size selection bias factors ---
export const ADDON = 0.000001;
export const BIAS1 = 6.0;
export const BIAS2 = 6.0;
export const BIAS3 = 10.0;
export const LONG_WAIT = 10;

// --- Failure limits ---
export const MXNCF = 10;      // max convergence failures per step
export const MXNEF = 7;       // max error test failures per step
export const MXNEF1 = 3;      // max error test failures on first step

// --- Internal action flags ---
export const DO_ERROR_TEST = 2;
export const PREDICT_AGAIN = 3;
export const TRY_AGAIN = 5;
export const FIRST_CALL = 6;
export const PREV_CONV_FAIL = 7;
export const PREV_ERR_FAIL = 9;
export const RHSFUNC_RECVR = 10;

// --- Tolerance types ---
export const CV_NN = 0;       // no tolerances set
export const CV_SS = 1;       // scalar-scalar
export const CV_SV = 2;       // scalar-vector

// --- Nonlinear solver failure types ---
export const CV_NO_FAILURES = 0;
export const CV_FAIL_BAD_J = 1;
export const CV_FAIL_OTHER = 2;

// --- Initial step size estimation ---
export const FUZZ_FACTOR = 100.0;
export const HLB_FACTOR = 100.0;
export const HUB_FACTOR = 0.1;
export const H_BIAS = 0.5;
export const MAX_ITERS = 4;
export const CORTES = 0.1;
export const TINY = 1.0e-10;

// --- Nonlinear solver constants ---
export const NLS_MAXCOR = 3;  // max corrector iterations
export const CRDOWN = 0.3;    // convergence rate damping
export const RDIV = 2.0;      // divergence threshold

// --- SUNNonlinearSolver return codes (simplified) ---
export const SUN_NLS_CONTINUE = 901;
export const SUN_NLS_CONV_RECVR = 902;

// --- Additional failure codes ---
export const CV_LSETUP_FAIL = -13;
export const CV_LSOLVE_FAIL = -14;
export const CV_UNREC_RHSFUNC_ERR = -10;
export const CV_REPTD_RHSFUNC_ERR = -11;
export const CV_LINIT_FAIL = -31;
export const CV_NLS_INIT_FAIL = -32;
export const CV_NLS_FAIL = -33;

// --- Machine epsilon ---
export const UROUND = 2.2204460492503131e-16;

// --- Rootfinding constants ---
export const RTFOUND = 1;
export const CLOSERT = 3;

// --- Linear solver constants (from cvode_ls_impl.h) ---
export const CVLS_MSBJ = 51;  // max steps between Jacobian evals
export const CVLS_DGMAX = 0.2; // gamma change threshold for Jacobian redo

// ============================================================================
// Function types
// ============================================================================

/** RHS function: y' = f(t, y). Returns 0 on success, non-zero on failure. */
export type CvodeRhsFn = (t: number, y: Float64Array, ydot: Float64Array, userData: any) => number;

/** Jacobian function: J = df/dy. Returns 0 on success, non-zero on failure. */
export type CvodeJacFn = (t: number, y: Float64Array, fy: Float64Array, J: Float64Array[], userData: any) => number;

/** Root function: g(t, y). Returns 0 on success, non-zero on failure. */
export type CvodeRootFn = (t: number, y: Float64Array, gout: Float64Array, userData: any) => number;

// ============================================================================
// CvLsMem - Linear solver memory (maps CVLsMemRec from cvode_ls_impl.h)
// ============================================================================

export class CvLsMem {
  /** N x N Jacobian matrix (row-major) */
  A: Float64Array[] = [];
  /** LU pivot indices */
  pivots: Int32Array = new Int32Array(0);
  /** Saved copy of Jacobian */
  savedJ: Float64Array[] = [];
  /** User-supplied Jacobian function */
  jacFn: CvodeJacFn | null = null;
  /** Use difference-quotient Jacobian by default */
  jacDQ: boolean = true;
  /** Number of Jacobian evaluations */
  nje: number = 0;
  /** Number of f evaluations for DQ Jacobian */
  nfeDQ: number = 0;
  /** nst at last Jacobian evaluation */
  nstlj: number = 0;
  /** t at last Jacobian evaluation */
  tnlj: number = 0;
  /** Max steps between Jacobian evaluations (CVLS_MSBJ) */
  msbj: number = CVLS_MSBJ;
  /** Gamma change threshold (CVLS_DGMAX) */
  dgmax_jbad: number = CVLS_DGMAX;
  /** Is Jacobian stale? */
  jbad: boolean = true;
}

// ============================================================================
// CvodeMem - Main solver state (maps CVodeMemRec from cvode_impl.h)
// ============================================================================

export class CvodeMem {
  // --- Problem specification ---
  /** RHS function: y' = f(t, y) */
  cv_f: CvodeRhsFn | null = null;
  /** User data passed to f, jac, root functions */
  cv_user_data: any = null;
  /** Linear multistep method: CV_ADAMS or CV_BDF */
  cv_lmm: number = 0;
  /** Tolerance type: CV_NN, CV_SS, or CV_SV */
  cv_itol: number = CV_NN;

  // --- Tolerances ---
  /** Relative tolerance */
  cv_reltol: number = 0;
  /** Scalar absolute tolerance */
  cv_Sabstol: number = 0;
  /** Vector absolute tolerance */
  cv_Vabstol: Float64Array = new Float64Array(0);
  /** Flag: is min(atol) == 0? */
  cv_atolmin0: boolean = false;

  // --- Nordsieck history array ---
  /** Nordsieck array: cv_zn[j][i], j=0..q, i=0..N-1 */
  cv_zn: Float64Array[] = [];

  // --- Work vectors ---
  /** Error weight vector */
  cv_ewt: Float64Array = new Float64Array(0);
  /** Solution vector (scratch for user output) */
  cv_y: Float64Array = new Float64Array(0);
  /** Accumulated correction */
  cv_acor: Float64Array = new Float64Array(0);
  /** Temporary vector */
  cv_tempv: Float64Array = new Float64Array(0);
  /** Temporary for f evaluations */
  cv_ftemp: Float64Array = new Float64Array(0);
  /** Temporary vector 1 */
  cv_vtemp1: Float64Array = new Float64Array(0);
  /** Temporary vector 2 */
  cv_vtemp2: Float64Array = new Float64Array(0);
  /** Temporary vector 3 */
  cv_vtemp3: Float64Array = new Float64Array(0);

  // --- Tstop ---
  /** Is tstop set? */
  cv_tstopset: boolean = false;
  /** Interpolate at tstop? */
  cv_tstopinterp: boolean = false;
  /** Stop time */
  cv_tstop: number = 0;

  // --- Step data ---
  /** Current method order */
  cv_q: number = 0;
  /** Order for next step */
  cv_qprime: number = 0;
  /** Next order to be used (during step) */
  cv_next_q: number = 0;
  /** Steps remaining before order change allowed */
  cv_qwait: number = 0;
  /** L = q + 1 */
  cv_L: number = 0;
  /** Initial step size (user-supplied) */
  cv_hin: number = 0;
  /** Current step size */
  cv_h: number = 0;
  /** Step size for next step */
  cv_hprime: number = 0;
  /** Next h (during step) */
  cv_next_h: number = 0;
  /** Step size ratio: hprime / h */
  cv_eta: number = 0;
  /** Step size at last Nordsieck rescale */
  cv_hscale: number = 0;
  /** Current internal time */
  cv_tn: number = 0;
  /** Last return time */
  cv_tretlast: number = 0;
  /** Step size history: tau[i], i=0..L_MAX */
  cv_tau: Float64Array = new Float64Array(L_MAX + 1);
  /** Error test quantities: tq[i], i=0..NUM_TESTS */
  cv_tq: Float64Array = new Float64Array(NUM_TESTS + 1);
  /** Polynomial coefficients: l[i], i=0..L_MAX-1 */
  cv_l: Float64Array = new Float64Array(L_MAX);
  /** 1 / l[1] */
  cv_rl1: number = 0;
  /** gamma = h * rl1 */
  cv_gamma: number = 0;
  /** gamma at last setup */
  cv_gammap: number = 0;
  /** gamma / gammap */
  cv_gamrat: number = 0;
  /** NLS convergence rate estimate */
  cv_crate: number = 0;
  /** Previous del value (NLS) */
  cv_delp: number = 0;
  /** WRMS norm of accumulated correction */
  cv_acnrm: number = 0;
  /** Current acnrm (during step) */
  cv_acnrmcur: boolean = false;
  /** NLS convergence coefficient */
  cv_nlscoef: number = 0;

  // --- Limits ---
  /** Max method order */
  cv_qmax: number = 0;
  /** Max internal steps per call */
  cv_mxstep: number = MXSTEP_DEFAULT;
  /** Max t+h==t warnings */
  cv_mxhnil: number = MXHNIL_DEFAULT;
  /** Max error test failures per step */
  cv_maxnef: number = MXNEF;
  /** Max convergence failures per step */
  cv_maxncf: number = MXNCF;
  /** Min step size */
  cv_hmin: number = HMIN_DEFAULT;
  /** 1/hmax (inverse of max step size) */
  cv_hmax_inv: number = HMAX_INV_DEFAULT;
  /** Max eta for step size increase */
  cv_etamax: number = 0;
  /** Min fixed eta threshold */
  cv_eta_min_fx: number = ETA_MIN_FX_DEFAULT;
  /** Max fixed eta threshold */
  cv_eta_max_fx: number = ETA_MAX_FX_DEFAULT;
  /** Max eta on first step */
  cv_eta_max_fs: number = ETA_MAX_FS_DEFAULT;
  /** Max eta on early steps */
  cv_eta_max_es: number = ETA_MAX_ES_DEFAULT;
  /** Max eta on general steps */
  cv_eta_max_gs: number = ETA_MAX_GS_DEFAULT;
  /** Min eta */
  cv_eta_min: number = ETA_MIN_DEFAULT;
  /** Min eta on error failure */
  cv_eta_min_ef: number = ETA_MIN_EF_DEFAULT;
  /** Max eta on error failure */
  cv_eta_max_ef: number = ETA_MAX_EF_DEFAULT;
  /** Eta on convergence failure */
  cv_eta_cf: number = ETA_CF_DEFAULT;
  /** Small step count threshold */
  cv_small_nst: number = SMALL_NST_DEFAULT;
  /** Small error failure threshold */
  cv_small_nef: number = SMALL_NEF_DEFAULT;

  // --- Counters ---
  /** Total internal steps taken */
  cv_nst: number = 0;
  /** Total RHS evaluations */
  cv_nfe: number = 0;
  /** Convergence failures */
  cv_ncfn: number = 0;
  /** Nonlinear iterations */
  cv_nni: number = 0;
  /** Nonlinear iteration failures */
  cv_nnf: number = 0;
  /** Error test failures */
  cv_netf: number = 0;
  /** Linear solver setups */
  cv_nsetups: number = 0;
  /** t+h==t warnings */
  cv_nhnil: number = 0;

  // --- Step size ratios for order selection ---
  /** Eta for order q-1 */
  cv_etaqm1: number = 0;
  /** Eta for order q */
  cv_etaq: number = 0;
  /** Eta for order q+1 */
  cv_etaqp1: number = 0;

  // --- Linear solver interface ---
  /** Linear solver init function */
  cv_linit: ((mem: CvodeMem) => number) | null = null;
  /** Linear solver setup function */
  cv_lsetup: ((mem: CvodeMem, convfail: number, ypred: Float64Array,
               fpred: Float64Array, jcurPtr: { value: boolean },
               tmp1: Float64Array, tmp2: Float64Array,
               tmp3: Float64Array) => number) | null = null;
  /** Linear solver solve function */
  cv_lsolve: ((mem: CvodeMem, b: Float64Array, weight: Float64Array,
               ycur: Float64Array, fcur: Float64Array) => number) | null = null;
  /** Linear solver free function */
  cv_lfree: ((mem: CvodeMem) => void) | null = null;
  /** Linear solver memory */
  cv_lmem: CvLsMem | null = null;
  /** Max steps between Jacobian setups */
  cv_msbp: number = MSBP_DEFAULT;
  /** Gamma change threshold for Jacobian redo */
  cv_dgmax_lsetup: number = DGMAX_LSETUP_DEFAULT;

  // --- Saved values ---
  /** Order used on last successful step */
  cv_qu: number = 0;
  /** nst at last Jacobian setup */
  cv_nstlp: number = 0;
  /** Initial step size actually used */
  cv_h0u: number = 0;
  /** Step size used on last successful step */
  cv_hu: number = 0;
  /** Saved tq[5] for order selection */
  cv_saved_tq5: number = 0;
  /** Is Jacobian current? */
  cv_jcur: boolean = false;
  /** Tolerance scale factor */
  cv_tolsf: number = 0;
  /** Max order allocated in Nordsieck array */
  cv_qmax_alloc: number = 0;
  /** Index of acor in zn array */
  cv_indx_acor: number = 0;

  // --- Initialization flags ---
  /** Has CvodeMem been allocated/initialized? */
  cv_MallocDone: boolean = false;
  /** Has vector abstol been allocated? */
  cv_VabstolMallocDone: boolean = false;

  // --- BDF stability limit detection ---
  /** Enable stability limit detection? */
  cv_sldeton: boolean = false;
  /** Stability data: ssdat[6][4] (6 rows, 4 cols) */
  cv_ssdat: number[][] = [];
  /** Steps since stability check */
  cv_nscon: number = 0;
  /** Order flag for stability */
  cv_nor: number = 0;

  // --- Rootfinding ---
  /** Root function */
  cv_gfun: CvodeRootFn | null = null;
  /** Number of root functions */
  cv_nrtfn: number = 0;
  /** Root direction info (which roots were found) */
  cv_iroots: Int32Array = new Int32Array(0);
  /** Root direction constraints */
  cv_rootdir: Int32Array = new Int32Array(0);
  /** t at low end of bracket */
  cv_tlo: number = 0;
  /** t at high end of bracket */
  cv_thi: number = 0;
  /** t at root */
  cv_trout: number = 0;
  /** g values at low end */
  cv_glo: Float64Array = new Float64Array(0);
  /** g values at high end */
  cv_ghi: Float64Array = new Float64Array(0);
  /** g values at root */
  cv_grout: Float64Array = new Float64Array(0);
  /** Output time saved for root check */
  cv_toutc: number = 0;
  /** Root time tolerance */
  cv_ttol: number = 0;
  /** Task flag for root check */
  cv_taskc: number = 0;
  /** Root found on previous step? */
  cv_irfnd: number = 0;
  /** Number of g evaluations */
  cv_nge: number = 0;
  /** Active root functions flags */
  cv_gactive: Uint8Array = new Uint8Array(0);
  /** Max steps with g inactive before declaring permanently inactive */
  cv_mxgnull: number = 1;

  // --- NLS ---
  /** Convergence failure type (CV_NO_FAILURES, CV_FAIL_BAD_J, CV_FAIL_OTHER) */
  convfail: number = CV_NO_FAILURES;

  // --- Error ---
  /** Error message string */
  cv_error: string = '';

  // --- Problem size ---
  /** Number of equations */
  cv_N: number = 0;

  // --- Resize flag ---
  /** First step after problem resize */
  first_step_after_resize: boolean = false;
}

// ============================================================================
// Utility functions
// ============================================================================

/**
 * Weighted root-mean-square norm.
 * Returns sqrt(sum((v[i]*w[i])^2) / n) for i = 0..n-1.
 */
export function wrmsNorm(n: number, v: Float64Array, w: Float64Array): number {
  let sum = 0.0;
  for (let i = 0; i < n; i++) {
    const vw = v[i] * w[i];
    sum += vw * vw;
  }
  return Math.sqrt(sum / n);
}

/**
 * Linear combination of two vectors: z[i] = a*x[i] + b*y[i] for i = 0..n-1.
 */
export function vLinearSum(a: number, x: Float64Array, b: number, y: Float64Array,
  z: Float64Array, n: number): void {
  for (let i = 0; i < n; i++)
    z[i] = a * x[i] + b * y[i];
}

/**
 * Scale a vector: z[i] = c * x[i] for i = 0..n-1.
 */
export function vScale(c: number, x: Float64Array, z: Float64Array, n: number): void {
  for (let i = 0; i < n; i++)
    z[i] = c * x[i];
}

/**
 * Set all elements of a vector to a constant: z[i] = c for i = 0..n-1.
 */
export function vConst(c: number, z: Float64Array, n: number): void {
  z.fill(c, 0, n);
}

/**
 * Dense output interpolation from the Nordsieck history array.
 *
 * Computes the k-th derivative of the interpolating polynomial at time t:
 *   dky = sum_{j=k}^{q} c(j,k) * s^(j-k) * h^(-k) * zn[j]
 * where s = (t - tn) / h, c(j,k) = j*(j-1)*...*(j-k+1).
 *
 * @returns 0 on success, -1 on invalid input.
 */
export function cvodeGetDky(
  mem: CvodeMem, t: number, k: number, dky: Float64Array,
): number {
  if (k < 0 || k > mem.cv_q) return -1;

  const N = mem.cv_N;
  const tfuzz = FUZZ_FACTOR * UROUND *
    (Math.abs(mem.cv_tn) + Math.abs(mem.cv_hu));
  const tp = mem.cv_tn - mem.cv_hu - Math.abs(tfuzz);
  const tn1 = mem.cv_tn + Math.abs(tfuzz);
  if ((t - tp) * (t - tn1) > 0) return -1;

  const s = (t - mem.cv_tn) / mem.cv_h;

  dky.fill(0.0, 0, N);
  for (let j = mem.cv_q; j >= k; j--) {
    let c = 1.0;
    for (let i = j; i >= j - k + 1; i--) c *= i;
    for (let i = 0; i < j - k; i++) c *= s;
    const znj = mem.cv_zn[j];
    for (let i = 0; i < N; i++) dky[i] += c * znj[i];
  }

  if (k > 0) {
    const r = Math.pow(mem.cv_h, -k);
    for (let i = 0; i < N; i++) dky[i] *= r;
  }

  return 0;
}
