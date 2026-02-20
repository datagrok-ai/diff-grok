/* LSODA — Livermore Solver for Ordinary Differential Equations with
   Automatic method switching between Adams-Moulton (non-stiff, orders 1–12)
   and BDF (stiff, orders 1–5), using the Nordsieck array representation.

   References:
   [1] Hindmarsh, A.C. and Petzold, L.R., "LSODA, Ordinary Differential
       Equation Solver for Stiff or Non-Stiff System", ODEPACK, 1983.
   [2] Byrne, G.D. and Hindmarsh, A.C., "A Polyalgorithm for the Numerical
       Solution of Ordinary Differential Equations", ACM TOMS, 1(1), 1975.
   [3] Hairer, E. and Wanner, G., "Solving Ordinary Differential Equations
       II: Stiff and Differential-Algebraic Problems", Springer, 1996.
   [4] Radhakrishnan, K. and Hindmarsh, A.C., "Description and Use of LSODE,
       the Livermore Solver for Ordinary Differential Equations", 1993.
*/

import {ODEs, max, abs, SAFETY, PSHRNK, PSGROW, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, EPS, jacobian, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {luDecomp, luSolve, solve1d2d} from './lin-alg-tools';

// ================================================================
//  MAXIMUM ORDERS
// ================================================================

const ADAMS_MAX_ORDER = 12;
const BDF_MAX_ORDER = 5;

// ================================================================
//  NEWTON ITERATION
// ================================================================

const MAX_NEWTON_ITER = 7;
const NEWTON_CONV_TOL = 0.33;

// ================================================================
//  STIFFNESS DETECTION
// ================================================================

const STIFF_DETECT_THRESHOLD = 5;
const NONSTIFF_TEST_INTERVAL = 100;

// ================================================================
//  STEP SIZE ADAPTATION
// ================================================================

const MAX_STEP_GROW = 5.0;
const MIN_STEP_SHRINK = 0.2;

// ================================================================
//  MODE SWITCH SAFEGUARDS
// ================================================================

const MAX_SWITCH_COOLDOWN = 3;
const MAX_SWITCH_GROW = 2.0;

// ================================================================
//  JACOBIAN REFRESH
// ================================================================

const MAX_JACOBIAN_AGE = 20;

// ================================================================
//  REPEATED FAILURE THRESHOLDS
// ================================================================

const MAX_CONSECUTIVE_FAILS = 30;
const FORCE_JAC_RECOMPUTE_FAILS = 3;
const ORDER_DROP_FAILS = 7;

// ================================================================
//  RKF45 BUTCHER TABLEAU (for bootstrap & trial steps)
// ================================================================

const RK_C2 = 1.0 / 4.0;
const RK_C3 = 3.0 / 8.0;
const RK_C4 = 12.0 / 13.0;
const RK_C6 = 1.0 / 2.0;

const RK_A21 = 1.0 / 4.0;
const RK_A31 = 3.0 / 32.0;
const RK_A32 = 9.0 / 32.0;
const RK_A41 = 1932.0 / 2197.0;
const RK_A42 = -7200.0 / 2197.0;
const RK_A43 = 7296.0 / 2197.0;
const RK_A51 = 439.0 / 216.0;
const RK_A52 = -8.0;
const RK_A53 = 3680.0 / 513.0;
const RK_A54 = -845.0 / 4104.0;
const RK_A61 = -8.0 / 27.0;
const RK_A62 = 2.0;
const RK_A63 = -3544.0 / 2565.0;
const RK_A64 = 1859.0 / 4104.0;
const RK_A65 = -11.0 / 40.0;

const RK_B1 = 25.0 / 216.0;
const RK_B3 = 1408.0 / 2565.0;
const RK_B4 = 2197.0 / 4104.0;
const RK_B5 = -1.0 / 5.0;

const RK_E1 = 1.0 / 360.0;
const RK_E3 = -128.0 / 4275.0;
const RK_E4 = -2197.0 / 75240.0;
const RK_E5 = 1.0 / 50.0;
const RK_E6 = 2.0 / 55.0;

// ================================================================
//  BINOMIAL COEFFICIENTS (precomputed for Pascal shift)
// ================================================================

const MAX_ROWS = ADAMS_MAX_ORDER + 1; // 13

/** Precomputed binomial coefficients C(n, k) for n = 0..MAX_ROWS-1 */
const BINOM: number[][] = [];
for (let n = 0; n < MAX_ROWS; n++) {
  BINOM[n] = new Array(n + 1);
  BINOM[n][0] = 1;
  BINOM[n][n] = 1;
  for (let k = 1; k < n; k++)
    BINOM[n][k] = BINOM[n - 1][k - 1] + BINOM[n - 1][k];
}

// ================================================================
//  COEFFICIENT COMPUTATION (from LSODE cfode subroutine)
// ================================================================

interface MethodCoeffs {
  /** Correction vectors l[order-1][0..order], one per order */
  l: number[][];
  /** Error coefficients for each order */
  errCoeff: number[];
}

/** Compute Adams-Moulton correction vectors and error coefficients
 *  for orders 1 through ADAMS_MAX_ORDER.
 *  Translated from LSODE cfode (meth=1).
 */
function computeAdamsCoeffs(): MethodCoeffs {
  const l: number[][] = [];
  const errCoeff: number[] = [];

  // Order 1: implicit Euler equivalent
  l.push([1, 1]);
  errCoeff.push(1.0 / 3.0);

  const pc = new Array(ADAMS_MAX_ORDER + 1).fill(0);
  pc[0] = 1;
  let rqfac = 1.0;

  for (let nq = 2; nq <= ADAMS_MAX_ORDER; nq++) {
    const rq1fac = rqfac;
    rqfac /= nq;
    const nqm1 = nq - 1;

    // Multiply polynomial p(x) by (x + nq - 1):
    //   p(x) accumulates (x+1)(x+2)...(x+nq-1)
    pc[nq - 1] = 0;
    for (let ib = 1; ib <= nqm1; ib++) {
      const idx = nq - ib;
      pc[idx] = pc[idx - 1] + nqm1 * pc[idx];
    }
    pc[0] = nqm1 * pc[0];

    // Integrate p(x) from -1 to 0
    let pint = pc[0];
    let tsign = 1.0;
    for (let i = 1; i < nq; i++) {
      tsign = -tsign;
      pint += tsign * pc[i] / (i + 1);
    }

    // Build the l vector
    const lVec = new Array(nq + 1);
    lVec[0] = pint * rqfac;
    lVec[1] = 1.0;
    for (let i = 2; i <= nq; i++)
      lVec[i] = rq1fac * pc[i - 1] / i;

    l.push(lVec);
    errCoeff.push(lVec[0] / (nq + 2));
  }

  return {l, errCoeff};
}

/** Compute BDF correction vectors and error coefficients
 *  for orders 1 through BDF_MAX_ORDER.
 *  Translated from LSODE cfode (meth=2).
 */
function computeBDFCoeffs(): MethodCoeffs {
  const l: number[][] = [];
  const errCoeff: number[] = [];

  const pc = new Array(BDF_MAX_ORDER + 2).fill(0);
  pc[0] = 1;

  for (let nq = 1; nq <= BDF_MAX_ORDER; nq++) {
    // Multiply polynomial by (x + nq):
    //   p(x) accumulates (x+1)(x+2)...(x+nq)
    pc[nq] = 0;
    for (let ib = 1; ib <= nq; ib++) {
      const idx = nq + 1 - ib;
      pc[idx] = pc[idx - 1] + ib * pc[idx];
    }
    pc[0] = nq * pc[0];

    // Normalize: l[i] = pc[i] / pc[1] so that l[1] = 1
    const lVec = new Array(nq + 1);
    for (let i = 0; i <= nq; i++)
      lVec[i] = pc[i] / pc[1];

    l.push(lVec);
    errCoeff.push(lVec[0] / (nq + 2));
  }

  return {l, errCoeff};
}

/** Precomputed Adams coefficients (indexed by order - 1) */
const ADAMS = computeAdamsCoeffs();

/** Precomputed BDF coefficients (indexed by order - 1) */
const BDF_COEFFS = computeBDFCoeffs();

// ================================================================
//  NORDSIECK ARRAY OPERATIONS
// ================================================================

/** Initialize Nordsieck array at order 1: z = [y0, h*f0, 0, 0, ...] */
function nordsieckInit(y0: Float64Array, f0: Float64Array, h: number,
  z: Float64Array, dim: number): void {
  z.fill(0);
  for (let j = 0; j < dim; j++) {
    z[j] = y0[j];
    z[dim + j] = h * f0[j];
  }
}

/** Predict: apply Pascal shift z_pred = P * z.
 *  z_pred[i] = sum_{k=i}^{q} C(k,i) * z[k] for each row i.
 */
function nordsieckPredict(z: Float64Array, zPred: Float64Array,
  q: number, dim: number): void {
  for (let i = 0; i <= q; i++) {
    const iOff = i * dim;
    // Start with the diagonal term C(i,i) * z[i] = 1 * z[i]
    for (let j = 0; j < dim; j++)
      zPred[iOff + j] = z[iOff + j];
    // Add off-diagonal terms
    for (let k = i + 1; k <= q; k++) {
      const c = BINOM[k][i];
      const kOff = k * dim;
      for (let j = 0; j < dim; j++)
        zPred[iOff + j] += c * z[kOff + j];
    }
  }
}

/** Correct: z_new[i] = z_pred[i] + l[i] * delta for i = 0..q */
function nordsieckCorrect(zPred: Float64Array, zNew: Float64Array,
  l: number[], delta: Float64Array, q: number, dim: number): void {
  for (let i = 0; i <= q; i++) {
    const iOff = i * dim;
    const li = l[i];
    for (let j = 0; j < dim; j++)
      zNew[iOff + j] = zPred[iOff + j] + li * delta[j];
  }
}

/** Rescale Nordsieck array after step size change: z[i] *= eta^i */
function nordsieckRescale(z: Float64Array, eta: number,
  q: number, dim: number): void {
  let etaPow = 1.0;
  for (let i = 1; i <= q; i++) {
    etaPow *= eta;
    const iOff = i * dim;
    for (let j = 0; j < dim; j++)
      z[iOff + j] *= etaPow;
  }
}

// ================================================================
//  MAIN SOLVER
// ================================================================

/** Solve initial value problem using the LSODA method with automatic
 *  stiffness detection and variable-order Adams/BDF switching.
 *  @param odes initial value problem for ordinary differential equations
 *  @param callback computations control callback
 *  @returns solution of the problem
 */
export function lsodaWeb(odes: ODEs, callback?: Callback): Float64Array[] {
  const f = odes.func;

  // --- Operating variables ---
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  const hDataframe = odes.arg.step;
  const hAdamsMax = hDataframe * 10;
  const tolerance = odes.tolerance;
  let h = hDataframe;

  const rowCount = Math.trunc((t1 - t0) / h) + 1;
  const dim = odes.initial.length;
  const dimSquared = dim * dim;

  // --- Result arrays ---
  const tArr = new Float64Array(rowCount);
  const yArrs = Array<Float64Array>(dim);
  for (let i = 0; i < dim; i++)
    yArrs[i] = new Float64Array(rowCount);

  // --- Loop control ---
  let timeDataframe = t0 + hDataframe;
  let t = t0;
  let tPrev = t0;
  let hNext = 0.0;
  let flag = true;
  let index = 1;
  let errmax = 0.0;
  let hTemp = 0.0;
  let tNew = 0.0;

  // ================================================================
  //  BUFFER PREALLOCATION
  // ================================================================

  // --- Common buffers ---
  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  // --- Nordsieck buffers ---
  const nordsieck = new Float64Array(MAX_ROWS * dim);
  const nordsieckPred = new Float64Array(MAX_ROWS * dim);
  const delta = new Float64Array(dim);
  const deltaPrev = new Float64Array(dim);

  // --- BDF / Newton buffers ---
  const Ident = new Float64Array(dimSquared);
  for (let i = 0; i < dim; i++)
    Ident[i + i * dim] = 1.0;

  const W = new Float64Array(dimSquared);
  const Jac = new Float64Array(dimSquared);
  const L = new Float64Array(dimSquared);
  const U = new Float64Array(dimSquared);
  const luBuf = new Float64Array(dim);
  const f0Buf = new Float64Array(dim);
  const f1Buf = new Float64Array(dim);
  const G = new Float64Array(dim);
  const yNewton = new Float64Array(dim);
  const fNewton = new Float64Array(dim);
  const toUseLU = dim > 2;

  // --- RKF45 buffers (bootstrap + trial steps) ---
  const k1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const k4 = new Float64Array(dim);
  const k5 = new Float64Array(dim);
  const k6 = new Float64Array(dim);

  // ================================================================
  //  STATE VARIABLES
  // ================================================================

  let isStiff = false;
  let q = 1; // current method order
  let stepsAtCurrentOrder = 0;
  let prevDeltaValid = false;
  let wasRejected = false;
  let consecutiveFails = 0;

  // --- Adams-specific ---
  let consecutiveAdamsRejects = 0;

  // --- BDF-specific ---
  let acceptedBdfSteps = 0;
  let jacobianValid = false;
  let jacobianAge = 0;
  let hAtJacobian = 0.0;
  let lastNewtonTheta = 0.0;

  // --- Mode switching ---
  let stepsSinceSwitch = MAX_SWITCH_COOLDOWN; // start with no cooldown

  // --- Step size tracking (for Nordsieck rescaling) ---
  let hNordsieck = h; // the h that the Nordsieck array is currently scaled to

  // ================================================================
  //  INITIALIZATION
  // ================================================================

  // Compute initial f-value
  f(t, y, dydt);

  // Estimate initial step size (Hairer-Wanner II.4 heuristic)
  {
    let d0 = 0;
    let d1 = 0;
    for (let i = 0; i < dim; i++) {
      const sc = max(abs(y[i]), tolerance);
      d0 = max(d0, abs(y[i]) / sc);
      d1 = max(d1, abs(dydt[i]) / sc);
    }

    let h0: number;
    if (d0 < 1e-5 || d1 < 1e-5)
      h0 = 1e-6;
    else
      h0 = 0.01 * d0 / d1;

    h0 = Math.min(h0, t1 - t0);

    // Explicit Euler probe step to estimate second derivative norm
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h0 * dydt[i];
    f(t0 + h0, yTemp, fNewton);

    let d2 = 0;
    for (let i = 0; i < dim; i++) {
      const sc = max(abs(y[i]), tolerance);
      d2 = max(d2, abs(fNewton[i] - dydt[i]) / sc);
    }
    d2 /= h0;

    let h1: number;
    if (max(d1, d2) <= 1e-15)
      h1 = max(1e-6, h0 * 1e-3);
    else
      h1 = (0.01 / max(d1, d2)) ** 0.5;

    h = Math.min(100 * h0, h1, hDataframe);
  }

  // Build Nordsieck array at order 1
  nordsieckInit(y, dydt, h, nordsieck, dim);
  hNordsieck = h;

  // 1. SOLUTION AT THE POINT t0
  tArr[0] = t0;
  for (let i = 0; i < dim; i++)
    yArrs[i][0] = y[i];

  // ================================================================
  //  MAIN INTEGRATION LOOP
  // ================================================================

  while (flag) {
    // Compute derivative and scale vector
    f(t, y, dydt);

    if (callback)
      callback.onIterationStart();

    for (let i = 0; i < dim; i++)
      yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;

    // Check end point
    if (t + h > t1) {
      h = t1 - t;
      flag = false;
    }

    // Rescale Nordsieck array if h differs from what it was scaled to
    if (abs(h - hNordsieck) > 1e-15 * abs(hNordsieck) && hNordsieck !== 0) {
      const eta = h / hNordsieck;
      nordsieckRescale(nordsieck, eta, q, dim);
      hNordsieck = h;
    }

    // ============================================================
    //  ADAPTIVE STEP LOOP
    // ============================================================

    while (true) {
      if (!isStiff) {
        // ======================================================
        //  ADAMS MODE (non-stiff) — PEC scheme
        // ======================================================

        if (q === 1 && stepsAtCurrentOrder === 0) {
          // --- BOOTSTRAP: use RKF45 for the first step to build history ---
          rkf45Step();
        } else {
          // --- NORDSIECK ADAMS PEC ---
          adamsPECStep();
        }
      } else {
        // ======================================================
        //  BDF MODE (stiff) — Newton iteration
        // ======================================================

        bdfNewtonStep();
      }

      // Check for step acceptance (errmax set by the step functions)
      if (errmax > 1.0) {
        // --- STEP REJECTED ---
        const pshrnk = -1.0 / (q + 1);
        hTemp = SAFETY * h * errmax ** pshrnk;
        const hOld = h;
        h = max(hTemp, REDUCE_COEF * h);
        tNew = t + h;
        if (tNew === t)
          throw new Error(ERROR_MSG.LSODA_FAILS);

        wasRejected = true;
        consecutiveFails++;

        // Rescale Nordsieck array for the new (smaller) step size
        const eta = h / hOld;
        nordsieckRescale(nordsieck, eta, q, dim);
        hNordsieck = h;

        // Stiffness detection in Adams mode
        if (!isStiff) {
          consecutiveAdamsRejects++;
          if (consecutiveAdamsRejects >= STIFF_DETECT_THRESHOLD)
            switchToStiff();
        } else {
          // BDF: force Jacobian recomputation on repeated failures
          if (consecutiveFails >= FORCE_JAC_RECOMPUTE_FAILS)
            jacobianValid = false;
          if (consecutiveFails >= ORDER_DROP_FAILS && q > 1) {
            q = 1;
            stepsAtCurrentOrder = 0;
            prevDeltaValid = false;
            // Reinitialize Nordsieck at order 1 since we dropped order
            f(t, y, dydt);
            nordsieckInit(y, dydt, h, nordsieck, dim);
            hNordsieck = h;
          }
        }

        if (consecutiveFails >= MAX_CONSECUTIVE_FAILS)
          throw new Error(ERROR_MSG.LSODA_FAILS);

        continue;
      }

      // --- STEP ACCEPTED ---

      // Update Nordsieck array with the converged correction
      const lVec = isStiff ? BDF_COEFFS.l[q - 1] : ADAMS.l[q - 1];
      nordsieckCorrect(nordsieckPred, nordsieck, lVec, delta, q, dim);

      // Extract y from updated Nordsieck row 0
      for (let j = 0; j < dim; j++)
        y[j] = nordsieck[j];

      t = t + h;
      consecutiveFails = 0;

      if (!isStiff)
        consecutiveAdamsRejects = 0;
      else {
        acceptedBdfSteps++;
        jacobianAge++;
      }

      stepsAtCurrentOrder++;

      // --- ORDER AND STEP SIZE SELECTION ---
      selectOrderAndStepSize();

      // --- POST-ACCEPT GUARDS ---
      applyPostAcceptGuards();

      // --- STIFFNESS RE-EVALUATION (BDF → Adams) ---
      // Only consider switching when h has grown to at least 1% of output step,
      // preventing false positives at microscopic step sizes.
      if (isStiff && acceptedBdfSteps >= NONSTIFF_TEST_INTERVAL && h > hDataframe * 0.01) {
        acceptedBdfSteps = 0;
        const _trial = trialExplicitStep();
        if (_trial)
          switchToAdams();
      }

      break;
    } // while (true) — adaptive step

    // ============================================================
    //  INTERPOLATION TO OUTPUT GRID
    // ============================================================

    while (timeDataframe < t) {
      const cLeft = (t - timeDataframe) / (t - tPrev);
      const cRight = 1.0 - cLeft;

      tArr[index] = timeDataframe;

      for (let j = 0; j < dim; j++)
        yArrs[j][index] = cRight * y[j] + cLeft * yPrev[j];

      timeDataframe += hDataframe;
      ++index;
    }

    h = hNext;
    tPrev = t;

    for (let i = 0; i < dim; i++)
      yPrev[i] = y[i];
  } // while (flag) — main loop

  // --- Final callback ---
  if (callback)
    callback.onComputationsCompleted();

  // 3. SOLUTION AT THE POINT t1
  tArr[rowCount - 1] = t1;
  for (let i = 0; i < dim; i++)
    yArrs[i][rowCount - 1] = y[i];

  // 4. PREPARE OUTPUT
  const solution = Array<Float64Array>(dim + 1);
  solution[0] = tArr;
  for (let i = 0; i < dim; i++)
    solution[i + 1] = yArrs[i];

  return solution;

  // ================================================================
  //  CLOSURE HELPER FUNCTIONS
  // ================================================================

  /** Perform a single RKF45 step (bootstrap for Adams mode). Sets errmax. */
  function rkf45Step(): void {
    // Stage 1
    f(t, y, k1);

    // Stage 2
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * RK_A21 * k1[i];
    f(t + RK_C2 * h, yTemp, k2);

    // Stage 3
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A31 * k1[i] + RK_A32 * k2[i]);
    f(t + RK_C3 * h, yTemp, k3);

    // Stage 4
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A41 * k1[i] + RK_A42 * k2[i] + RK_A43 * k3[i]);
    f(t + RK_C4 * h, yTemp, k4);

    // Stage 5
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A51 * k1[i] + RK_A52 * k2[i] + RK_A53 * k3[i] + RK_A54 * k4[i]);
    f(t + h, yTemp, k5);

    // Stage 6
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A61 * k1[i] + RK_A62 * k2[i] + RK_A63 * k3[i] + RK_A64 * k4[i] + RK_A65 * k5[i]);
    f(t + RK_C6 * h, yTemp, k6);

    // 4th-order solution
    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_B1 * k1[i] + RK_B3 * k3[i] + RK_B4 * k4[i] + RK_B5 * k5[i]);

    // Error estimate
    for (let i = 0; i < dim; i++)
      yErr[i] = h * (RK_E1 * k1[i] + RK_E3 * k3[i] + RK_E4 * k4[i] + RK_E5 * k5[i] + RK_E6 * k6[i]);

    // Compute error norm
    errmax = 0;
    for (let i = 0; i < dim; i++)
      errmax = max(errmax, abs(yErr[i] / yScale[i]));
    errmax /= tolerance;

    if (errmax <= 1.0) {
      // On acceptance, build fresh Nordsieck array at the new point
      f(t + h, yTemp, fNewton);

      // Set Nordsieck: row 0 = y, row 1 = h*f, rest = 0
      nordsieck.fill(0);
      for (let j = 0; j < dim; j++) {
        nordsieck[j] = yTemp[j];
        nordsieck[dim + j] = h * fNewton[j];
      }

      // Copy to nordsieckPred so the main accept logic (nordsieckCorrect) is a no-op
      for (let j = 0; j < MAX_ROWS * dim; j++)
        nordsieckPred[j] = nordsieck[j];

      // delta = 0 (correction already applied via RKF45)
      delta.fill(0);
    }
  }

  /** Perform Adams PEC step using Nordsieck array. Sets errmax. */
  function adamsPECStep(): void {
    const lVec = ADAMS.l[q - 1];
    const errC = ADAMS.errCoeff[q - 1];

    // Predict
    nordsieckPredict(nordsieck, nordsieckPred, q, dim);

    // Evaluate f at predicted y
    // yTemp = predicted y (row 0 of nordsieckPred)
    for (let j = 0; j < dim; j++)
      yTemp[j] = nordsieckPred[j];

    f(t + h, yTemp, fNewton);

    // Compute derivative mismatch delta = h*f - predicted h*y'
    for (let j = 0; j < dim; j++)
      delta[j] = h * fNewton[j] - nordsieckPred[dim + j];

    // Corrected y
    for (let j = 0; j < dim; j++)
      yTemp[j] = nordsieckPred[j] + lVec[0] * delta[j];

    // Second evaluation at corrected point (PECE)
    f(t + h, yTemp, fNewton);

    // Recompute delta at corrected point
    for (let j = 0; j < dim; j++)
      delta[j] = h * fNewton[j] - nordsieckPred[dim + j];

    // Error estimate
    for (let j = 0; j < dim; j++)
      yErr[j] = errC * delta[j];

    errmax = 0;
    for (let j = 0; j < dim; j++)
      errmax = max(errmax, abs(yErr[j] / yScale[j]));
    errmax /= tolerance;
  }

  /** Perform BDF step with Newton iteration. Sets errmax. */
  function bdfNewtonStep(): void {
    const lVec = BDF_COEFFS.l[q - 1];
    const errC = BDF_COEFFS.errCoeff[q - 1];
    const gamma = lVec[0]; // l[0] = BDF gamma

    // --- Jacobian update check ---
    const needJacobian = !jacobianValid ||
      jacobianAge > MAX_JACOBIAN_AGE ||
      abs(h - hAtJacobian) > 0.5 * abs(hAtJacobian);

    if (needJacobian) {
      jacobian(t, y, f, EPS, f0Buf, f1Buf, Jac);
      jacobianValid = true;
      jacobianAge = 0;
      hAtJacobian = h;
    }

    // --- Form iteration matrix W = I - gamma*h*J ---
    const ghCoeff = gamma * h;
    for (let i = 0; i < dimSquared; i++)
      W[i] = Ident[i] - ghCoeff * Jac[i];

    // LU decomposition
    if (toUseLU)
      luDecomp(W, L, U, dim);

    // --- Predict ---
    nordsieckPredict(nordsieck, nordsieckPred, q, dim);

    // Initial Newton guess: predicted y
    for (let j = 0; j < dim; j++)
      yNewton[j] = nordsieckPred[j];

    // --- Newton iteration ---
    let newtonConverged = false;
    lastNewtonTheta = 0;
    let prevDeltaNorm = 0;

    for (let iter = 0; iter < MAX_NEWTON_ITER; iter++) {
      // Evaluate f at current iterate
      f(t + h, yNewton, fNewton);

      // Derivative mismatch
      for (let j = 0; j < dim; j++)
        delta[j] = h * fNewton[j] - nordsieckPred[dim + j];

      // Accumulated correction e = yNewton - yPredicted
      // Residual: G = -(e - l[0] * delta)
      for (let j = 0; j < dim; j++)
        G[j] = -(yNewton[j] - nordsieckPred[j] - gamma * delta[j]);

      // Solve W * deltaE = G
      if (toUseLU)
        luSolve(L, U, G, luBuf, yTemp, dim);
      else
        solve1d2d(W, G, yTemp);

      // Update Newton iterate: yNewton += deltaE
      for (let j = 0; j < dim; j++)
        yNewton[j] += yTemp[j];

      // Convergence check (max norm of deltaE, using tolerance-based weight)
      let deltaNorm = 0;
      for (let j = 0; j < dim; j++) {
        const ewt = tolerance * abs(yNewton[j]) + tolerance;
        deltaNorm = max(deltaNorm, abs(yTemp[j]) / ewt);
      }

      // Contraction factor
      if (iter > 0 && prevDeltaNorm > 0)
        lastNewtonTheta = deltaNorm / prevDeltaNorm;

      prevDeltaNorm = deltaNorm;

      if (deltaNorm < NEWTON_CONV_TOL) {
        newtonConverged = true;
        break;
      }

      // Early divergence detection
      if (iter > 0 && lastNewtonTheta > NEWTON_CONV_TOL)
        break;
    }

    if (!newtonConverged) {
      // Signal rejection — use errmax=4 for gentle step shrink (factor ~0.5)
      // instead of 1e6 which causes catastrophic h collapse
      errmax = 4.0;
      jacobianValid = false;
      return;
    }

    // Compute final delta at converged point
    f(t + h, yNewton, fNewton);
    for (let j = 0; j < dim; j++)
      delta[j] = h * fNewton[j] - nordsieckPred[dim + j];

    // Error estimate (using tolerance-based weight for BDF — ewt already includes tolerance)
    errmax = 0;
    for (let j = 0; j < dim; j++) {
      const ewt = tolerance * abs(yNewton[j]) + tolerance;
      yErr[j] = errC * delta[j];
      errmax = max(errmax, abs(yErr[j]) / ewt);
    }
  }

  /** Select new order and step size after an accepted step. Sets hNext. */
  function selectOrderAndStepSize(): void {
    const maxOrder = isStiff ? BDF_MAX_ORDER : ADAMS_MAX_ORDER;
    const coeffs = isStiff ? BDF_COEFFS : ADAMS;
    const psgrow = -1.0 / (q + 1);

    // --- Eta for current order q ---
    let etaSame: number;
    if (errmax > ERR_CONTR)
      etaSame = Math.min(SAFETY * errmax ** psgrow, MAX_STEP_GROW);
    else
      etaSame = MAX_STEP_GROW;

    let qNew = q;
    let etaBest = etaSame;

    // --- Eta for order q-1 ---
    if (q > 1) {
      // Error estimate from Nordsieck component z[q]
      let errDown = 0;
      for (let j = 0; j < dim; j++) {
        const val = nordsieck[q * dim + j];
        const w = isStiff ? (tolerance * abs(y[j]) + tolerance) : yScale[j];
        errDown = max(errDown, abs(val) / w);
      }
      if (!isStiff) errDown /= tolerance;

      let etaDown: number;
      if (errDown > TINY) {
        const p = -1.0 / q;
        etaDown = Math.min(SAFETY * errDown ** p, MAX_STEP_GROW);
      } else
        etaDown = MAX_STEP_GROW;


      if (etaDown > etaBest) {
        qNew = q - 1;
        etaBest = etaDown;
      }
    }

    // --- Eta for order q+1 ---
    if (q < maxOrder && stepsAtCurrentOrder >= q + 1 && !wasRejected && prevDeltaValid) {
      // Error estimate from difference of consecutive corrections
      const errCoeffUp = (q < maxOrder) ? coeffs.errCoeff[q] : 1.0; // errCoeff for order q+1
      let errUp = 0;
      for (let j = 0; j < dim; j++) {
        const val = errCoeffUp * (delta[j] - deltaPrev[j]);
        const w = isStiff ? (tolerance * abs(y[j]) + tolerance) : yScale[j];
        errUp = max(errUp, abs(val) / w);
      }
      if (!isStiff) errUp /= tolerance;

      let etaUp: number;
      if (errUp > TINY) {
        const p = -1.0 / (q + 2);
        etaUp = Math.min(SAFETY * errUp ** p, MAX_STEP_GROW);
      } else
        etaUp = MAX_STEP_GROW;


      if (etaUp > etaBest) {
        qNew = q + 1;
        etaBest = etaUp;
      }
    }

    // Clamp eta
    etaBest = Math.min(etaBest, MAX_STEP_GROW);
    etaBest = Math.max(etaBest, MIN_STEP_SHRINK);

    // Save current delta for next step's order-up estimate
    for (let j = 0; j < dim; j++)
      deltaPrev[j] = delta[j];
    prevDeltaValid = true;

    // --- Apply order change ---
    if (qNew !== q) {
      if (qNew > q) {
        // Increasing order: initialize the new Nordsieck row
        // z[qNew][j] = errCoeff_new * delta[j] (rough estimate of higher derivative)
        const newOff = qNew * dim;
        const eC = coeffs.errCoeff[qNew - 1];
        for (let j = 0; j < dim; j++)
          nordsieck[newOff + j] = eC * delta[j];
      }
      // For order decrease: the extra row is simply not used

      q = qNew;
      stepsAtCurrentOrder = 0;
      prevDeltaValid = false;
    }

    // --- Apply step size change ---
    hNext = h * etaBest;
    if (!isStiff && hNext > hAdamsMax)
      hNext = hAdamsMax;

    // Rescale Nordsieck array for new step size
    const finalEta = hNext / h;
    if (abs(finalEta - 1.0) > 1e-12) {
      nordsieckRescale(nordsieck, finalEta, q, dim);
      hNordsieck = hNext;
    }
  }

  /** Post-acceptance guards on step growth. */
  function applyPostAcceptGuards(): void {
    // After rejection, don't grow
    if (wasRejected) {
      hNext = Math.min(hNext, h);
      wasRejected = false;
    }

    // Post-switch cooldown
    if (stepsSinceSwitch < MAX_SWITCH_COOLDOWN)
      hNext = Math.min(hNext, MAX_SWITCH_GROW * h);
    stepsSinceSwitch++;

    // Slow Newton convergence: invalidate Jacobian rather than shrinking h
    if (isStiff && lastNewtonTheta > 0.5)
      jacobianValid = false;
  }

  /** Switch to stiff (BDF) mode. */
  function switchToStiff(): void {
    isStiff = true;
    q = 1;
    stepsAtCurrentOrder = 0;
    prevDeltaValid = false;
    jacobianValid = false;
    acceptedBdfSteps = 0;
    consecutiveAdamsRejects = 0;
    consecutiveFails = 0;
    f(t, y, dydt);
    nordsieckInit(y, dydt, h, nordsieck, dim);
    hNordsieck = h;
    stepsSinceSwitch = 0;
  }

  /** Switch to non-stiff (Adams) mode. */
  function switchToAdams(): void {
    isStiff = false;
    q = 1;
    stepsAtCurrentOrder = 0;
    prevDeltaValid = false;
    consecutiveAdamsRejects = 0;
    consecutiveFails = 0;
    f(t, y, dydt);
    nordsieckInit(y, dydt, h, nordsieck, dim);
    hNordsieck = h;
    stepsSinceSwitch = 0;
  }

  /** Trial RKF45 step to check if problem is no longer stiff.
   *  Does NOT advance the solution.
   *  @returns true if the explicit step succeeds (non-stiff).
   */
  function trialExplicitStep(): boolean {
    f(t, y, k1);

    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * RK_A21 * k1[i];
    f(t + RK_C2 * h, yTemp, k2);

    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A31 * k1[i] + RK_A32 * k2[i]);
    f(t + RK_C3 * h, yTemp, k3);

    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A41 * k1[i] + RK_A42 * k2[i] + RK_A43 * k3[i]);
    f(t + RK_C4 * h, yTemp, k4);

    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A51 * k1[i] + RK_A52 * k2[i] + RK_A53 * k3[i] + RK_A54 * k4[i]);
    f(t + h, yTemp, k5);

    for (let i = 0; i < dim; i++)
      yTemp[i] = y[i] + h * (RK_A61 * k1[i] + RK_A62 * k2[i] + RK_A63 * k3[i] + RK_A64 * k4[i] + RK_A65 * k5[i]);
    f(t + RK_C6 * h, yTemp, k6);

    // Error estimate
    for (let i = 0; i < dim; i++)
      yErr[i] = h * (RK_E1 * k1[i] + RK_E3 * k3[i] + RK_E4 * k4[i] + RK_E5 * k5[i] + RK_E6 * k6[i]);

    let trialErr = 0;
    for (let i = 0; i < dim; i++)
      trialErr = max(trialErr, abs(yErr[i] / yScale[i]));
    trialErr /= tolerance;

    return trialErr <= 1.0;
  }
} // lsodaWeb

// ================================================================
//  FIXED-ORDER LSODA (Adams-4 + BDF2)
// ================================================================

// Adams-Bashforth 4-step coefficients (pre-divided by 24)
const AB0 = 55.0 / 24.0;
const AB1 = -59.0 / 24.0;
const AB2 = 37.0 / 24.0;
const AB3 = -9.0 / 24.0;

// Adams-Moulton 3-step coefficients (pre-divided by 24)
const AM_NEW = 9.0 / 24.0;
const AM0 = 19.0 / 24.0;
const AM1 = -5.0 / 24.0;
const AM2 = 1.0 / 24.0;

// Milne error estimation coefficient for the AB4-AM3 pair
const MILNE = 19.0 / 270.0;

// Number of uniformly-spaced history points required for AB4
const ADAMS_HIST_LEN = 4;

// BDF2 coefficients: y_{n+1} = (4/3)*y_n - (1/3)*y_{n-1} + (2/3)*h*f(t_{n+1}, y_{n+1})
const BDF2_A0 = 4.0 / 3.0;
const BDF2_A1 = -1.0 / 3.0;
const BDF2_GAMMA = 2.0 / 3.0;

// BDF adaptive step exponents (order 2)
const BDF_PSHRNK = -0.5;
const BDF_PSGROW = -1.0 / 3.0;

/** Solve initial value problem using the LSODA method with automatic
 *  stiffness detection. Uses fixed-order Adams-Bashforth-Moulton 4 (non-stiff)
 *  and BDF2 (stiff) with automatic switching.
 *  @param odes initial value problem for ordinary differential equations
 *  @param callback computations control callback
 *  @returns solution of the problem
 */
export function lsoda(odes: ODEs, callback?: Callback): Float64Array[] {
  const f = odes.func;

  // operating variables
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  let h = odes.arg.step;
  const hDataframe = h;
  const hMax = hDataframe * 10;
  const tolerance = odes.tolerance;

  const rowCount = Math.trunc((t1 - t0) / h) + 1;
  const dim = odes.initial.length;
  const dimSquared = dim * dim;

  // result arrays
  const tArr = new Float64Array(rowCount);
  const yArrs = Array<Float64Array>(dim);
  for (let i = 0; i < dim; ++i)
    yArrs[i] = new Float64Array(rowCount);

  // method routine
  let timeDataframe = t0 + hDataframe;
  let t = t0;
  let tPrev = t0;
  let hNext = 0.0;
  let flag = true;
  let index = 1;
  let errmax = 0;
  let hTemp = 0;
  let tNew = 0;

  // --- Common buffers ---
  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  // --- Adams buffers ---
  const fBuf: Float64Array[] = [];
  for (let i = 0; i < ADAMS_HIST_LEN; ++i)
    fBuf[i] = new Float64Array(dim);
  const yCorr = new Float64Array(dim);
  const fPred = new Float64Array(dim);

  // --- RKF45 stage buffers (bootstrap + trial steps) ---
  const k1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const k4 = new Float64Array(dim);
  const k5 = new Float64Array(dim);
  const k6 = new Float64Array(dim);

  // --- BDF buffers ---
  const yBdfPrev = new Float64Array(dim); // y_{n-1} for BDF2
  const yNewton = new Float64Array(dim);
  const fNewton = new Float64Array(dim);
  const delta = new Float64Array(dim);
  const G = new Float64Array(dim);
  const yBDF1 = new Float64Array(dim);

  // --- Linear algebra buffers ---
  const Ident = new Float64Array(dimSquared);
  for (let i = 0; i < dim; ++i)
    Ident[i + i * dim] = 1.0;

  const W = new Float64Array(dimSquared);
  const Jac = new Float64Array(dimSquared);
  const L = new Float64Array(dimSquared);
  const U = new Float64Array(dimSquared);
  const luBuf = new Float64Array(dim);
  const f0Buf = new Float64Array(dim);
  const f1Buf = new Float64Array(dim);
  const toUseLU = dim > 2;

  // --- State variables ---
  let isStiff = false;
  let adamsHistCount = 0;
  let consecutiveAdamsRejects = 0;
  let acceptedBdfSteps = 0;
  let bdfHistCount = 0;
  let jacobianValid = false;
  let jacobianAge = 0;
  let lastAdamsH = 0;
  let lastBdfH = 0;

  // store initial f-value
  f(t, y, fBuf[0]);
  adamsHistCount = 1;

  // 1. SOLUTION AT THE POINT t0
  tArr[0] = t0;
  for (let i = 0; i < dim; ++i)
    yArrs[i][0] = y[i];

  // 2. COMPUTE NUMERICAL SOLUTION
  while (flag) {
    f(t, y, dydt);

    if (callback)
      callback.onIterationStart();

    for (let i = 0; i < dim; ++i)
      yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;

    if (t + h > t1) {
      h = t1 - t;
      flag = false;
    }

    // adaptive step loop
    while (true) {
      if (!isStiff) {
        // ============ ADAMS MODE (non-stiff) ============

        // history invalidation on step size change
        if (adamsHistCount > 1 && h !== lastAdamsH) {
          adamsHistCount = 1;
          f(t, y, fBuf[0]);
        }

        if (adamsHistCount < ADAMS_HIST_LEN) {
          // BOOTSTRAP PHASE: RKF45 single step
          f(t, y, k1);

          for (let i = 0; i < dim; ++i)
            yTemp[i] = y[i] + h * RK_A21 * k1[i];
          f(t + RK_C2 * h, yTemp, k2);

          for (let i = 0; i < dim; ++i)
            yTemp[i] = y[i] + h * (RK_A31 * k1[i] + RK_A32 * k2[i]);
          f(t + RK_C3 * h, yTemp, k3);

          for (let i = 0; i < dim; ++i)
            yTemp[i] = y[i] + h * (RK_A41 * k1[i] + RK_A42 * k2[i] + RK_A43 * k3[i]);
          f(t + RK_C4 * h, yTemp, k4);

          for (let i = 0; i < dim; ++i)
            yTemp[i] = y[i] + h * (RK_A51 * k1[i] + RK_A52 * k2[i] + RK_A53 * k3[i] + RK_A54 * k4[i]);
          f(t + h, yTemp, k5);

          for (let i = 0; i < dim; ++i) {
            yTemp[i] = y[i] + h * (RK_A61 * k1[i] + RK_A62 * k2[i] + RK_A63 * k3[i] +
              RK_A64 * k4[i] + RK_A65 * k5[i]);
          }
          f(t + RK_C6 * h, yTemp, k6);

          // 4th-order solution
          for (let i = 0; i < dim; ++i)
            yTemp[i] = y[i] + h * (RK_B1 * k1[i] + RK_B3 * k3[i] + RK_B4 * k4[i] + RK_B5 * k5[i]);

          // Error estimate
          for (let i = 0; i < dim; ++i)
            yErr[i] = h * (RK_E1 * k1[i] + RK_E3 * k3[i] + RK_E4 * k4[i] + RK_E5 * k5[i] + RK_E6 * k6[i]);

          errmax = 0;
          for (let i = 0; i < dim; ++i)
            errmax = max(errmax, abs(yErr[i] / yScale[i]));
          errmax /= tolerance;

          if (errmax > 1) {
            hTemp = SAFETY * h * errmax**PSHRNK;
            h = max(hTemp, REDUCE_COEF * h);
            tNew = t + h;
            if (tNew == t)
              throw new Error(ERROR_MSG.LSODA_FAILS);
            if (adamsHistCount > 1) {
              adamsHistCount = 1;
              f(t, y, fBuf[0]);
            }
            consecutiveAdamsRejects++;
            if (consecutiveAdamsRejects >= STIFF_DETECT_THRESHOLD) {
              switchToStiff();
              continue;
            }
          } else {
            t = t + h;
            for (let i = 0; i < dim; ++i)
              y[i] = yTemp[i];

            // shift history and store new f-value
            for (let j = Math.min(adamsHistCount, ADAMS_HIST_LEN - 1); j > 0; --j)
              fBuf[j].set(fBuf[j - 1]);
            f(t, y, fBuf[0]);
            adamsHistCount++;
            lastAdamsH = h;
            consecutiveAdamsRejects = 0;

            hNext = h;
            if (hNext > hMax)
              hNext = hMax;
            break;
          }
        } else {
          // AB4-AM3 PECE step

          // Predict (Adams-Bashforth 4-step)
          for (let i = 0; i < dim; ++i) {
            yTemp[i] = y[i] + h * (AB0 * fBuf[0][i] + AB1 * fBuf[1][i] +
              AB2 * fBuf[2][i] + AB3 * fBuf[3][i]);
          }

          // Evaluate f at predicted point
          f(t + h, yTemp, fPred);

          // Correct (Adams-Moulton 3-step)
          for (let i = 0; i < dim; ++i) {
            yCorr[i] = y[i] + h * (AM_NEW * fPred[i] + AM0 * fBuf[0][i] +
              AM1 * fBuf[1][i] + AM2 * fBuf[2][i]);
          }

          // Error estimate (Milne's device)
          for (let i = 0; i < dim; ++i)
            yErr[i] = MILNE * (yCorr[i] - yTemp[i]);

          errmax = 0;
          for (let i = 0; i < dim; ++i)
            errmax = max(errmax, abs(yErr[i] / yScale[i]));
          errmax /= tolerance;

          if (errmax > 1) {
            // step rejected: shrink step and re-bootstrap
            hTemp = SAFETY * h * errmax**PSHRNK;
            h = max(hTemp, REDUCE_COEF * h);
            tNew = t + h;
            if (tNew == t)
              throw new Error(ERROR_MSG.LSODA_FAILS);
            adamsHistCount = 1;
            f(t, y, fBuf[0]);
            consecutiveAdamsRejects++;
            if (consecutiveAdamsRejects >= STIFF_DETECT_THRESHOLD) {
              switchToStiff();
              continue;
            }
          } else {
            if (errmax > ERR_CONTR)
              hNext = SAFETY * h * errmax**PSGROW;
            else
              hNext = GROW_COEF * h;

            if (hNext > hMax)
              hNext = hMax;

            t = t + h;
            for (let i = 0; i < dim; ++i)
              y[i] = yCorr[i];

            // shift history and store f at corrected point
            const recycled = fBuf[ADAMS_HIST_LEN - 1];
            for (let j = ADAMS_HIST_LEN - 1; j > 0; --j)
              fBuf[j] = fBuf[j - 1];
            fBuf[0] = recycled;
            f(t, y, fBuf[0]);
            lastAdamsH = h;
            consecutiveAdamsRejects = 0;

            break;
          }
        }
      } else {
        // ============ BDF MODE (stiff) ============

        const gamma = bdfHistCount >= 2 ? BDF2_GAMMA : 1.0;

        // Jacobian update
        if (!jacobianValid || jacobianAge > 20 || h !== lastBdfH) {
          jacobian(t, y, f, EPS, f0Buf, f1Buf, Jac);
          jacobianValid = true;
          jacobianAge = 0;
          lastBdfH = h;
        }

        // Form W = I - gamma*h*J
        const ghCoeff = gamma * h;
        for (let i = 0; i < dimSquared; ++i)
          W[i] = Ident[i] - ghCoeff * Jac[i];

        // LU decompose
        if (toUseLU)
          luDecomp(W, L, U, dim);

        // Newton predictor
        if (bdfHistCount >= 2) {
          // BDF2 extrapolation: y^(0) = 2*y_n - y_{n-1}
          for (let i = 0; i < dim; ++i)
            yNewton[i] = 2.0 * y[i] - yBdfPrev[i];
        } else {
          // Explicit Euler: y^(0) = y_n + h*f(t_n, y_n)
          f(t, y, fNewton);
          for (let i = 0; i < dim; ++i)
            yNewton[i] = y[i] + h * fNewton[i];
        }

        // Newton iteration
        let newtonConverged = false;
        for (let iter = 0; iter < MAX_NEWTON_ITER; ++iter) {
          f(t + h, yNewton, fNewton);

          // Compute residual G = -(y^(k) - BDF_rhs - gamma*h*fNewton)
          if (bdfHistCount >= 2) {
            for (let i = 0; i < dim; ++i)
              G[i] = -(yNewton[i] - BDF2_A0 * y[i] - BDF2_A1 * yBdfPrev[i] - gamma * h * fNewton[i]);
          } else {
            for (let i = 0; i < dim; ++i)
              G[i] = -(yNewton[i] - y[i] - gamma * h * fNewton[i]);
          }

          // Solve W * delta = G
          if (toUseLU)
            luSolve(L, U, G, luBuf, delta, dim);
          else
            solve1d2d(W, G, delta);

          // Update
          for (let i = 0; i < dim; ++i)
            yNewton[i] += delta[i];

          // Convergence check: ||delta/yScale|| < tolerance
          let deltaNorm = 0;
          for (let i = 0; i < dim; ++i)
            deltaNorm = max(deltaNorm, abs(delta[i] / yScale[i]));

          if (deltaNorm < tolerance) {
            newtonConverged = true;
            break;
          }
        }

        if (!newtonConverged) {
          // Newton failure: treat as step rejection, halve step, invalidate Jacobian
          h = h * 0.5;
          tNew = t + h;
          if (tNew == t)
            throw new Error(ERROR_MSG.LSODA_FAILS);
          jacobianValid = false;
          if (bdfHistCount > 1)
            bdfHistCount = 1;
          continue;
        }

        // Error estimation: compare BDF2 vs BDF1 approximation
        for (let i = 0; i < dim; ++i)
          yBDF1[i] = y[i] + h * fNewton[i];

        for (let i = 0; i < dim; ++i)
          yErr[i] = (yNewton[i] - yBDF1[i]) / 3.0;

        errmax = 0;
        for (let i = 0; i < dim; ++i)
          errmax = max(errmax, abs(yErr[i] / yScale[i]));
        errmax /= tolerance;

        if (errmax > 1) {
          // Step rejected
          hTemp = SAFETY * h * errmax**BDF_PSHRNK;
          h = max(hTemp, REDUCE_COEF * h);
          tNew = t + h;
          if (tNew == t)
            throw new Error(ERROR_MSG.LSODA_FAILS);
          jacobianValid = false;
          if (bdfHistCount > 1)
            bdfHistCount = 1;
        } else {
          // Step accepted
          if (errmax > ERR_CONTR)
            hNext = SAFETY * h * errmax**BDF_PSGROW;
          else
            hNext = GROW_COEF * h;

          // no hMax cap for stiff problems

          t = t + h;

          // Update BDF history
          for (let i = 0; i < dim; ++i)
            yBdfPrev[i] = y[i];
          for (let i = 0; i < dim; ++i)
            y[i] = yNewton[i];

          if (bdfHistCount < 2)
            bdfHistCount++;
          lastBdfH = h;
          jacobianAge++;
          acceptedBdfSteps++;

          // Stiffness re-evaluation: BDF -> Adams
          if (acceptedBdfSteps >= NONSTIFF_TEST_INTERVAL) {
            acceptedBdfSteps = 0;
            if (trialRKF45Step())
              switchToAdams();
          }

          break;
        }
      }
    } // while (true) — adaptive step

    // compute linearly interpolated results and store them in dataframe
    while (timeDataframe < t) {
      const cLeft = (t - timeDataframe) / (t - tPrev);
      const cRight = 1.0 - cLeft;

      tArr[index] = timeDataframe;

      for (let j = 0; j < dim; ++j)
        yArrs[j][index] = cRight * y[j] + cLeft * yPrev[j];

      timeDataframe += hDataframe;
      ++index;
    }

    h = hNext;
    tPrev = t;

    for (let i = 0; i < dim; ++i)
      yPrev[i] = y[i];
  } // while (flag)

  // perform final callback actions
  if (callback)
    callback.onComputationsCompleted();

  // 3. solution at the point t1
  tArr[rowCount - 1] = t1;
  for (let i = 0; i < dim; ++i)
    yArrs[i][rowCount - 1] = y[i];

  // 4. prepare output
  const solution = Array<Float64Array>(dim + 1);
  solution[0] = tArr;
  for (let i = 0; i < dim; ++i)
    solution[i + 1] = yArrs[i];

  return solution;

  // --- Helper functions (closures) ---

  /** Switch to stiff (BDF) mode */
  function switchToStiff(): void {
    isStiff = true;
    acceptedBdfSteps = 0;
    jacobianValid = false;
    // Initialize yBdfPrev from yPrev (already available)
    for (let i = 0; i < dim; ++i)
      yBdfPrev[i] = yPrev[i];
    bdfHistCount = 2; // both y_n and y_{n-1} available
  }

  /** Switch to non-stiff (Adams) mode */
  function switchToAdams(): void {
    isStiff = false;
    consecutiveAdamsRejects = 0;
    adamsHistCount = 1;
    f(t, y, fBuf[0]); // re-bootstrap needed
  }

  /** Trial RKF45 step to check if problem is no longer stiff.
   *  Does NOT advance the solution (diagnostic only).
   *  @returns true if the explicit step succeeds (non-stiff)
   */
  function trialRKF45Step(): boolean {
    f(t, y, k1);

    for (let i = 0; i < dim; ++i)
      yTemp[i] = y[i] + h * RK_A21 * k1[i];
    f(t + RK_C2 * h, yTemp, k2);

    for (let i = 0; i < dim; ++i)
      yTemp[i] = y[i] + h * (RK_A31 * k1[i] + RK_A32 * k2[i]);
    f(t + RK_C3 * h, yTemp, k3);

    for (let i = 0; i < dim; ++i)
      yTemp[i] = y[i] + h * (RK_A41 * k1[i] + RK_A42 * k2[i] + RK_A43 * k3[i]);
    f(t + RK_C4 * h, yTemp, k4);

    for (let i = 0; i < dim; ++i)
      yTemp[i] = y[i] + h * (RK_A51 * k1[i] + RK_A52 * k2[i] + RK_A53 * k3[i] + RK_A54 * k4[i]);
    f(t + h, yTemp, k5);

    for (let i = 0; i < dim; ++i) {
      yTemp[i] = y[i] + h * (RK_A61 * k1[i] + RK_A62 * k2[i] + RK_A63 * k3[i] +
        RK_A64 * k4[i] + RK_A65 * k5[i]);
    }
    f(t + RK_C6 * h, yTemp, k6);

    for (let i = 0; i < dim; ++i)
      yErr[i] = h * (RK_E1 * k1[i] + RK_E3 * k3[i] + RK_E4 * k4[i] + RK_E5 * k5[i] + RK_E6 * k6[i]);

    let trialErr = 0;
    for (let i = 0; i < dim; ++i)
      trialErr = max(trialErr, abs(yErr[i] / yScale[i]));
    trialErr /= tolerance;

    return trialErr <= 1.0;
  }
} // lsoda
