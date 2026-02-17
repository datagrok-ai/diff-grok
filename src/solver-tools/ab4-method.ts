/* The Adams-Bashforth-Moulton predictor-corrector method of order 4 (AB4).

   Uses the 4-step Adams-Bashforth formula as predictor and the 3-step
   Adams-Moulton formula as corrector in PECE mode (Predict-Evaluate-
   Correct-Evaluate), with Milne's device for local error estimation.
   Bootstrapping is performed using the Runge-Kutta-Fehlberg 4(5) method.

   References:
   [1] E. Hairer, S.P. Norsett, G. Wanner, "Solving Ordinary Differential
       Equations I: Nonstiff Problems", Springer, 2nd ed., 1993.
   [2] L.F. Shampine, M.K. Gordon, "Computer Solution of Ordinary Differential
       Equations: The Initial Value Problem", W.H. Freeman, 1975.
*/

import {ODEs, max, abs, SAFETY, PSHRNK, PSGROW, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';

// Adams-Bashforth 4-step coefficients (pre-divided by 24)
// y_{n+1} = y_n + h * (AB0*f_n + AB1*f_{n-1} + AB2*f_{n-2} + AB3*f_{n-3})
const AB0 = 55.0 / 24.0;
const AB1 = -59.0 / 24.0;
const AB2 = 37.0 / 24.0;
const AB3 = -9.0 / 24.0;

// Adams-Moulton 3-step coefficients (pre-divided by 24)
// y_{n+1} = y_n + h * (AM_NEW*f_{n+1} + AM0*f_n + AM1*f_{n-1} + AM2*f_{n-2})
const AM_NEW = 9.0 / 24.0;
const AM0 = 19.0 / 24.0;
const AM1 = -5.0 / 24.0;
const AM2 = 1.0 / 24.0;

// Milne error estimation coefficient for the AB4-AM3 pair: 19/270
const MILNE = 19.0 / 270.0;

// Number of uniformly-spaced history points required for AB4
const HIST_LEN = 4;

// Fehlberg 4(5) Butcher tableau coefficients (used during bootstrap phase)
const C2 = 1.0 / 4.0;
const C3 = 3.0 / 8.0;
const C4 = 12.0 / 13.0;
const C6 = 1.0 / 2.0;

const A21 = 1.0 / 4.0;
const A31 = 3.0 / 32.0;
const A32 = 9.0 / 32.0;
const A41 = 1932.0 / 2197.0;
const A42 = -7200.0 / 2197.0;
const A43 = 7296.0 / 2197.0;
const A51 = 439.0 / 216.0;
const A52 = -8.0;
const A53 = 3680.0 / 513.0;
const A54 = -845.0 / 4104.0;
const A61 = -8.0 / 27.0;
const A62 = 2.0;
const A63 = -3544.0 / 2565.0;
const A64 = 1859.0 / 4104.0;
const A65 = -11.0 / 40.0;

const B1 = 25.0 / 216.0;
const B3 = 1408.0 / 2565.0;
const B4 = 2197.0 / 4104.0;
const B5 = -1.0 / 5.0;

const RKE1 = 1.0 / 360.0;
const RKE3 = -128.0 / 4275.0;
const RKE4 = -2197.0 / 75240.0;
const RKE5 = 1.0 / 50.0;
const RKE6 = 2.0 / 55.0;

/** Solve initial value problem using the Adams-Bashforth-Moulton predictor-corrector method of order 4
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the AB4 method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, ab4 } from 'diff-grok';
 *
 * const task: ODEs = {
 *   name: 'Example', // name of your model
 *   arg: {
 *       name: 't',  // name of the argument
 *       start: 0,   // initial value of the argument
 *       finish: 2,  // final value of the argument
 *       step: 0.01, // solution grid step
 *   },
 *   initial: [1, -1], // initial values
 *   func: (t: number, y: Float64Array, output: Float64Array) => { // right-hand side of the system
 *     output[0] = y[0] + y[1] - t; // 1-st equation
 *     output[1] = y[0] * y[1] + t; // 2-nd equation
 *   },
 *   tolerance: 1e-7, // tolerance
 *   solutionColNames: ['x', 'y'], // names of solution functions
 * };
 *
 * try {
 *   // Solve the problem using the AB4 method
 *   const solution = ab4(task);
 *
 *   // Print a table with the results
 *   console.log(task.arg.name, task.solutionColNames[0], task.solutionColNames[1]);
 *   const length = solution[0].length;
 *   for (let i = 0; i < length; ++i)
 *     console.log(solution[0][i], solution[1][i], solution[2][i]);
 * } catch (err) {
 *   console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * }
*/
export function ab4(odes: ODEs, callback?: Callback): Float64Array[] {
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

  // working buffers
  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  // AB4-AM3 history: fBuf[0] = f(t_n, y_n), fBuf[3] = f(t_{n-3}, y_{n-3})
  const fBuf: Float64Array[] = [];
  for (let i = 0; i < HIST_LEN; ++i)
    fBuf[i] = new Float64Array(dim);
  let histCount = 0;
  let lastH = 0;

  // AB4-AM3 corrector buffers
  const yCorr = new Float64Array(dim);
  const fPred = new Float64Array(dim);

  // RKF45 stage buffers (bootstrap)
  const k1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const k4 = new Float64Array(dim);
  const k5 = new Float64Array(dim);
  const k6 = new Float64Array(dim);

  // store initial f-value
  f(t, y, fBuf[0]);
  histCount = 1;

  // 1. SOLUTION AT THE POINT t0
  tArr[0] = t0;
  for (let i = 0; i < dim; ++i)
    yArrs[i][0] = y[i];

  // 2. COMPUTE NUMERICAL SOLUTION FOR THE POINTS FROM THE INTERVAL (t0, t1)
  while (flag) {
    // compute derivative
    f(t, y, dydt);

    // check whether to go on computations
    if (callback)
      callback.onIterationStart();

    // compute scale vector
    for (let i = 0; i < dim; ++i)
      yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;

    // check end point
    if (t + h > t1) {
      h = t1 - t;
      flag = false;
    }

    // if step size changed, history is invalid (AB4 needs uniform spacing)
    if (histCount > 1 && h !== lastH) {
      histCount = 1;
      f(t, y, fBuf[0]);
    }

    // adaptive step loop
    while (true) {
      if (histCount < HIST_LEN) {
        // BOOTSTRAP PHASE: RKF45 single step to build up history

        // Stage 1: k1 = f(t, y)
        f(t, y, k1);

        // Stage 2: k2 = f(t + C2*h, y + h*A21*k1)
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * A21 * k1[i];
        f(t + C2 * h, yTemp, k2);

        // Stage 3: k3 = f(t + C3*h, y + h*(A31*k1 + A32*k2))
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
        f(t + C3 * h, yTemp, k3);

        // Stage 4: k4 = f(t + C4*h, y + h*(A41*k1 + A42*k2 + A43*k3))
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
        f(t + C4 * h, yTemp, k4);

        // Stage 5: k5 = f(t + h, y + h*(A51*k1 + A52*k2 + A53*k3 + A54*k4))
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
        f(t + h, yTemp, k5);

        // Stage 6: k6 = f(t + C6*h, y + h*(A61*k1 + A62*k2 + A63*k3 + A64*k4 + A65*k5))
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
        f(t + C6 * h, yTemp, k6);

        // 4th-order solution
        for (let i = 0; i < dim; ++i)
          yTemp[i] = y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i]);

        // Error estimate: difference between 4th and 5th order solutions
        for (let i = 0; i < dim; ++i)
          yErr[i] = h * (RKE1 * k1[i] + RKE3 * k3[i] + RKE4 * k4[i] + RKE5 * k5[i] + RKE6 * k6[i]);

        // estimating error
        errmax = 0;
        for (let i = 0; i < dim; ++i)
          errmax = max(errmax, abs(yErr[i] / yScale[i]));
        errmax /= tolerance;

        // processing the error obtained
        if (errmax > 1) {
          hTemp = SAFETY * h * errmax**PSHRNK;
          h = max(hTemp, REDUCE_COEF * h);
          tNew = t + h;
          if (tNew == t)
            throw new Error(ERROR_MSG.AB4_FAILS);
          // step size changed, reset history if needed
          if (histCount > 1) {
            histCount = 1;
            f(t, y, fBuf[0]);
          }
        } else {
          // step accepted
          t = t + h;

          for (let i = 0; i < dim; ++i)
            y[i] = yTemp[i];

          // shift history and store new f-value
          for (let j = Math.min(histCount, HIST_LEN - 1); j > 0; --j)
            fBuf[j].set(fBuf[j - 1]);
          f(t, y, fBuf[0]);
          histCount++;
          lastH = h;

          // during bootstrap, keep h constant for uniform spacing
          hNext = h;
          if (hNext > hMax)
            hNext = hMax;

          break;
        }
      } else {
        // AB4-AM3 PHASE: Predictor-Corrector (PECE)

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

        // estimating error
        errmax = 0;
        for (let i = 0; i < dim; ++i)
          errmax = max(errmax, abs(yErr[i] / yScale[i]));
        errmax /= tolerance;

        // processing the error obtained
        if (errmax > 1) {
          // step rejected: shrink step and re-bootstrap
          hTemp = SAFETY * h * errmax**PSHRNK;
          h = max(hTemp, REDUCE_COEF * h);
          tNew = t + h;
          if (tNew == t)
            throw new Error(ERROR_MSG.AB4_FAILS);
          histCount = 1;
          f(t, y, fBuf[0]);
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

          // shift history and store f at corrected point (PECE final evaluate)
          const recycled = fBuf[HIST_LEN - 1];
          for (let j = HIST_LEN - 1; j > 0; --j)
            fBuf[j] = fBuf[j - 1];
          fBuf[0] = recycled;
          f(t, y, fBuf[0]);
          lastH = h;

          break;
        }
      }
    } // while (true)

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
  const solution = Array<Float64Array>(dim);

  // independent variable
  solution[0] = tArr;

  // functions
  for (let i = 0; i < dim; ++i)
    solution[i + 1] = yArrs[i];

  return solution;
} // ab4
