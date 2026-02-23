/* The modified Rosenbrock triple method (MTR).

   References:
   [1] https://doi.org/10.1137/S1064827594276424
   [2] https://doi.org/10.1016/S0898-1221(00)00175-9
*/

import {ODEs, max, min, abs, SAFETY, TINY, EPS, tDerivative, jacobian, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {luDecomp, luSolve, solve1d2d} from './lin-alg-tools';

// Quantities used in Rosenbrock method (see [1], [2] for more details)
const D = 1.0 - Math.sqrt(2.0) / 2.0;
const E32 = 6.0 + Math.sqrt(2.0);
const TWO_D = 2.0 * D;
const ONE_MINUS_2D = 1.0 - TWO_D;

// Adaptive step size control bounds (see MRT.md, Section 6)
const ALPHA_MAX = 5;
const ALPHA_MIN = 0.2;

/** Compute the initial step size using the Hairer-Norsett-Wanner algorithm.
 *  @param t0    initial time
 *  @param y0    initial state vector
 *  @param f     right-hand side function
 *  @param tol   scalar tolerance (plays the role of rtol)
 *  @param dim   dimension of the system
 *  @param p     order of the method (p = 2 for MRT)
 *  @param tEnd  end of the integration interval (used to cap h0)
 *  @returns     estimated initial step size
 */
function initialStepSize(
  t0: number,
  y0: Float64Array,
  f: (t: number, y: Float64Array, out: Float64Array) => void,
  tol: number,
  dim: number,
  p: number,
  tEnd: number,
): number {
  const f0 = new Float64Array(dim);
  const f1 = new Float64Array(dim);
  const y1 = new Float64Array(dim);
  const sc = new Float64Array(dim);

  // Step 1: scale vector (sc[i] = |y0[i]| + TINY)
  for (let i = 0; i < dim; ++i)
    sc[i] = Math.abs(y0[i]) + TINY;

  // Step 2: evaluate the right-hand side at the initial point
  f(t0, y0, f0);

  // Step 3: compute scaled norms d0 = ||y0/sc||, d1 = ||f0/sc||
  let d0 = 0.0;
  let d1 = 0.0;

  for (let i = 0; i < dim; ++i) {
    d0 = Math.max(d0, Math.abs(y0[i]) / sc[i]);
    d1 = Math.max(d1, Math.abs(f0[i]) / sc[i]);
  }

  // Step 4: first guess
  let h0: number;

  if (d0 < 1e-5 || d1 < 1e-5)
    h0 = 1e-6;
  else
    h0 = 0.01 * d0 / d1;

  // Cap h0 to the integration interval
  h0 = Math.min(h0, Math.abs(tEnd - t0));

  // Step 5: explicit Euler probe to estimate the second derivative
  for (let i = 0; i < dim; ++i)
    y1[i] = y0[i] + h0 * f0[i];

  f(t0 + h0, y1, f1);

  let d2 = 0.0;

  for (let i = 0; i < dim; ++i)
    d2 = Math.max(d2, Math.abs(f1[i] - f0[i]) / sc[i]);

  d2 /= h0;

  // Step 6: refine using the method order
  //   h1 = (tol / max(d1, d2))^(1/(p+1))
  const maxD = Math.max(d1, d2);
  let h1: number;

  if (maxD <= 1e-15)
    h1 = Math.max(1e-6, h0 * 1e-3);
  else
    h1 = Math.pow(tol / maxD, 1.0 / (p + 1));

  // Step 7: final initial step size
  h0 = Math.min(100.0 * h0, h1);

  // Safety: clip to the integration interval
  h0 = Math.min(h0, Math.abs(tEnd - t0));

  return h0;
}

/** Solve initial value problem the modified Rosenbrock triple (MRT) method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the MRT method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, mrt } from 'diff-grok';
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
 *   // Solve the problem using the MRT method
 *   const solution = mrt(task);
 *
 *   // Print a table with the results
 *   console.log(task.arg.name, task.solutionColNames[0], task.solutionColNames[1]);
 *   const length = solution[0].length;
 *   for (let i = 0; i < length; ++i)
 *     console.log(solution[0][i], solution[1][i], solution[2][i]);
 * } catch (err) {
 *   console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * }
 * ```
*/
export function mrt(odes: ODEs, callback?: Callback): Float64Array[] {
  /** right-hand side of the IVP solved */
  const f = odes.func;
  // operating variables
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  const hDataframe = odes.arg.step;
  const tolerance = odes.tolerance;

  /** dimension of the problem */
  const dim = odes.initial.length;

  // compute the initial step size automatically
  let h = initialStepSize(t0, new Float64Array(odes.initial), f, tolerance, dim, 2, t1);

  //console.log(`MRT, the problem "${odes.name}", grid step: ${hDataframe}, initial step size: ${h}.`);

  /** number of solution dataframe rows */
  const rowCount = Math.trunc((t1 - t0) / hDataframe) + 1;
  const dimSquared = dim * dim;

  /** independent variable values */
  const tArr = new Float64Array(rowCount);

  /** arrays of solution values */
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
  let tNew = 0;

  // 0 BUFFERS & TEMP STRUCTURES

  /** identity matrix */
  const I = new Float64Array(dim * dim);

  // compute identity matrix
  for (let i = 0; i < dim; ++i) {
    for (let j = 0; j < dim; ++j)
      I[j + i * dim] = (i === j) ? 1 : 0;
  }

  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  const W = new Float64Array(dimSquared);

  const f0 = new Float64Array(dim);
  const k1 = new Float64Array(dim);
  const f1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const f2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const yDer = new Float64Array(dim);
  const hdT = new Float64Array(dim);

  const f0Buf = new Float64Array(dim);
  const f1Buf = new Float64Array(dim);
  let hd = 0;
  let hDivNum = 0;

  const L = new Float64Array(dimSquared);
  const U = new Float64Array(dimSquared);
  const luBuf = new Float64Array(dim);
  const toUseLU = dim > 2;

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

    // call of adaptive step modified Rosenbrok triple method
    // computation of solution (y), time (t) and next step (hNext)
    while (true) {
      // one stage of the modified Rosenbrok triple approach
      // hdT = h * d * T(t, y, EPS);
      tDerivative(t, y, f, EPS, f0Buf, f1Buf, hdT);
      hd = h * D;
      for (let i = 0; i < dim; ++i)
        hdT[i] *= hd;

      // The main computations

      // f0 = f(t, y);
      f(t, y, f0);

      // W = I - h * d * J(t, y, EPS);
      jacobian(t, y, f, EPS, f0Buf, f1Buf, W);
      for (let i = 0; i < dimSquared; ++i)
        W[i] = I[i] - hd * W[i];

      // compute LU-decomposition
      if (toUseLU)
        luDecomp(W, L, U, dim);

      // compute k1: solve the system W * k1 = f0 + hdT
      for (let i = 0; i < dim; ++i)
        f0Buf[i] = f0[i] + hdT[i];

      if (toUseLU)
        luSolve(L, U, f0Buf, luBuf, k1, dim);
      else
        solve1d2d(W, f0Buf, k1);

      hDivNum = 0.5 * h;

      // yDer = y + 0.5 * h * k1;
      for (let i = 0; i < dim; ++i)
        yDer[i] = y[i] + hDivNum * k1[i];

      // f1 = f(t + 0.5 * h, yDer);
      f(t + hDivNum, yDer, f1);

      // compute k2: solve the system W * (k2 - k1) = f1 - k1
      for (let i = 0; i < dim; ++i)
        f1Buf[i] = f1[i] - k1[i];

      if (toUseLU)
        luSolve(L, U, f1Buf, luBuf, k2, dim);
      else
        solve1d2d(W, f1Buf, k2);

      for (let i = 0; i < dim; ++i)
        k2[i] = k2[i] + k1[i];

      // yOut = y + k2 * h; <--> yTemp
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * k2[i];

      // f2 = f(t + h, yOut);
      f(t + h, yTemp, f2);

      // compute k3: solve the system W * k3 = f2 - e32 * (k2 - f1) - 2.0 * (k1 - f0) + hdT
      for (let i = 0; i < dim; ++i)
        f1Buf[i] = f2[i] - E32 * (k2[i] - f1[i]) - 2.0 * (k1[i] - f0[i]) + hdT[i];

      if (toUseLU)
        luSolve(L, U, f1Buf, luBuf, k3, dim);
      else
        solve1d2d(W, f1Buf, k3);

      // yErr = (k1 - 2.0 * k2 + k3) * h / 6;
      hDivNum = h / 6;

      for (let i = 0; i < dim; ++i)
        yErr[i] = (k1[i] - 2.0 * k2[i] + k3[i]) * hDivNum;

      // estimating error
      errmax = 0;
      for (let i = 0; i < dim; ++i)
        errmax = max(errmax, abs(yErr[i] / yScale[i]));
      errmax /= tolerance;

      // adaptive step size control (see MRT.md, Section 6)
      if (errmax > 1) {
        // step rejected: reduce step size, but not below ALPHA_MIN factor
        h = h * max(ALPHA_MIN, SAFETY * errmax ** (-1 / 3));
        tNew = t + h;
        if (tNew == t)
          throw new Error(ERROR_MSG.MRT_FAILS);
      } else {
        // step accepted: grow step size, but not above ALPHA_MAX factor
        hNext = h * min(ALPHA_MAX, SAFETY * errmax ** (-1 / 3));
        t = t + h;

        for (let i = 0; i < dim; ++i)
          y[i] = yTemp[i];

        break;
      }
    } // while (true)

    // compute dense output using Shampine-Reichelt continuous extension
    const hStep = t - tPrev;
    while (timeDataframe < t) {
      const s = (timeDataframe - tPrev) / hStep;
      const b1 = s * (1.0 - s) / ONE_MINUS_2D;
      const b2 = s * (s - TWO_D) / ONE_MINUS_2D;

      tArr[index] = timeDataframe;

      for (let j = 0; j < dim; ++j)
        yArrs[j][index] = yPrev[j] + hStep * (b1 * k1[j] + b2 * k2[j]);

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

  // independent variable
  solution[0] = tArr;

  // functions
  for (let i = 0; i < dim; ++i)
    solution[i + 1] = yArrs[i];

  return solution;
  //return [tArr].concat(yArrs);
} // MTR
