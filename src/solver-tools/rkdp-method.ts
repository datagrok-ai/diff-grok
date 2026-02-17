/* The Dormand-Prince 5(4) method (RKDP).

   An explicit adaptive Runge-Kutta method of order 5 with an embedded
   4th-order error estimate for step size control. Uses the FSAL (First
   Same As Last) property for efficiency.

   Reference:
   [1] J.R. Dormand and P.J. Prince, "A family of embedded Runge-Kutta formulae",
       Journal of Computational and Applied Mathematics, 6(1):19-26, 1980.
*/

import {ODEs, max, abs, SAFETY, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';

// Step size control exponents for a 5th-order method
const PSHRNK = -0.2; // -1/p = -1/5
const PSGROW = -1.0 / 6.0; // -1/(p+1) = -1/6

// Dormand-Prince 5(4) Butcher tableau coefficients

// Nodes (c_i)
const C2 = 1.0 / 5.0;
const C3 = 3.0 / 10.0;
const C4 = 4.0 / 5.0;
const C5 = 8.0 / 9.0;
// C6 = 1
// C7 = 1

// Matrix coefficients (a_ij)
const A21 = 1.0 / 5.0;

const A31 = 3.0 / 40.0;
const A32 = 9.0 / 40.0;

const A41 = 44.0 / 45.0;
const A42 = -56.0 / 15.0;
const A43 = 32.0 / 9.0;

const A51 = 19372.0 / 6561.0;
const A52 = -25360.0 / 2187.0;
const A53 = 64448.0 / 6561.0;
const A54 = -212.0 / 729.0;

const A61 = 9017.0 / 3168.0;
const A62 = -355.0 / 33.0;
const A63 = 46732.0 / 5247.0;
const A64 = 49.0 / 176.0;
const A65 = -5103.0 / 18656.0;

// A7i = Bi (FSAL property: 7th-stage coefficients equal the 5th-order weights)

// 5th-order weights (b_i) — used to advance the solution
const B1 = 35.0 / 384.0;
// B2 = 0
const B3 = 500.0 / 1113.0;
const B4 = 125.0 / 192.0;
const B5 = -2187.0 / 6784.0;
const B6 = 11.0 / 84.0;
// B7 = 0

// Error coefficients: (b_i - bHat_i), where bHat are the 4th-order weights
const E1 = 71.0 / 57600.0;
// E2 = 0
const E3 = -71.0 / 16695.0;
const E4 = 71.0 / 1920.0;
const E5 = -17253.0 / 339200.0;
const E6 = 22.0 / 525.0;
const E7 = -1.0 / 40.0;

/** Solve initial value problem using the Dormand-Prince 5(4) method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the RKDP method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, rkdp } from 'diff-grok';
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
 *   // Solve the problem using the RKDP method
 *   const solution = rkdp(task);
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
export function rkdp(odes: ODEs, callback?: Callback): Float64Array[] {
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

  // buffers
  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  const k1 = new Float64Array(dim);
  const k2 = new Float64Array(dim);
  const k3 = new Float64Array(dim);
  const k4 = new Float64Array(dim);
  const k5 = new Float64Array(dim);
  const k6 = new Float64Array(dim);
  const k7 = new Float64Array(dim);

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

    // adaptive step loop
    while (true) {
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

      // Stage 5: k5 = f(t + C5*h, y + h*(A51*k1 + A52*k2 + A53*k3 + A54*k4))
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
      f(t + C5 * h, yTemp, k5);

      // Stage 6: k6 = f(t + h, y + h*(A61*k1 + A62*k2 + A63*k3 + A64*k4 + A65*k5))
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
      f(t + h, yTemp, k6);

      // 5th-order solution
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]);

      // Stage 7: k7 = f(t + h, yTemp) — needed for error estimate
      f(t + h, yTemp, k7);

      // Error estimate: difference between 5th and 4th order solutions
      for (let i = 0; i < dim; ++i)
        yErr[i] = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);

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
          throw new Error(ERROR_MSG.RKDP_FAILS);
      } else {
        if (errmax > ERR_CONTR)
          hNext = SAFETY * h * errmax**PSGROW;
        else
          hNext = GROW_COEF * h;

        if (hNext > hMax)
          hNext = hMax;

        t = t + h;

        for (let i = 0; i < dim; ++i)
          y[i] = yTemp[i];

        break;
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
} // rkdp
