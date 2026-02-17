/* The Bogacki-Shampine 3(2) method (RK3).

   An explicit adaptive Runge-Kutta method of order 3 with an embedded
   2nd-order error estimate for step size control. Uses the FSAL (First
   Same As Last) property for efficiency.

   Reference:
   [1] P. Bogacki and L.F. Shampine, "A 3(2) pair of Runge-Kutta formulas",
       Applied Mathematics Letters, 2(4):321-325, 1989.
*/

import {ODEs, max, abs, SAFETY, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';

// Step size control exponents for a 3rd-order method
const PSHRNK = -1.0 / 3.0; // -1/p = -1/3
const PSGROW = -0.25; // -1/(p+1) = -1/4

// Bogacki-Shampine 3(2) Butcher tableau coefficients

// Nodes (c_i)
const C2 = 1.0 / 2.0;
const C3 = 3.0 / 4.0;
// C4 = 1 (FSAL)

// Matrix coefficients (a_ij)
const A21 = 1.0 / 2.0;

const A32 = 3.0 / 4.0;

// A4i = Bi (FSAL property: 4th-stage coefficients equal the 3rd-order weights)

// 3rd-order weights (b_i) — used to advance the solution
const B1 = 2.0 / 9.0;
const B2 = 1.0 / 3.0;
const B3 = 4.0 / 9.0;
// B4 = 0

// Error coefficients: (b_i - bHat_i), where bHat are the 2nd-order weights
const E1 = -5.0 / 72.0;
const E2 = 1.0 / 12.0;
const E3 = 1.0 / 9.0;
const E4 = -1.0 / 8.0;

/** Solve initial value problem using the Bogacki-Shampine 3(2) method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the RK3 method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, rk3 } from 'diff-grok';
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
 *   // Solve the problem using the RK3 method
 *   const solution = rk3(task);
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
export function rk3(odes: ODEs, callback?: Callback): Float64Array[] {
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

      // Stage 3: k3 = f(t + C3*h, y + h*A32*k2)
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * A32 * k2[i];
      f(t + C3 * h, yTemp, k3);

      // 3rd-order solution
      for (let i = 0; i < dim; ++i)
        yTemp[i] = y[i] + h * (B1 * k1[i] + B2 * k2[i] + B3 * k3[i]);

      // Stage 4: k4 = f(t + h, yTemp) — needed for error estimate (FSAL)
      f(t + h, yTemp, k4);

      // Error estimate: difference between 3rd and 2nd order solutions
      for (let i = 0; i < dim; ++i)
        yErr[i] = h * (E1 * k1[i] + E2 * k2[i] + E3 * k3[i] + E4 * k4[i]);

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
          throw new Error(ERROR_MSG.RK3_FAILS);
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
} // rk3
