/* CVODE library wrapper method.

   Wraps the CVODE BDF solver into the standard SolverMethod interface.
   Uses variable-order BDF method for stiff problems with dense direct linear solver.

   References:
   [1] Hindmarsh, A.C. et al. "SUNDIALS: Suite of Nonlinear and Differential/Algebraic
       Equation Solvers." ACM TOMS, 31(3), 363-396, 2005.
   [2] Cohen, S.D. and Hindmarsh, A.C. "CVODE, a Stiff/Nonstiff ODE Solver in C."
       Computers in Physics, 10(2), 138-143, 1996.
*/

import {ODEs, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {Cvode} from './cvode/cvode_class';

/** Solve initial value problem using the CVODE BDF method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the CVODE method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, cvode } from 'diff-grok';
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
 *   // Solve the problem using the CVODE method
 *   const solution = cvode(task);
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
export function cvode(odes: ODEs, callback?: Callback): Float64Array[] {
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  const step = odes.arg.step;
  const tolerance = odes.tolerance;
  const dim = odes.initial.length;

  const solver = new Cvode(odes.func, dim, t0, Float64Array.from(odes.initial), {
    lmm: 'bdf',
    rtol: tolerance,
    atol: tolerance,
    maxSteps: 50000,
  });

  // Warmup: a tiny initial step bootstraps the BDF method for extremely stiff problems.
  // Try 1e-k for k from 4 to 10, picking the largest step that succeeds.
  const base = Math.min(step, 1.0);
  let warmupOk = false;
  for (let k = 4; k <= 15; k++) {
    const warmupTout = t0 + base * Math.pow(10, -k);
    if (warmupTout <= t0 || warmupTout >= t1)
      continue;

    const wr = solver.solve(warmupTout);
    if (wr.flag >= 0) {
      warmupOk = true;
      break;
    }

    // No reset available — recreate solver to retry with a smaller warmup step
    solver.reInit(t0, Float64Array.from(odes.initial));
  }

  if (!warmupOk)
    throw new Error(ERROR_MSG.CVODE_FAILS);

  // Build output grid
  const gridPoints = Math.trunc((t1 - t0) / step) + 1;
  const solution = new Array<Float64Array>(dim + 1);
  for (let i = 0; i <= dim; i++)
    solution[i] = new Float64Array(gridPoints);

  // Initial values
  solution[0][0] = t0;
  for (let j = 0; j < dim; j++)
    solution[j + 1][0] = odes.initial[j];

  // Integrate to each grid point using CV_NORMAL mode (interpolates to exact tout)
  for (let i = 1; i < gridPoints; i++) {
    // check whether to go on computations
    if (callback)
      callback.onIterationStart();

    const tout = (i < gridPoints - 1) ? t0 + i * step : t1;
    const result = solver.solve(tout);

    if (result.flag < 0)
      throw new Error(ERROR_MSG.CVODE_FAILS);

    solution[0][i] = result.t;
    for (let j = 0; j < dim; j++)
      solution[j + 1][i] = result.y[j];
  }

  if (callback)
    callback.onComputationsCompleted();

  return solution;
}
