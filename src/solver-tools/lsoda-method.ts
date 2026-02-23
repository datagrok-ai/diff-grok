/* LSODA library wrapper method.

   Wraps the variable-order Nordsieck-based LSODA solver into the standard SolverMethod interface.
   Automatically switches between Adams (non-stiff) and BDF (stiff) methods.

   References:
   [1] Hindmarsh, A.C. "ODEPACK, A Systematized Collection of ODE Solvers."
   [2] Petzold, L.R. "Automatic Selection of Methods for Solving Stiff and Nonstiff
       Systems of Ordinary Differential Equations."
*/

import {ODEs, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {Lsoda} from './lsoda';

/** Wraps an ODEs func into the OdeFunction signature compatible with Lsoda (returns 0). */
function wrapFunc(f: ODEs['func']): (t: number, y: Float64Array, ydot: Float64Array) => number {
  return (t, y, ydot) => {
    f(t, y, ydot);
    return 0;
  };
}

/** Solve initial value problem using the LSODA library method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
 * @example
 * // Example. Solve the following initial value problem using the LSODA method:
 * //   dx/dt = x + y - t
 * //   dy/dt = x * y + t
 * //   x(0) = 1
 * //   y(0) = -1
 * // on [0, 2] with the step 0.01.
 * import { ODEs, lsoda } from 'diff-grok';
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
 *   // Solve the problem using the LSODA method
 *   const solution = lsoda(task);
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
export function lsoda(odes: ODEs, callback?: Callback): Float64Array[] {
  console.log('Solving with LSODA method...');
  // 1. Extract problem parameters
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  const step = odes.arg.step;
  const tolerance = odes.tolerance;
  const dim = odes.initial.length;
  const mxstep = 50000;

  // 2. Create LSODA solver with dense output enabled
  const rtol = new Float64Array(dim).fill(tolerance);
  const atol = new Float64Array(dim).fill(tolerance);

  const solver = new Lsoda(wrapFunc(odes.func), dim, {
    rtol, atol, itask: 1, dense: true, mxstep,
  });

  // 3. Integrate through checkpoints, collecting Nordsieck snapshots for dense output.
  //    Checkpoints are spaced coarsely (up to 1000) so the solver handles adaptive
  //    stepping internally; mxstep per call stays within bounds.
  if (callback)
    callback.onIterationStart();

  let y: ArrayLike<number> = [...odes.initial];
  let t = t0;

  // Warmup: a tiny initial step bootstraps the BDF method for extremely stiff problems
  const warmupTout = t0 + Math.min(step, 1.0) * 1e-5;
  if (warmupTout > t0 && warmupTout < t1) {
    const wr = solver.solve(y, t, warmupTout);
    y = wr.y;
    t = wr.t;
    if (solver.state <= 0)
      throw new Error(ERROR_MSG.LSODA_FAILS);
  }

  const gridPoints = Math.trunc((t1 - t0) / step) + 1;
  const numCheckpoints = Math.min(gridPoints, 1000);
  const cpStep = (t1 - t0) / numCheckpoints;

  for (let i = 1; i <= numCheckpoints; i++) {
    const tout = (i < numCheckpoints) ? t0 + i * cpStep : t1;
    const result = solver.solve(y, t, tout);
    y = result.y;
    t = result.t;

    if (solver.state <= 0)
      throw new Error(ERROR_MSG.LSODA_FAILS);
  }

  // 4. Use dense output to interpolate on the uniform output grid
  const dense = solver.getDenseOutput();
  const grid = dense.solveOnGrid(t0, t1, step);

  if (callback)
    callback.onComputationsCompleted();

  // 5. Return solution: [tArr, y0Arr, y1Arr, ...]
  const rowCount = grid.t.length;
  const solution = Array<Float64Array>(dim + 1);
  solution[0] = grid.t;
  for (let i = 0; i < dim; ++i)
    solution[i + 1] = grid.y[i];

  // Patch last y values with the exact endpoint from the solver
  for (let i = 0; i < dim; ++i)
    solution[i + 1][rowCount - 1] = y[i];

  return solution;
} // lsoda
