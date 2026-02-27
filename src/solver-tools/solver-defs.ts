// Solver definitions

import {Callback} from './callbacks/callback-base';

/**
 * Represents the right-hand side function **F(t, y)** of an ODE system.
 *
 * @param t - The independent variable.
 * @param y - The current values of the system's state variables.
 * @param output - A preallocated array where the computed derivatives
 *   should be written. This approach avoids unnecessary memory allocations
 *   and improves performance.
 */
export type Func = (t: number, y: Float64Array, output: Float64Array) => void;

/**
 * Optional configuration settings of the numerical solver of an initial value problem (ODE).
 *
 * @property maxIterations - The maximum number of iterations allowed.
 * @property maxTimeMs - The maximum time in milliseconds allowed for the solver.
 * @property method - The name of the numerical method to use.
 */
export type SolverOptions = {
  maxIterations: number,
  maxTimeMs: number,
  method: string,
};

/**
 * Specifies a single-stage ODE model.
 *
 * @property name - The name of the model.
 * @property initial - The initial values of the state variables. Can be
 *   an array or typed array.
 * @property solutionColNames - Names of the columns in the model's output.
 * @property tolerance - Tolerance for the numerical solver.
 * @property func - The function representing the right-hand side F(t, y)
 *   of the ODE system.
 * @property arg - Specification of the independent variable:
 *   - name: the variable's name
 *   - start: starting value
 *   - finish: ending value
 *   - step: step size for the grid
 * * Example: Solve a simple initial value problem (IVP)
 *
 * @example
 * // Example. Solve the following initial value problem:
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
export type ODEs = {
  name: string,
  arg: {
    name: string,
    start: number,
    finish: number,
    step: number
  },
  initial: number[] | Float32Array | Float64Array | Int32Array,
  func: Func,
  tolerance: number,
  solutionColNames: string[],
};

/** The abs function */
export const abs = (x: number) => (x > 0) ? x : -x;

/** The max function */
export const max = (x: number, y: number) => (x > y) ? x : y;

/** The max function */
export const min = (x: number, y: number) => (x < y) ? x : y;

/** Routine constants of adaptive step method */
export const SAFETY = 0.9;
export const PSHRNK = -0.25;
export const PSGROW = -0.2;
export const REDUCE_COEF = 0.25;
export const GROW_COEF = 4.0;
export const ERR_CONTR = 1.89e-4;

/** Misc */
export const TINY = 1e-20;
export const EPS = 1.0e-10;

/** Returns derivative with respect to t. */
export function tDerivative(t: number, y: Float64Array, f: Func, eps: number,
  f0Buf: Float64Array, f1Buf: Float64Array, output: Float64Array): void {
  const size = y.length;
  f(t, y, f0Buf);
  f(t + eps, y, f1Buf);

  for (let i = 0; i < size; ++i)
    output[i] = (f1Buf[i] - f0Buf[i]) / eps;
}

/** Returns Jacobian. */
export function jacobian(t: number, y: Float64Array, f: Func, eps: number,
  f0Buf: Float64Array, f1Buf: Float64Array, output: Float64Array): void {
  const size = y.length;
  f(t, y, f0Buf);

  for (let j = 0; j < size; ++j) {
    y[j] += eps;
    f(t, y, f1Buf);

    for (let i = 0; i < size; ++i)
      output[j + i * size] = (f1Buf[i] - f0Buf[i]) / eps;

    y[j] -= eps;
  }
}

/** Error messeges
 * @internal
*/
export enum ERROR_MSG {
  MRT_FAILS = 'The modified Rosenbrock triple method fails',
  ROS3PRW_FAILS = 'The ROS3PRw method fails',
  ROS34PRW_FAILS = 'The ROS34PRw method fails',
  RK4_FAILS = 'The Runge-Kutta-Fehlberg 4(5) method fails',
  AB5_FAILS = 'The Adams-Bashforth-Moulton 5 method fails',
  AB4_FAILS = 'The Adams-Bashforth-Moulton 4 method fails',
  RKDP_FAILS = 'The Dormand-Prince 5(4) method fails',
  RK3_FAILS = 'The Bogacki-Shampine 3(2) method fails',
  LSODA_FAILS = 'The LSODA method fails',
  CVODE_FAILS = 'The CVODE method fails',
};

/** Callback action
 * @internal
 */
export class CallbackAction extends Error {
  constructor(msg: string) {
    super(msg);
  }
}

/** Default options of the solver
 * @internal
 */
export enum DEFAULT_OPTIONS {
  SCRIPTING = '{maxIterations: 1}',
  NO_CHECKS = '{ }',
}

/** IVP solving method
 * @internal
 */
export type SolverMethod = (odes: ODEs, callback?: Callback) => Float64Array[];
