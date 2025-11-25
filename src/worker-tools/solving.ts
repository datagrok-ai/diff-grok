// Tools for solving IVPs in web-workers

import {SolverMethod, mrt, ros3prw, ros34prw, ODEs, getCallback} from '../solver-tools';
import {Func} from '../solver-tools/solver-defs';
import {IVP2WebWorker, getFunc} from './scripting';

const ARG_INIT = 0;
const ARG_FINAL = 1;
const ARG_STEP = 2;
const ARG_INP_COUNT = 3;

/** Return method for solving IVP */
function getMethod(name: string | undefined): SolverMethod {
  switch (name) {
  case 'mrt':
    return mrt;

  case 'ros3prw':
    return ros3prw;

  default:
    return ros34prw;
  }
}

/** Check the correspondence of inputs to IVP */
function checkInputs(ivp: IVP2WebWorker, inputs: Float64Array): void {
  const expected = ARG_INP_COUNT + ivp.deqsCount + ivp.paramNames.length;

  if (expected > inputs.length)
    throw new Error(`Incorrect inputs count, expected: ${expected}, current: ${inputs.length}`);
}

/** Return solution of IVP with a specified inputs
 * @param ivp - initial value problem to be solved
 * @param inputs - vector of input values
 * @returns solution of the given initial value problem
 */
export function solveIvp(ivp: IVP2WebWorker, inputs: Float64Array): Float64Array[] {
  checkInputs(ivp, inputs);
  const nonParamInputsCount = inputs.length - ivp.paramNames.length;

  const funcCode = getFunc(ivp, inputs.slice(nonParamInputsCount));
  const func = new Function(funcCode) as Func;

  const odes: ODEs = {
    name: '',
    arg: {
      name: '',
      start: inputs[ARG_INIT],
      finish: inputs[ARG_FINAL],
      step: inputs[ARG_STEP],
    },
    initial: inputs.slice(ARG_STEP + 1, nonParamInputsCount),
    func: func,
    tolerance: ivp.tolerance,
    solutionColNames: [],
  };

  const method = getMethod(ivp.solverOpts.method);
  const callback = getCallback(ivp.solverOpts);

  return method(odes, callback);
}
