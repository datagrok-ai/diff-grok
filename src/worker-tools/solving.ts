import {SolverMethod, mrt, ros3prw, ros34prw, ODEs, getCallback} from '../solver-tools';
import {Func} from '../solver-tools/solver-defs';
import {IVP2WebWorker, getFunc} from './scripting';

const ARG_INIT = 0;
const ARG_FINAL = 1;
const ARG_STEP = 2;

function getMethod(name: string): SolverMethod {
    switch (name) {
        case 'mrt':
            return mrt;
        
        case 'ros3prw':
            return ros3prw;
    
        default:
            return ros34prw;
    }
}

function checkInputs(ivp: IVP2WebWorker, inputs: Float64Array): void {
    const expected = ivp.arg.vals.length + ivp.initVals.vals.length + ivp.params.vals.length;

    if (expected !== inputs.length)
      throw new Error(`Incorrect inputs count, expected: ${expected}, current: ${inputs.length}`);
}

export function solveIvp(ivp: IVP2WebWorker, inputs: Float64Array, paramsCount: number): Float64Array[] {
  checkInputs(ivp, inputs);
  const nonParamInputsCount = inputs.length - paramsCount;

  const funcCode = getFunc(ivp, inputs.slice(nonParamInputsCount));
  const func = new Function(funcCode) as Func;

  const odes: ODEs = {
    name: ivp.name,
    arg: {
      name: ivp.arg.name,
      start: inputs[ARG_INIT],
      finish: inputs[ARG_FINAL],
      step: inputs[ARG_STEP],
    },
    initial: inputs.slice(ARG_STEP + 1, nonParamInputsCount),
    func: func,
    tolerance: ivp.tolerance,
    solutionColNames: ivp.initVals.names,
  };

  const method = getMethod(ivp.solverOpts.method);
  const callback = getCallback(ivp.solverOpts);

  return method(odes, callback);
}
