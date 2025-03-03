import {IVP} from '../scripting-tools';
import {ARG_INP_COUNT, argName2IdxMap} from './constants';

export function getInputVector(inputs: Record<string, number>, ivp: IVP): Float64Array {
  const eqsCount = ivp.deqs.solutionNames.length;
  const paramsCount = (ivp.params !== null) ? ivp.params.size : 0;
  const size = ARG_INP_COUNT + eqsCount + paramsCount;
  const inputVector = new Float64Array(size);

  // Argument
  argName2IdxMap.forEach((index, key) => inputVector[index] = inputs[key]);

  let idx = ARG_INP_COUNT;

  // Init values
  ivp.deqs.solutionNames.forEach((name) => {
    inputVector[idx] = inputs[name];
    ++idx;
  });

  // Parameters

  if (ivp.params === null)
    return inputVector;

  ivp.params.forEach((_, name) => {
    inputVector[idx] = inputs[name];
    ++idx;
  });

  return inputVector;
} // getInputVector
