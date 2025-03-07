import {IVP} from '../scripting-tools';
import {ARG_INP_COUNT, argName2IdxMap, LOOP_PARAMS_CONT, LOOP_COUNT_NAME} from './constants';

/** Return inputs in a form of a vector */
export function getInputVector(inputs: Record<string, number>, ivp: IVP): Float64Array {
  const eqsCount = ivp.deqs.solutionNames.length;
  const paramsCount = (ivp.params !== null) ? ivp.params.size : 0;
  const loopParamsCount = (ivp.loop !== null) ? LOOP_PARAMS_CONT : 0;
  const size = ARG_INP_COUNT + eqsCount + paramsCount;
  const inputVector = new Float64Array(size);

  // Argument
  argName2IdxMap.forEach((index, key) => {
    if (key in inputs)
      inputVector[index] = inputs[key];
    else
      throw new Error(`Inconsistent inputs: "${key}" is missing`);
  });

  let idx = ARG_INP_COUNT;

  // Init values
  ivp.deqs.solutionNames.forEach((name) => {
    if (name in inputs)
      inputVector[idx] = inputs[name];
    else
      throw new Error(`Inconsistent inputs: "${name}" is missing`);

    ++idx;
  });

  // Parameters
  if (ivp.params !== null) {
    ivp.params.forEach((_, name) => {
      if (name in inputs)
        inputVector[idx] = inputs[name];
      else
        throw new Error(`Inconsistent inputs: "${name}" is missing`);

      ++idx;
    });
  }

  // Loop params
  if (ivp.loop !== null) {
    if (LOOP_COUNT_NAME in inputs)
      inputVector[idx] = inputs[LOOP_COUNT_NAME];
    else
      throw new Error(`Inconsistent inputs, loop repetitions count is not defined: "${LOOP_COUNT_NAME}" is missing`);
  }

  return inputVector;
} // getInputVector
