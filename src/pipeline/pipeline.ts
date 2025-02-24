import {IVP2WebWorker, solveIvp} from '../worker-tools';

/** Solution step wrapper */
export type Wrapper = {
    preproc: string | null,
    out: string | null,
    postproc: string | null,
};

/** Solution pipeline */
export type Pipeline = {
    wrappers: Wrapper[],
    out: string | null,
};

/** Return concatenated arrays */
function concat(arr1: Float64Array, arr2: Float64Array): Float64Array {
  const result = new Float64Array(arr1.length + arr2.length);
  result.set(arr1, 0);
  result.set(arr2, arr1.length);
  return result;
}

/** Apply pipeline to initial value problem with the source inputs, and return a solution */
export function applyPipeline(pipeline: Pipeline, ivp: IVP2WebWorker, sourceInputs: Float64Array): Float64Array[] {
  let solution: Float64Array[] = [];
  let currentSolution: Float64Array[];
  let func: Function;
  let inputs = sourceInputs;
  const wrappers = pipeline.wrappers;
  let step: Wrapper;
  let isFirstSolverCall = true;

  const stepsCount = wrappers.length;

  if (stepsCount < 1)
    throw new Error(`Incorrect solution pipeline: ${stepsCount} steps.`);

  // Perform pipeline steps
  for (let idx = 0; idx < stepsCount; ++idx) {
    step = wrappers[idx];

    // Pre-processing
    if (step.preproc !== null) {
      func = new Function(step.preproc);
      inputs = func(solution, inputs);
    }

    // Solving a problem
    currentSolution = solveIvp(ivp, inputs);

    // Output from solution
    if (step.out !== null) {
      func = new Function(step.out);
      currentSolution = func(currentSolution, inputs);
    }

    if (isFirstSolverCall) {
      solution = currentSolution;
      isFirstSolverCall = false;
    } else
      currentSolution.forEach((arr, kdx) => solution[kdx] = concat(solution[kdx], arr));

    // Post-processing
    if (step.postproc !== null) {
      func = new Function(step.postproc);
      inputs = func(solution, inputs);
    }
  }

  // Compute final output
  const out = pipeline.out;

  if (out !== null) {
    func = new Function(out);

    return func(solution, sourceInputs);
  }

  return solution;
} // performPipeline
