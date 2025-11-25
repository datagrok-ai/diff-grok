/* eslint-disable max-len */
import {IVP2WebWorker, solveIvp} from '../worker-tools';
import {PipelineCreator} from './pipeline-creator';
import {BasicModelPipelineCreator} from './basic-pipeline-creator';
import {IVP} from '../scripting-tools';
import {UpdatesModelPipelineCreator} from './updates-pipeline-creator';
import {CyclicModelPipelineCreator} from './loops-pipeline-creator';


/**
 * Represents a wrapper for a single modeling step in a multi-stage simulation.
 *
 * @property preproc - A preprocessing script for the step's inputs.
 * @property postproc - A postprocessing script applied after the step is executed.
 * @property out - A script that generates the output of the step.
 */
export type Wrapper = {
    preproc: string | null,
    out: string | null,
    postproc: string | null,
};

/**
 * Represents a computational pipeline for multi-stage simulations.
 *
 * @property wrappers - An array of modeling steps (wrappers) that compose the pipeline.
 * @property out - A script or specification that generates the final output of the pipeline.
 */
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

/** Apply pipeline to initial value problem with the source inputs, and return a solution
 * @param pipeline - computational pipeleine
 * @param ivp - initial value problem
 * @param sourceInputs - source inputs of the model
 * @returns the solution of the initial value problem
 * @example
 *  // Here, we consider modeling queues: https://en.wikipedia.org/wiki/M/M/c_queue
 * import * as DGL from 'diff-grok';
 *
 * // 1. Model specification: the M|M|2|2 model
 * const model = `#name: M|M|2|2
 * #description: Modelling the system with two servers & two waiting places (EXPLAINED)
 *
 * #equations:
 *   dp0/dt = -la * p0 + mu * p1
 *   dp1/dt = la * p0 - (la + mu) * p1 + 2 * mu * p2
 *   dp2/dt = la * p1 - (la + 2 * mu) * p2 + 2 * mu * p3
 *   dp3/dt = la * p2 - (la + 2 * mu) * p3 + 2 * mu * p4
 *
 * #expressions:
 *   la = 1 / arrival
 *   mu = 1 / service
 *   p4 = 1 - p0 - p1 - p2 - p3
 *   busy = 1 - p0 - p1
 *   queue = p3 + 2 * p4
 *
 * #argument: t
 *   _t0 = 0     {min: 0;      max: 10;   caption: start;   category: Time; units: min}  [Initial time of simulation]
 *   _t1 = 60    {min: 20;     max: 100;  caption: finish;  category: Time; units: min}  [Final time of simulation]
 *   _h  = 1     {min: 0.01;   max: 0.1;  caption: step;    category: Time; units: min}  [Time step of simulation]
 *
 * #inits:
 *   p0 = 1   {min: 0; max: 1; category: Initial state; caption: empty}   [Probability that initially there are NO customers]
 *   p1 = 0   {min: 0; max: 1; category: Initial state; caption: single}  [Probability that initially there is ONE customer]
 *   p2 = 0   {min: 0; max: 1; category: Initial state; caption: two}     [Probability that initially there are TWO customers]
 *   p3 = 0   {min: 0; max: 1; category: Initial state; caption: three}   [Probability that initially there are THREE customers]
 *
 * #output:
 *   t     {caption: Time, min}
 *   p0    {caption: P(Empty)}
 *   busy  {caption: P(Full load)}
 *   queue {caption: E(Queue)}
 *
 * #parameters:
 *   arrival = 10   {min: 1; max: 100; category: Means; caption: arrival;  units: min} [Mean arrival time]
 *   service = 100  {min: 1; max: 100; category: Means; caption: service; units: min} [Mean service time]`;
 *
 * // 2. Generate IVP-objects: for the main thread & for computations in webworkers
 * const ivp = DGL.getIVP(model);
 * const ivpWW = DGL.getIvp2WebWorker(ivp);
 *
 *  // 3. Perform computations
 * try {
 *   // 3.1) Extract names of outputs
 *   const outputNames = DGL.getOutputNames(ivp);
 *   const outSize = outputNames.length;
 *
 *   // 3.2) Set model inputs
 *   const inputs = {
 *     _t0: 0, // Initial time of simulation
 *     _t1: 60, // Final time of simulation
 *     _h: 1, // Time step of simulation
 *     p0: 1, // Probability that initially there are NO customers
 *     p1: 0, // Probability that initially there is ONE customer
 *     p2: 0, // Probability that initially there are TWO customers
 *     p3: 0, // Probability that initially there are THREE customers
 *     arrival: 10, // Mean arrival time
 *     service: 100, // Mean service time
 *   };
 *   const inputVector = DGL.getInputVector(inputs, ivp);
 *
 *  // 3.3) Create a pipeline
 *   const creator = DGL.getPipelineCreator(ivp);
 *   const pipeline = creator.getPipeline(inputVector);
 *
 *   // 3.4) Apply pipeline to perform computations
 *   const solution = DGL.applyPipeline(pipeline, ivpWW, inputVector);
 *   // 3.5) Print results
 *
 *   // 3.5.1) Table header
 *   let line = '';
 *   outputNames.forEach((name) => line += name + '      ');
 *   console.log(line);
 *
 *   // 3.5.2) Table with solution
 *   const length = solution[0].length;
 *   for (let i = 0; i < length; ++i) {
 *     line = '';
 *
 *     for (let j = 0; j < outSize; ++j)
 *       line += solution[j][i].toFixed(8) + '     ';
 *     console.log(line);
 *   }
 * } catch (err) {
 *   console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * }
 */
export function applyPipeline(pipeline: Pipeline, ivp: IVP2WebWorker, sourceInputs: Float64Array): Float64Array[] {
  let solution: Float64Array[] = [];
  let currentSolution: Float64Array[];
  let func: Function;
  let inputs = new Float64Array(sourceInputs);
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
      //console.log('Pre');
    }

    //console.log(inputs);

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
      //console.log('Pro');
    }

    //console.log(inputs);
  }

  // Compute final output
  const out = pipeline.out;

  if (out !== null) {
    func = new Function(out);

    return func(solution, sourceInputs);
  }

  return solution;
} // performPipeline

/** Return pipeline creator specified by the initial value problem
 * @param ivp - initial value problem
 * @returns a creator of computational pipeline corresponding to the given model
 * @example
 * // Here, we consider modeling queues: https://en.wikipedia.org/wiki/M/M/c_queue
 * import * as DGL from 'diff-grok';
 *
 * // 1. Model specification: the M|M|2|2 model
 * const model = `#name: M|M|2|2
 * #description: Modelling the system with two servers & two waiting places (EXPLAINED)
 *
 * #equations:
 *   dp0/dt = -la * p0 + mu * p1
 *   dp1/dt = la * p0 - (la + mu) * p1 + 2 * mu * p2
 *   dp2/dt = la * p1 - (la + 2 * mu) * p2 + 2 * mu * p3
 *   dp3/dt = la * p2 - (la + 2 * mu) * p3 + 2 * mu * p4
 *
 * #expressions:
 *   la = 1 / arrival
 *   mu = 1 / service
 *   p4 = 1 - p0 - p1 - p2 - p3
 *   busy = 1 - p0 - p1
 *   queue = p3 + 2 * p4
 *
 * #argument: t
 *   _t0 = 0     {min: 0;      max: 10;   caption: start;   category: Time; units: min}  [Initial time of simulation]
 *   _t1 = 60    {min: 20;     max: 100;  caption: finish;  category: Time; units: min}  [Final time of simulation]
 *   _h  = 1     {min: 0.01;   max: 0.1;  caption: step;    category: Time; units: min}  [Time step of simulation]
 *
 * #inits:
 *   p0 = 1   {min: 0; max: 1; category: Initial state; caption: empty}   [Probability that initially there are NO customers]
 *   p1 = 0   {min: 0; max: 1; category: Initial state; caption: single}  [Probability that initially there is ONE customer]
 *   p2 = 0   {min: 0; max: 1; category: Initial state; caption: two}     [Probability that initially there are TWO customers]
 *   p3 = 0   {min: 0; max: 1; category: Initial state; caption: three}   [Probability that initially there are THREE customers]
 *
 * #output:
 *   t     {caption: Time, min}
 *   p0    {caption: P(Empty)}
 *   busy  {caption: P(Full load)}
 *   queue {caption: E(Queue)}
 *
 * #parameters:
 *   arrival = 10   {min: 1; max: 100; category: Means; caption: arrival;  units: min} [Mean arrival time]
 *   service = 100  {min: 1; max: 100; category: Means; caption: service; units: min} [Mean service time]`;
 *
 * // 2. Generate IVP-objects: for the main thread & for computations in webworkers
 * const ivp = DGL.getIVP(model);
 * const ivpWW = DGL.getIvp2WebWorker(ivp);
 *
 *  // 3. Perform computations
 * try {
 *   // 3.1) Extract names of outputs
 *   const outputNames = DGL.getOutputNames(ivp);
 *   const outSize = outputNames.length;
 *
 *   // 3.2) Set model inputs
 *   const inputs = {
 *     _t0: 0, // Initial time of simulation
 *     _t1: 60, // Final time of simulation
 *     _h: 1, // Time step of simulation
 *     p0: 1, // Probability that initially there are NO customers
 *     p1: 0, // Probability that initially there is ONE customer
 *     p2: 0, // Probability that initially there are TWO customers
 *     p3: 0, // Probability that initially there are THREE customers
 *     arrival: 10, // Mean arrival time
 *     service: 100, // Mean service time
 *   };
 *   const inputVector = DGL.getInputVector(inputs, ivp);
 *
 *  // 3.3) Create a pipeline
 *   const creator = DGL.getPipelineCreator(ivp);
 *   const pipeline = creator.getPipeline(inputVector);
 *
 *   // 3.4) Apply pipeline to perform computations
 *   const solution = DGL.applyPipeline(pipeline, ivpWW, inputVector);
 *   // 3.5) Print results
 *
 *   // 3.5.1) Table header
 *   let line = '';
 *   outputNames.forEach((name) => line += name + '      ');
 *   console.log(line);
 *
 *   // 3.5.2) Table with solution
 *   const length = solution[0].length;
 *   for (let i = 0; i < length; ++i) {
 *     line = '';
 *
 *     for (let j = 0; j < outSize; ++j)
 *       line += solution[j][i].toFixed(8) + '     ';
 *     console.log(line);
 *   }
 * } catch (err) {
 *   console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * }
 */
export function getPipelineCreator(ivp: IVP): PipelineCreator {
  if (ivp.updates !== null)
    return new UpdatesModelPipelineCreator(ivp);

  if (ivp.loop !== null)
    return new CyclicModelPipelineCreator(ivp);

  return new BasicModelPipelineCreator(ivp);
}
