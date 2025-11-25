/* eslint-disable max-len */
// Scripting tools for solving IVPs in web-workers

import {IVP, MATH_CONSTS, POW_IDX, MATH_FUNCS} from '../scripting-tools/scripting-tools';
import {SolverOptions} from '../solver-tools/solver-defs';

/**
 * Represents a fully parsed simulation prepared for execution in Web Workers.
 *
 * @property deqsCount - The number of differential equations in the system.
 * @property funcMainBody - The main body of the function implementing the ODE system.
 * @property paramNames - Names of the model parameters.
 * @property solverOpts - Configuration options for the numerical solver (partial).
 * @property tolerance - Tolerance for the numerical solver.
 */
export type IVP2WebWorker = {
    tolerance: number,
    solverOpts: Partial<SolverOptions>,
    funcMainBody: string,
    paramNames: string[],
    deqsCount: number,
};

const DEFAULT_METHOD = 'ros34prw';
const DEFAULT_MAX_TIME_MS = 50000;

/** Return solver options */
function getSolverOpts(opts: string): Partial<SolverOptions> {
  try {
    const json = JSON.parse(opts
      .replaceAll(`'`, `"`)
      .replace('method', '"method"')
      .replace('maxTimeMs', '"maxTimeMs"')
      .replace('maxIterations', '"maxIterations"')
      .replace(';', ','),
    );

    return json;
  } catch (e) {
    return {
      method: DEFAULT_METHOD,
      maxTimeMs: DEFAULT_MAX_TIME_MS,
    };
  }
}

/** Return main body of a right-hand side function of IVP */
function getFuncMainBody(ivp: IVP): string {
  const lines: string[] = [];

  // Used math funcs
  ivp.usedMathFuncs.forEach((idx) => {
    if (idx !== POW_IDX)
      lines.push(`const ${MATH_FUNCS[idx]} = (x) => Math.${MATH_FUNCS[idx]}(x);`);
    else
      lines.push(`const pow = (x, y) => Math.pow(x, y);`);
  });

  // Used math consts
  ivp.usedMathConsts.forEach((idx) => lines.push(`const ${MATH_CONSTS[idx]} = Math.${MATH_CONSTS[idx]};`));

  // Model constants
  if (ivp.consts !== null)
    ivp.consts.forEach((input, name) => lines.push(`const ${name} = ${input.value};`));

  // extract arg & function values
  lines.push(`const ${ivp.arg.name} = arguments[0];`);
  ivp.deqs.solutionNames.forEach((name, idx) => lines.push(`const ${name} = arguments[1][${idx}];`));

  // evaluate expressions
  if (ivp.exprs !== null)
    ivp.exprs.forEach(((expr, name) => lines.push(`const ${name} = ${expr};`)));

  // compute output
  ivp.deqs.solutionNames.forEach((name, idx) => lines.push(`arguments[2][${idx}] = ${ivp.deqs.equations.get(name)};`));

  return lines.join('\n');
}

/** Return IVP-structure to be passed to web-worker
 * @param ivp - initial value problem to be solved
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
export function getIvp2WebWorker(ivp: IVP): IVP2WebWorker {
  const paramNames: string[] = [];
  if (ivp.params !== null)
    ivp.params.forEach((_, name) => paramNames.push(name));

  const ivpWW: IVP2WebWorker = {
    tolerance: Number(ivp.tolerance),
    solverOpts: getSolverOpts(ivp.solverSettings),
    funcMainBody: getFuncMainBody(ivp),
    paramNames: paramNames,
    deqsCount: ivp.deqs.solutionNames.length,
  };

  return ivpWW;
}

/** Return a right-hand side function of IVP with the specified values of parameters */
export function getFunc(ivp: IVP2WebWorker, paramVals: Float64Array): string {
  const lines: string[] = [];

  ivp.paramNames.forEach((name, idx) => lines.push(`const ${name} = ${paramVals[idx]};`));

  return lines.concat(ivp.funcMainBody).join('\n');
}
