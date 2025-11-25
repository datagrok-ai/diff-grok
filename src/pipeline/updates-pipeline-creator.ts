/* eslint-disable max-len */

import {IVP} from '../scripting-tools';
import {MATH_CONSTS, MATH_FUNCS, POW_IDX} from '../scripting-tools/scripting-tools';
import {getOutputCode} from './output';
import {Pipeline, Wrapper} from './pipeline';
import {PipelineCreator} from './pipeline-creator';

import {argName2IdxMap, ARG, ARG_COL_IDX, CORRECTION_FACTOR, DEFAULT_CORRECTION} from './constants';

/** Pipeline creator for model with updates
 * @example
 * // This example shows how to apply pipelines and models with updates.
 * // This approach can be used for in-webworkers analysis of models.
 * // Here, we consider gluconic acid (GA) production by Aspergillus niger modeling.
 *
 * import * as DGL from 'diff-grok';
 *
 * // 1. Model specification
 * const model = `#name: GA-production
 * #tags: model
 * #description: Gluconic acid (GA) production by Aspergillus niger modeling
 *
 * #equations:
 *    dX/dt = rX
 *    dS/dt = -gamma * rX - lambda * X
 *    dO/dt = Kla * (Cod - O) - delta * rX - phi * X
 *    dP/dt = alpha * rX + beta * X
 *
 * #expressions:
 *    mu = muM * S / (Ks + S) * O / (Ko + O)
 *    rX = mu * X
 *
 * #argument: t, 1-st stage
 *    _t0 = 0 {units: h; caption: initial; category: Misc} [Start of the process]
 *    _t1 = 60 {units: h; caption: 1-st stage; category: Durations; min: 20; max: 80} [Duration of the 1-st stage]
 *    step = 0.1 {units: h; caption: step; category: Misc; min: 0.01; max: 1} [Time step of simulation]
 *
 * #update: 2-nd stage
 *    duration = overall - _t1
 *    S += 70
 *
 * #inits:
 *    X = 5 {units: kg/m³; caption: biomass; category: Initial concentrations; min: 1; max: 10} [Aspergillus niger biomass]
 *    S = 150 {units: kg/m³; caption: glucose; category: Initial concentrations; min: 50; max: 200} [Glucose]
 *    O = 7 {units: kg/m³; caption: oxygen; category: Initial concentrations; min: 1; max: 10} [Dissolved oxygen]
 *    P = 0 {units: kg/m³; caption: acid; category: Initial concentrations; min: 0; max: 0.1} [Gluconic acid]
 *
 * #output:
 *    t {caption: time}
 *    X {caption: biomass}
 *    S {caption: glucose}
 *    O {caption: oxygen}
 *    P {caption: acid}
 *
 * #parameters:
 *    overall = 100 {units: h; category: Durations; min: 100; max: 140} [Overall duration]
 *    muM = 0.668 {units: 1/h; category: Parameters} [Monod type model parameter]
 *    alpha = 2.92 {category: Parameters} [Monod type model parameter]
 *    beta = 0.131 {units: 1/h; category: Parameters} [Monod type model parameter]
 *    gamma = 2.12 {category: Parameters} [Monod type model parameter]
 *    lambda = 0.232 {units: 1/h; category: Parameters} [Monod type model parameter]
 *    delta = 0.278 {category: Parameters} [Monod type model parameter]
 *    phi = 4.87e-3 {units: 1/h; category: Parameters} [Monod type model parameter]
 *    Ks = 1.309e2 {units: g/L; category: Parameters} [Monod type model parameter]
 *    Ko = 3.63e-4 {units: g/L; category: Parameters} [Monod type model parameter]
 *    Kla = 1.7e-2 {units: 1/s; category: Parameters} [Volumetric mass transfer coefficient]
 *    Cod = 15 {units: kg/m³; category: Parameters} [Liquid phase dissolved oxygen saturation concentration]
 *
 *  #tolerance: 1e-9`;
 *
 * // 2. Generate IVP-objects: for the main thread & for computations in webworkers
 * const ivp = DGL.getIVP(model);
 * const ivpWW = DGL.getIvp2WebWorker(ivp);
 *
 * // 3. Perform computations
 * try {
 *   // 3.1) Extract names of outputs
 *   const outputNames = DGL.getOutputNames(ivp);
 *   const outSize = outputNames.length;
 *
 *   // 3.2) Set model inputs
 *   const inputs = {
 *     _t0: 0,
 *     _t1: 60,
 *     _h: 1,
 *     X: 5,
 *     S: 150,
 *     O: 7,
 *     P: 0,
 *     overall: 100,
 *     muM: 0.668,
 *     alpha: 2.92,
 *     beta: 0.131,
 *     gamma: 2.12,
 *     lambda: 0.232,
 *     delta: 0.278,
 *     phi: 4.87e-3,
 *     Ks: 1.309e2,
 *     Ko: 3.63e-4,
 *     Kla: 1.7e-2,
 *     Cod: 15,
 *   };
 *
 *   const inputVector = DGL.getInputVector(inputs, ivp);
 *
 *   // 3.3) Create a pipeline
 *   const creator = DGL.getPipelineCreator(ivp);
 *   const pipeline = creator.getPipeline(inputVector);
 *
 *   // 3.4) Apply pipeline to perform computations
 *   const solution = DGL.applyPipeline(pipeline, ivpWW, inputVector);
 *
 *   // 3.5) Print results
 *
 *   // 3.5.1) Table header
 *   let line = '     ';
 *   outputNames.forEach((name) => line += name + '         ');
 *   console.log(line);
 *
 *   // 3.5.2) Table with solution  const length = solution[0].length;
 *   for (let i = 0; i < length; ++i) {
 *       line = '';
 *
 *       for (let j = 0; j < outSize; ++j)
 *         line += solution[j][i].toFixed(8) + '     ';
 *
 *     console.log(line);
 *   }
 * } catch (err) {
 *   console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * }
 */
export class UpdatesModelPipelineCreator extends PipelineCreator {
  /**
   * @param ivp - initial value problem to be solved
   */
  constructor(ivp: IVP) {
    super(ivp);

    if (ivp.updates === null)
      throw new Error('Incorrect use of UpdatesModelPipelineCreator: model has no updates');
  }

  /** Return computations pipeline corresponding to the given input vector
   * @param inputs - input vector of the model
  */
  public getPipeline(inputs: Float64Array): Pipeline {
    if (this.ivp.updates === null)
      throw new Error('Incorrect use of UpdatesModelPipelineCreator: model has no updates');

    const wrappers: Wrapper[] = [{
      preproc: null,
      out: null,
      postproc: null,
    }];

    this.ivp.updates.forEach((update) => {
      const lines: string[] = [];

      // Used math funcs
      this.ivp.usedMathFuncs.forEach((idx) => {
        if (idx !== POW_IDX)
          lines.push(`const ${MATH_FUNCS[idx]} = (x) => Math.${MATH_FUNCS[idx]}(x);`);
        else
          lines.push(`const pow = (x, y) => Math.pow(x, y);`);
      });

      // Used math consts
      this.ivp.usedMathConsts.forEach((idx) => lines.push(`const ${MATH_CONSTS[idx]} = Math.${MATH_CONSTS[idx]};`));

      // Model constants
      if (this.ivp.consts !== null)
        this.ivp.consts.forEach((input, name) => lines.push(`const ${name} = ${input.value};`));

      lines.push('const inputs = new Float64Array(arguments[1]);');
      lines.push('const solution = arguments[0];');

      // Arguments
      argName2IdxMap.forEach((idx, name) => lines.push(`let ${name} = inputs[${idx}];`));

      // Functions
      lines.push('const lastIdx = solution[0].length - 1;');

      this.ivp.deqs.solutionNames.forEach((name, idx) => {
        lines.push(`let ${name} = solution[${idx + 1}][lastIdx];`);
      });

      // Parameters
      if (this.ivp.params !== null) {
        const shift = argName2IdxMap.size + this.ivp.deqs.solutionNames.length;
        let idx = 0;

        this.ivp.params.forEach((_, name) => {
          lines.push(`let ${name} = inputs[${idx + shift}];`);
          ++idx;
        });
      }

      // Updates
      update.updates.forEach((expr) => lines.push(`${expr};`));

      // Argument update
      lines.push(`const duration = ${update.durationFormula};`);
      lines.push(`${ARG.START} = ${ARG.FINISH};`);
      lines.push(`${ARG.FINISH} += duration;`);

      // Input vector update
      argName2IdxMap.forEach((idx, name) => lines.push(`inputs[${idx}] = ${name};`));

      this.ivp.deqs.solutionNames.forEach((name, idx) => {
        lines.push(`inputs[${idx + argName2IdxMap.size}] = ${name};`);
      });

      if (this.ivp.params !== null) {
        const shift = argName2IdxMap.size + this.ivp.deqs.solutionNames.length;
        let idx = 0;

        this.ivp.params.forEach((_, name) => {
          lines.push(`inputs[${idx + shift}] = ${name};`);
          ++idx;
        });
      }

      lines.push('return inputs;');

      wrappers.push({
        preproc: lines.join('\n'),
        out: null,
        postproc: null,
      });
    });

    wrappers[wrappers.length - 1].postproc = this.getArgumentCorrectionCode();

    return {
      wrappers: wrappers,
      out: getOutputCode(this.ivp),
    };
  } // getPipeline

  /** Return code for correction of an argument */
  private getArgumentCorrectionCode(): string {
    return `const arg = arguments[0][${ARG_COL_IDX}];
    const length = arg.length;
    const step = arg[1] - arg[0];
    const correction = Math.min(step * ${CORRECTION_FACTOR}, ${DEFAULT_CORRECTION});
    for (let i = 0; i < length - 1; ++i) {
      if (arg[i + 1] - arg[i] < step)
        arg[i] -= correction;
    }
    return arguments[0];`;
  } // getArgumentCorrectionCode
} // UpdatesModelPipelineCreator
