/* eslint-disable max-len */
import {IVP} from '../scripting-tools';
import {MATH_CONSTS, MATH_FUNCS, POW_IDX} from '../scripting-tools/scripting-tools';
import {getOutputCode} from './output';
import {Pipeline, Wrapper} from './pipeline';
import {PipelineCreator} from './pipeline-creator';

import {argName2IdxMap, ARG, ARG_COL_IDX, CORRECTION_FACTOR, DEFAULT_CORRECTION, LOOP_MIN_COUNT} from './constants';

/** Pipeline creator for cyclic model
 * @example
 * // This example shows how to apply pipelines and cyclic models.
 * // This approach can be used for in-webworkers analysis of models.
 * // Here, we consider pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model.
 * import * as DGL from 'diff-grok';
 *
 * // 1. Model specification
 * const model = `#name: PK-PD
 * #tags: model
 * #description: Pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model
 * #equations:
 *   d(depot)/dt = -KA * depot
 *   d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
 *   d(peri)/dt  = Q * C2 - Q * C3
 *   d(eff)/dt  = Kin - Kout * (1 - C2/(EC50 + C2)) * eff
 *
 * #expressions:
 *   C2 = centr / V2
 *   C3 = peri / V3
 *
 * #loop:
 *   _count = 10 {caption: count; category: Dosing; min: 1; max: 20} [Number of doses]
 *   depot += dose
 *
 * #argument: t
 *   _t0 = 0 {units: h; caption: begin; category: Dosing; min: 0; max: 1} [Begin of dosing interval]
 *   _t1 = 12 {units: h; caption: end; category: Dosing; min: 5; max: 15} [End of dosing interval]
 *   _h = 1 {units: h; caption: step; category: Dosing; min: 0.01; max: 0.1} [Time step of simulation]
 *
 * #inits:
 *   depot = 0 {category: Initial values}
 *   centr = 0 {category: Initial values} [Central]
 *   peri = 0 {category: Initial values} [Peripheral]
 *   eff = 0.2 {category: Initial values} [Effective compartment rate]
 *
 * #parameters:
 *   dose = 1e4 {category: Dosing; min: 1e3; max: 2e4; step: 1e3} [Dosage]
 *   KA = 0.3 {caption: rate constant; category: Parameters; min: 0.1; max: 1}
 *   CL = 2 {caption: clearance; category: Parameters; min: 1; max: 5}
 *   V2 = 4 {caption: central volume; category: Parameters; min: 1; max: 10} [Central compartment volume]
 *   Q = 1 {caption: inter rate; category: Parameters; min: 0.1; max: 1} [Intercompartmental rate]
 *   V3 = 30 {caption: peri volume; category: Parameters; min: 20; max: 40} [Peripheral compartment volume]
 *   EC50 = 8 {caption: effect; category: Parameters; min: 1; max: 10}
 *   Kin = 0.2 {caption: Kin; category: Parameters; min: 0.1; max: 0.5} [The first-order production constant]
 *   Kout = 0.2 {caption: Kout; category: Parameters; min: 0.1; max: 0.5} [The first-order dissipation rate constant]
 *
 * #tolerance: 1e-9`;
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
 *     _count: 10,
 *     _t0: 0,
 *     _t1: 16,
 *     _h: 1,
 *     depot: 0,
 *     centr: 0,
 *     peri: 0,
 *     eff: 0.2,
 *     dose: 10000,
 *     KA: 0.63,
 *     CL: 3.2,
 *     V2: 6.58,
 *     Q: 0.622,
 *     V3: 35.6,
 *     EC50: 5.41,
 *     Kin: 0.272,
 *     Kout: 0.276,
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
 *   let line = '         ';
 *   outputNames.forEach((name) => line += name + '           ');
 *   console.log(line);
 *
 *   // 3.5.2) Table with solution
 *   const length = solution[0].length;
 *   for (let i = 0; i < length; ++i) {
 *     line = '';
 *
 *     for (let j = 0; j < outSize; ++j)
 *       line += solution[j][i].toFixed(8) + '     ';
 *
 *     console.log(line);
 *   }
 * } catch (err) {
 *   console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
 * } */
export class CyclicModelPipelineCreator extends PipelineCreator {
  /**
   * @param ivp - initial value problem to be solved
   */
  constructor(ivp: IVP) {
    super(ivp);

    if (ivp.loop === null)
      throw new Error('Incorrect use of CyclicModelPipelineCreator: model has no loop');
  }

  /** Return computations pipeline corresponding to the given input vector
   * @param inputs - input vector of the model
   */
  public getPipeline(inputs: Float64Array): Pipeline {
    if (this.ivp.loop === null)
      throw new Error('Incorrect use of CyclicModelPipelineCreator: model has no loop');

    const argSize = argName2IdxMap.size;
    const restricted = inputs.length - 1;
    const repetitions = inputs[restricted];

    if (repetitions < LOOP_MIN_COUNT)
      throw new Error(`Incorrect repetitions count: ${repetitions}, expected: > 0`);

    const preProcLines: string[] = [];

    // Used math funcs
    this.ivp.usedMathFuncs.forEach((idx) => {
      if (idx !== POW_IDX)
        preProcLines.push(`const ${MATH_FUNCS[idx]} = (x) => Math.${MATH_FUNCS[idx]}(x);`);
      else
        preProcLines.push(`const pow = (x, y) => Math.pow(x, y);`);
    });

    // Used math consts
    this.ivp.usedMathConsts.forEach((idx) => {
      preProcLines.push(`const ${MATH_CONSTS[idx]} = Math.${MATH_CONSTS[idx]};`);
    });

    // Model constants
    if (this.ivp.consts !== null)
      this.ivp.consts.forEach((input, name) => preProcLines.push(`const ${name} = ${input.value};`));

    preProcLines.push(`const inputs = arguments[1].slice(0, ${restricted});`);

    // Extract inputs
    argName2IdxMap.forEach((idx, name) => preProcLines.push(`let ${name} = inputs[${idx}];`));

    this.ivp.deqs.solutionNames.forEach((name, idx) => {
      preProcLines.push(`let ${name} = inputs[${idx + argSize}];`);
    });

    if (this.ivp.params !== null) {
      const shift = argSize + this.ivp.deqs.solutionNames.length;
      let idx = 0;

      this.ivp.params.forEach((_, name) => {
        preProcLines.push(`let ${name} = inputs[${idx + shift}];`);
        ++idx;
      });
    }

    // Updates
    this.ivp.loop.updates.forEach((expr) => preProcLines.push(`${expr};`));

    // Input vector update
    argName2IdxMap.forEach((idx, name) => preProcLines.push(`inputs[${idx}] = ${name};`));

    this.ivp.deqs.solutionNames.forEach((name, idx) => {
      preProcLines.push(`inputs[${idx + argSize}] = ${name};`);
    });

    if (this.ivp.params !== null) {
      const shift = argSize + this.ivp.deqs.solutionNames.length;
      let idx = 0;

      this.ivp.params.forEach((_, name) => {
        preProcLines.push(`inputs[${idx + shift}] = ${name};`);
        ++idx;
      });
    }

    preProcLines.push(`return inputs;`);

    const postProcLines = [
      'const solution = arguments[0];',
      'const lastIdx = solution[0].length - 1;',
      'const inputs = new Float64Array(arguments[1]);',
      `const duration = inputs[${argName2IdxMap.get(ARG.FINISH)}] - inputs[${argName2IdxMap.get(ARG.START)}];`,
      `inputs[${argName2IdxMap.get(ARG.START)}] = inputs[${argName2IdxMap.get(ARG.FINISH)}];`,
      `inputs[${argName2IdxMap.get(ARG.FINISH)}] += duration;`,
    ];

    this.ivp.deqs.solutionNames.forEach((_, idx) => {
      postProcLines.push(`inputs[${idx + argSize}] = solution[${idx + 1}][lastIdx];`);
    });

    postProcLines.push('return inputs;');

    const preCode = preProcLines.join('\n');
    const postCode = postProcLines.join('\n');

    const wrappers = new Array<Wrapper>(repetitions);

    for (let i = 0; i < repetitions; ++i) {
      wrappers[i] = {
        preproc: preCode,
        out: null,
        postproc: postCode,
      };
    }

    wrappers[repetitions - 1].postproc = this.getArgumentCorrectionCode();

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
} // CyclicModelPipelineCreator
