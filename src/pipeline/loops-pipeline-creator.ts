import {IVP} from '../scripting-tools';
import {MATH_CONSTS, MATH_FUNCS, POW_IDX} from '../scripting-tools/scripting-tools';
import {getOutputCode} from './output';
import {Pipeline, Wrapper} from './pipeline';
import {PipelineCreator} from './pipeline-creator';

import {argName2IdxMap, ARG, ARG_COL_IDX, CORRECTION_FACTOR, DEFAULT_CORRECTION, LOOP_MIN_COUNT} from './constants';

/** Pipeline creator for cyclic model */
export class CyclicModelPipelineCreator extends PipelineCreator {
  constructor(ivp: IVP) {
    super(ivp);

    if (ivp.loop === null)
      throw new Error('Incorrect use of CyclicModelPipelineCreator: model has no loop');
  }

  /** Return computations pipeline */
  public getPipeline(inputs: Float64Array): Pipeline {
    if (this.ivp.loop === null)
      throw new Error('Incorrect use of CyclicModelPipelineCreator: model has no loop');

    const argSize = argName2IdxMap.size;

    const wrappers: Wrapper[] = [];

    const restricted = inputs.length - 1;
    const repetitions = inputs[restricted];

    if (repetitions < LOOP_MIN_COUNT)
      throw new Error(`Incorrect repetitions count: ${repetitions}, expected: > 0`);

    for (let i = 0; i < repetitions; ++i) {
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

      preProcLines.push('const inputs = new Float64Array(arguments[1]);');

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

      preProcLines.push(`return inputs.slice(0, ${restricted});`);

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

      wrappers.push({
        preproc: preProcLines.join('\n'),
        out: null,
        postproc: postProcLines.join('\n'),
      });
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
