import {IVP} from '../scripting-tools';
import {MATH_CONSTS, MATH_FUNCS, POW_IDX} from '../scripting-tools/scripting-tools';
import {getOutputCode} from './output';
import {Pipeline, Wrapper} from './pipeline';
import {PipelineCreator} from './pipeline-creator';

import {argName2IdxMap, ARG, ARG_COL_IDX, CORRECTION_FACTOR, DEFAULT_CORRECTION} from './constants';

/** Pipeline creator for model with updates */
export class UpdatesModelPipelineCreator extends PipelineCreator {
  constructor(ivp: IVP) {
    super(ivp);

    if (ivp.updates === null)
      throw new Error('Incorrect use of UpdatesModelPipelineCreator: model has no updates');
  }

  /** Return computations pipeline */
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
