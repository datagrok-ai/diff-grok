// Features for creating output code for a pipeline

import {IVP, MATH_CONSTS, POW_IDX, MATH_FUNCS} from '../scripting-tools/scripting-tools';
import {ARG_INP_COUNT} from './constants';

const SHIFT = 1;

/** Return output for the model WITHOUT expressions */
function getOutputNoExprs(ivp: IVP): string {
  const outputs = ivp.outputs;
  if (outputs === null)
    throw new Error('Grok Lib issue: no outputs in the initial value problem');

  if (exprsInOutputs(ivp))
    throw new Error('Grok Lib issue: no expressions are allowed');


  const lines: string[] = ['const solution = [];'];

  outputs.forEach((_, name) => {
    if (name === ivp.arg.name)
      lines.push('solution.push(arguments[0][0]);');
    else {
      const idx = ivp.deqs.solutionNames.indexOf(name);

      if (idx > -1)
        lines.push(`solution.push(arguments[0][${idx + 1}]);`);
    }
  });

  lines.push('return solution;');

  return lines.join('\n');
} // getOutputNoExprs

/** Return output for the model WITH expressions */
function getOutputWithExprs(ivp: IVP): string {
  const outputs = ivp.outputs;
  if (outputs === null)
    throw new Error('Grok Lib issue: no outputs in the initial value problem');

  const exprs = ivp.exprs;

  if (!exprsInOutputs(ivp) || (exprs === null))
    throw new Error('Grok Lib issue: expressions must be given');


  const argName = ivp.arg.name;
  const lines: string[] = ['const solution = [];'];

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

  // Model parameters
  if (ivp.params !== null) {
    let idx = ARG_INP_COUNT + ivp.deqs.solutionNames.length;
    ivp.params.forEach((_, name) => {
      lines.push(`const ${name} = arguments[1][${idx}];`);
      ++idx;
    });
  }

  lines.push('');

  lines.push('const length = arguments[0][0].length;');
  const funcNames = ivp.deqs.solutionNames;

  lines.push(`let ${argName} = 0;`);
  funcNames.forEach((name) => lines.push(`let ${name} = 0;`));

  lines.push('');

  exprs.forEach((_, name) => {
    lines.push(`let ${name} = 0;`);

    if (outputs.has(name))
      lines.push(`let ${name}Raw = new Float64Array(length);`);
  });

  lines.push('');

  lines.push('for (let k = 0; k < length; ++k) {');
  lines.push(`  ${argName} = arguments[0][0][k];`);
  funcNames.forEach((name, idx) => lines.push(`  ${name} = arguments[0][${SHIFT + idx}][k];`));
  lines.push('');
  exprs.forEach((expr, name) => {
    lines.push(`  ${name} = ${expr};`);

    if (outputs.has(name))
      lines.push(`  ${name}Raw[k] = ${name};`);
  });
  lines.push('}');

  lines.push('');

  outputs.forEach((_, name) => {
    if (name === argName)
      lines.push('solution.push(arguments[0][0]);');
    else {
      const idx = ivp.deqs.solutionNames.indexOf(name);

      if (idx > -1)
        lines.push(`solution.push(arguments[0][${idx + 1}]);`);
      else if (exprs.has(name))
        lines.push(`solution.push(${name}Raw);`);
    }
  });

  lines.push('return solution;');

  return lines.join('\n');
} // getOutputWithExprs

/** Return a code for output extraction
 * @internal
 */
export function getOutputCode(ivp: IVP): string | null {
  const outputs = ivp.outputs;
  if (outputs === null)
    return null;

  checkApplicability(ivp);

  if (exprsInOutputs(ivp))
    return getOutputWithExprs(ivp);
  else
    return getOutputNoExprs(ivp);
} // getOutputCode

/** Check whether outputs have items from expressions */
function exprsInOutputs(ivp: IVP): boolean {
  const outputs = ivp.outputs;
  if (outputs === null)
    return false;

  const exprs = ivp.exprs;
  if (exprs === null)
    return false;

  let flag = false;

  exprs.forEach((_, name) => flag ||= outputs.has(name));

  return flag;
}

/** Check applicability of model outputs */
function checkApplicability(ivp: IVP): void {
  const outputs = ivp.outputs;

  if (outputs === null)
    throw new Error('Model has no outputs');

  const exprs = ivp.exprs;
  if (exprs === null)
    return;

  if (ivp.loop !== null) {
    if (exprsInOutputs(ivp))
      throw new Error('Non-supported model: expressions in output & loops');

    return;
  }

  if (ivp.updates !== null) {
    if (exprsInOutputs(ivp))
      throw new Error('Non-supported model: expressions in output and updates');
  }
} // checkApplicability

/** Return names of the model outputs
 * @param ivp - initial value problem
 * @returns an array of names of the model outputs
 */
export function getOutputNames(ivp: IVP): string[] {
  const outputs = ivp.outputs;
  if (outputs === null)
    return [ivp.arg.name].concat(ivp.deqs.solutionNames);

  const names = new Array<string>(outputs.size);

  let idx = 0;
  outputs.forEach((output) => {
    names[idx] = output.caption;
    ++idx;
  });

  return names;
} // getOutputNames
