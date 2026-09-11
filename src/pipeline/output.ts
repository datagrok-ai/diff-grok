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

/** Output-expression names, in #expressions order — the columns a stage appends. */
function outputExprNames(ivp: IVP): string[] {
  const outputs = ivp.outputs;
  const exprs = ivp.exprs;
  if (outputs === null || exprs === null)
    return [];

  const names: string[] = [];
  exprs.forEach((_, name) => {
    if (outputs.has(name))
      names.push(name);
  });

  return names;
} // outputExprNames

/** Lines that recompute the expressions row-by-row from the argument & solution columns of
    arguments[0], reading parameters from the input vector arguments[1], and fill a Float64Array
    <name>Raw for each output expression. Shared by the basic recompute and the per-stage code. */
function getExprComputeLines(ivp: IVP): string[] {
  const outputs = ivp.outputs!;
  const exprs = ivp.exprs!;
  const argName = ivp.arg.name;
  const funcNames = ivp.deqs.solutionNames;
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

  // Model parameters (per-stage input vector)
  if (ivp.params !== null) {
    let idx = ARG_INP_COUNT + funcNames.length;
    ivp.params.forEach((_, name) => {
      lines.push(`const ${name} = arguments[1][${idx}];`);
      ++idx;
    });
  }

  lines.push('');
  // __len / __i are underscore-prefixed so they cannot collide with a model constant or parameter
  // (e.g. a model with a constant named `k` would otherwise be shadowed by the loop counter).
  lines.push('const __len = arguments[0][0].length;');

  lines.push(`let ${argName} = 0;`);
  funcNames.forEach((name) => lines.push(`let ${name} = 0;`));

  lines.push('');

  exprs.forEach((_, name) => {
    lines.push(`let ${name} = 0;`);

    if (outputs.has(name))
      lines.push(`let ${name}Raw = new Float64Array(__len);`);
  });

  lines.push('');

  lines.push('for (let __i = 0; __i < __len; ++__i) {');
  lines.push(`  ${argName} = arguments[0][0][__i];`);
  funcNames.forEach((name, idx) => lines.push(`  ${name} = arguments[0][${SHIFT + idx}][__i];`));
  lines.push('');
  exprs.forEach((expr, name) => {
    lines.push(`  ${name} = ${expr};`);

    if (outputs.has(name))
      lines.push(`  ${name}Raw[__i] = ${name};`);
  });
  lines.push('}');

  lines.push('');

  return lines;
} // getExprComputeLines

/** Return output for a basic (non-cyclic, non-multistage) model WITH expressions.
    Recomputes the output expressions once, at top level, over the full solution. */
function getOutputWithExprs(ivp: IVP): string {
  const outputs = ivp.outputs;
  const exprs = ivp.exprs;

  if (outputs === null || exprs === null || !exprsInOutputs(ivp))
    throw new Error('Grok Lib issue: expressions must be given');

  const argName = ivp.arg.name;
  const lines: string[] = getExprComputeLines(ivp);

  lines.push('const solution = [];');

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

/** Return per-stage output code for cyclic/multistage models: recompute the output expressions
    for this stage (using this stage's parameters), then APPEND them as columns to the stage
    solution. Returns null when the model has no expressions in its outputs. */
export function getStageOutputCode(ivp: IVP): string | null {
  if (!exprsInOutputs(ivp))
    return null;

  const lines: string[] = getExprComputeLines(ivp);

  // Keep [arg, f0..fN] intact and append the output-expression columns after them.
  lines.push('const solution = arguments[0];');
  outputExprNames(ivp).forEach((name) => lines.push(`solution.push(${name}Raw);`));
  lines.push('return solution;');

  return lines.join('\n');
} // getStageOutputCode

/** Return top-level output code for cyclic/multistage models with output expressions: the stages
    already appended the expression columns, so this only selects columns into #output order. */
function getOutputSelectWithExprCols(ivp: IVP): string {
  const outputs = ivp.outputs!;
  const argName = ivp.arg.name;
  const funcNames = ivp.deqs.solutionNames;
  const exprNames = outputExprNames(ivp);
  const lines: string[] = ['const solution = [];'];

  outputs.forEach((_, name) => {
    if (name === argName) {
      lines.push('solution.push(arguments[0][0]);');
      return;
    }

    const funcIdx = funcNames.indexOf(name);
    if (funcIdx > -1) {
      lines.push(`solution.push(arguments[0][${funcIdx + 1}]);`);
      return;
    }

    const exprIdx = exprNames.indexOf(name);
    if (exprIdx > -1)
      lines.push(`solution.push(arguments[0][${1 + funcNames.length + exprIdx}]);`);
  });

  lines.push('return solution;');

  return lines.join('\n');
} // getOutputSelectWithExprCols

/** Return a code for output extraction
 * @internal
 */
export function getOutputCode(ivp: IVP): string | null {
  const outputs = ivp.outputs;
  if (outputs === null)
    return null;

  if (!exprsInOutputs(ivp))
    return getOutputNoExprs(ivp);

  // Cyclic (#loop) & multistage (#update) models compute expression columns per stage; the
  // top level only selects them. Basic models recompute the expressions here.
  if (ivp.loop !== null || ivp.updates !== null)
    return getOutputSelectWithExprCols(ivp);

  return getOutputWithExprs(ivp);
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
