// Scripting tools for solving IVPs in web-workers

import {IVP, MATH_CONSTS, POW_IDX, MATH_FUNCS} from '../scripting-tools/scripting-tools';
import {SolverOptions} from '../solver-tools/solver-defs';

/** Initial value problem structure to be passed to web-worker */
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

/** Return IVP-structure to be passed to web-worker */
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
