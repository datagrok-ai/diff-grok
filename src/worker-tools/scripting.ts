import {IVP, MATH_CONSTS, POW_IDX, MATH_FUNCS} from '../scripting-tools/scripting-tools';

import {SolverOptions} from "../solver-tools/solver-defs";

export type IVP2WebWorker = {
    name: string,
    equations: {rightHandParts: string[], names: string[]},
    expressions: {rightHandParts: string[], names: string[]},
    arg: {name: string, vals: number[]},
    initVals: {names: string[], vals: number[]},
    consts: {names: string[], vals: number[]},
    params: {names: string[], vals: number[]},
    tolerance: number,
    solverOpts: SolverOptions,
    usedMathFuncs: number[],
    usedMathConsts: number[],
    funcMainBody: string,
};

const DEFAULT_METHOD = 'ros34prw';
const DEFAULT_MAX_TIME = -1;

function getSolverOpts(opts: string): SolverOptions {
    try {
      const json = JSON.parse(opts
        .replaceAll(`'`, `"`)
        .replace('method','"method"')
        .replace('maxTimeMs', '"maxTimeMs"')
        .replace('maxIterations', '"maxIterations"')
        .replace(';', ',')
      );
      
      return json;
    } catch(e) {
      return {
        method: DEFAULT_METHOD,
        maxTimeMs: DEFAULT_MAX_TIME,
        maxIterations: 1,
      };
    }
}

function getFuncMainBody(ivp: IVP2WebWorker): string {
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
    ivp.consts.names.forEach((name, idx) => lines.push(`const ${name} = ${ivp.consts.vals[idx]};`));    

    // extract arg & function values
    lines.push(`const ${ivp.arg.name} = arguments[0];`);
    ivp.equations.names.forEach((name, idx) => lines.push(`const ${name} = arguments[1][${idx}];`))

    // evaluate expressions
    ivp.expressions.names.forEach(((name, idx) => lines.push(`const ${name} = ${ivp.expressions.rightHandParts[idx]};`)));

    // compute output
    ivp.equations.rightHandParts.forEach(((expr, idx) => lines.push(`arguments[2][${idx}] = ${expr};`)));

    return lines.join('\n');
}

export function getIvp2WebWorker(ivp: IVP): IVP2WebWorker {
    
    const arg = ivp.arg;
    let names: string[] = [];
    let vals: number[] = [];
    let exprs: string[] = [];

    // equations
    names = [];
    exprs = [];

    ivp.deqs.equations.forEach((val, key) => {
        names.push(key);
        exprs.push(val);
    });

    const equations = {rightHandParts: exprs, names: names};

    // expressions
    names = [];
    exprs = [];

    if (ivp.exprs !== null) {
        ivp.exprs.forEach((val, key) => {
            names.push(key);
            exprs.push(val);
        });
    }

    const expressions = {rightHandParts: exprs, names: names};

    // initial values
    names = [];
    vals = [];

    ivp.deqs.solutionNames.forEach((name) => {
        names.push(name);
        vals.push(ivp.inits.get(name)?.value || 0);
    });

    const initVals = {names: names, vals: vals};

    // costants
    names = [];
    vals = [];

    if (ivp.consts !== null) {
        ivp.consts.forEach((val, key) => {
            names.push(key);
            vals.push(val.value);
      });
    }
    const consts = {names: names, vals: vals};

    // parameters
    names = [];
    vals = [];
    
    if (ivp.params !== null) {
        ivp.params.forEach((val, key) => {
            names.push(key);
            vals.push(val.value);
      });
    }
    const params = {names: names, vals: vals};

    const ivpWW: IVP2WebWorker = {
        name: ivp.name,
        equations: equations,
        expressions: expressions,
        initVals: initVals,
        arg: {name: arg.name, vals: [arg.initial.value, arg.final.value, arg.step.value]},
        consts: consts,
        params: params,
        tolerance: Number(ivp.tolerance),
        solverOpts: getSolverOpts(ivp.solverSettings),
        usedMathConsts: ivp.usedMathConsts,
        usedMathFuncs: ivp.usedMathFuncs,
        funcMainBody: '',
    };

    ivpWW.funcMainBody = getFuncMainBody(ivpWW);

    return ivpWW;
}

export function getFunc(ivp: IVP2WebWorker, paramVals: Float64Array): string {
    const lines: string[] = [];

    ivp.params.names.forEach((name, idx) => lines.push(`const ${name} = ${paramVals[idx]};`));
  
    return lines.concat(ivp.funcMainBody).join('\n');
}