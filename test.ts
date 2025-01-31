import {getIVP, IVP, mrt, ros34prw, ros3prw} from './index';

const model = `#name: Extended
#tags: model
#description: 2D ordinary differential equations system sample
#comment:
  This is an extended template. It has additional scripting annotations.

#equations:
  dx/dt = E1 * y + sin(t)
  dy/dt = E2 * x - pow(t, 5)

#expressions:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2

#constants:
  C1 = 1
  C2 = 3

#parameters:
  P1 = 1 {category: Parameters; min: 1; max: 10} [P1 parameter]
  P2 = -1 {category: Parameters; min: -10; max: -1} [P2 parameter]

#inits:  
  y = 0 {category: Initial values; min: -2; max: 2} [Initial value of y]
  x = 2 {category: Initial values; min: 0; max: 5} [Initial value of x]

#argument: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10} [Initial time of simulation]
  finish = 10 {caption: Final; category: Time; min: 10; max: 20} [Final time of simulation]
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]

#tolerance: 5e-5

#meta.solver: {method: 'ros34prw'}`;

const ivp = getIVP(model);

console.log(ivp);

type IVP2WebWorker = {
    name: string,
    equations: {rightHandParts: string[], names: string[]},
    expressions: {rightHandParts: string[], names: string[]},
    arg: {name: string, vals: number[]},
    initVals: {names: string[], vals: number[]},
    consts: {names: string[], vals: number[]},
    params: {names: string[], vals: number[]},
    tolerance: number,
    method: string,
    maxTimeMs: number,
    usedMathFuncs: number[],
    usedMathConsts: number[],    
};

const DEFAULT_METHOD = 'ros34prw';
const DEFAULT_MAX_TIME = -1;

function getMethodOpts(opts: string) {
    try {
        return JSON.parse(opts
            .replaceAll(`'`, `"`)
            .replace('method','"method"')
            .replace('maxTimeMs', '"maxTimeMs"')
            .replace('maxIterations', '"maxIterations"')
        );
    } catch(e) {
      return {
        method: DEFAULT_METHOD,
        maxTimeMs: DEFAULT_MAX_TIME,
      };
    }
}

function getIvp2WebWorker(ivp: IVP): IVP2WebWorker {
    
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

    // solver settings
    const solverOpts = getMethodOpts(ivp.solverSettings);

    return {
        name: ivp.name,
        equations: equations,
        expressions: expressions,
        initVals: initVals,
        arg: {name: arg.name, vals: [arg.initial.value, arg.final.value, arg.step.value]},
        consts: consts,
        params: params,
        tolerance: Number(ivp.tolerance),
        method: solverOpts.method ?? DEFAULT_METHOD,
        maxTimeMs: solverOpts.maxTimeMs ?? DEFAULT_MAX_TIME,
        usedMathConsts: ivp.usedMathConsts,
        usedMathFuncs: ivp.usedMathFuncs,
    };
}

console.log('=========================================================================================');

const ivpWW = getIvp2WebWorker(ivp);

console.log('=========================================================================================');

const MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh'];
const POW_IDX = MATH_FUNCS.indexOf('pow');
const MATH_CONSTS = ['PI', 'E'];

function getFunc4worker(ivp: IVP2WebWorker, paramVals: Float64Array): string {
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

    // Model parameters
    ivp.params.names.forEach((name, idx) => lines.push(`const ${name} = ${paramVals[idx]};`));

    // extract arg & function values
    lines.push(`const ${ivp.arg.name} = arguments[0];`);
    ivp.equations.names.forEach((name, idx) => lines.push(`const ${name} = arguments[1][${idx}];`))

    // evaluate expressions
    ivp.expressions.names.forEach(((name, idx) => lines.push(`const ${name} = ${ivp.expressions.rightHandParts[idx]};`)));

    // compute output
    ivp.equations.rightHandParts.forEach(((expr, idx) => lines.push(`arguments[2][${idx}] = ${expr};`)));

    return lines.join('\n');
}

function getMethod(name: string) {
    switch (name) {
        case 'mrt':
            return mrt;
        
        case 'ros3prw':
            return ros3prw;
    
        default:
            return ros34prw;
    }
}

const funcCode = getFunc4worker(ivpWW, new Float64Array([1, -1]));

type Func = (t: number, y: Float64Array, output: Float64Array) => void;

const func = new Function(funcCode) as Func;

console.log(funcCode);

let task = {
    name: 'Extended',
    arg: {name: 't', start: 0, finish: 10, step: 0.01},
    initial: [2, 0],
    func: func,
    tolerance: 5e-5,
    solutionColNames: ['x', 'y']
};

try {
  // Solve the problem
  const solution = getMethod(ivpWW.method)(task);

  // Output results
  console.log(task.arg.name, '    ', task.solutionColNames[0], '  ', task.solutionColNames[1]);

  const length = solution[0].length;

  for (let i = 0; i < length; ++i)
    console.log(solution[0][i], '    ', solution[1][i], '  ', solution[2][i]);
} catch (err) {
  console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}

function solveIvpWW(func: Func, inputs: Float64Array, paramsCount: number, method, callback) {
    
}
