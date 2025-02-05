import {getIVP, IVP, mrt, ros34prw, ros3prw, SolverOptions, getCallback, SolverMethod, ODEs} from './index';

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
  step = 0.1 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]

#tolerance: 5e-5

#meta.solver: {method: 'mrt'; maxTimeMs: 100}`;

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

      console.log(json);
      
      return json;
    } catch(e) {
      return {
        method: DEFAULT_METHOD,
        maxTimeMs: DEFAULT_MAX_TIME,
        maxIterations: 1,
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

    return {
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
}

console.log('=========================================================================================');

const ivpWW = getIvp2WebWorker(ivp);

console.log('=========================================================================================');

const MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh'];
const POW_IDX = MATH_FUNCS.indexOf('pow');
const MATH_CONSTS = ['PI', 'E'];

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

function getFunc(ivp: IVP2WebWorker, paramVals: Float64Array): string {
  const lines: string[] = [];
  // Model parameters
  ivp.params.names.forEach((name, idx) => lines.push(`const ${name} = ${paramVals[idx]};`));

  return lines.concat(ivp.funcMainBody).join('\n');
}

function getMethod(name: string): SolverMethod {
    switch (name) {
        case 'mrt':
            return mrt;
        
        case 'ros3prw':
            return ros3prw;
    
        default:
            return ros34prw;
    }
}

ivpWW.funcMainBody = getFuncMainBody(ivpWW);

type Func = (t: number, y: Float64Array, output: Float64Array) => void;

console.log(ivpWW);

/*try {
  // Solve the problem
  const method = getMethod(ivpWW.solverOpts.method);
  const callback = getCallback(ivpWW.solverOpts);
  const solution = method(task, callback);

  // Output results
  console.log(task.arg.name, '    ', task.solutionColNames[0], '  ', task.solutionColNames[1]);

  const length = solution[0].length;

  for (let i = 0; i < length; ++i)
    console.log(solution[0][i], '    ', solution[1][i], '  ', solution[2][i]);
} catch (err) {
  console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}*/

function checkInputs(ivp: IVP2WebWorker, inputs: Float64Array): void {
  // check inputs count
  const expected = ivp.arg.vals.length + ivp.initVals.vals.length + ivp.params.vals.length;
  if (expected !== inputs.length)
    throw new Error(`Incorrect inputs count, expected: ${expected}, current: ${inputs.length}`);
}

const inputs = new Float64Array([0, 10, 0.1, 2, 0, 10, -10]);
const paramsCount = 2;

checkInputs(ivpWW, inputs);

console.log(inputs.slice(inputs.length - paramsCount));

const ARG_INIT = 0;
const ARG_FINAL = 1;
const ARG_STEP = 2;

function solveIvp(ivp: IVP2WebWorker, inputs: Float64Array, paramsCount: number): Float64Array[] {
  checkInputs(ivp, inputs);
  const nonParamInputsCount = inputs.length - paramsCount;

  const funcCode = getFunc(ivp, inputs.slice(nonParamInputsCount));
  const func = new Function(funcCode) as Func;

  const odes: ODEs = {
    name: ivp.name,
    arg: {
      name: ivp.arg.name,
      start: inputs[ARG_INIT],
      finish: inputs[ARG_FINAL],
      step: inputs[ARG_STEP],
    },
    initial: inputs.slice(ARG_STEP + 1, nonParamInputsCount),
    func: func,
    tolerance: ivp.tolerance,
    solutionColNames: ivp.initVals.names,
  };

  const method = getMethod(ivpWW.solverOpts.method);
  const callback = getCallback(ivpWW.solverOpts);

  return method(odes, callback);
}

try {
  // Solve the problem  
  const solution = solveIvp(ivpWW, inputs, paramsCount);

  // Output results
  console.log(ivpWW.arg.name, '    ', ivpWW.initVals.names[0], '  ', ivpWW.initVals.names[1]);

  const length = solution[0].length;

  for (let i = 0; i < length; ++i)
    console.log(solution[0][i], '    ', solution[1][i], '  ', solution[2][i]);
} catch (err) {
  console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}
