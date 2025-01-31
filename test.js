"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var index_1 = require("./index");
var model = "#name: Extended\n#tags: model\n#description: 2D ordinary differential equations system sample\n#comment:\n  This is an extended template. It has additional scripting annotations.\n\n#equations:\n  dx/dt = E1 * y + sin(t)\n  dy/dt = E2 * x - pow(t, 5)\n\n#expressions:\n  E1 = C1 * exp(-t) + P1\n  E2 = C2 * cos(2 * t) + P2\n\n#constants:\n  C1 = 1\n  C2 = 3\n\n#parameters:\n  P1 = 1 {category: Parameters; min: 1; max: 10} [P1 parameter]\n  P2 = -1 {category: Parameters; min: -10; max: -1} [P2 parameter]\n\n#inits:  \n  y = 0 {category: Initial values; min: -2; max: 2} [Initial value of y]\n  x = 2 {category: Initial values; min: 0; max: 5} [Initial value of x]\n\n#argument: t\n  start = 0 {caption: Initial; category: Time; min: 0; max: 10} [Initial time of simulation]\n  finish = 10 {caption: Final; category: Time; min: 10; max: 20} [Final time of simulation]\n  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]\n\n#tolerance: 5e-5\n\n#meta.solver: {method: 'ros34prw'}";
var ivp = (0, index_1.getIVP)(model);
console.log(ivp);
var DEFAULT_METHOD = 'ros34prw';
var DEFAULT_MAX_TIME = -1;
function getMethodOpts(opts) {
    try {
        return JSON.parse(opts
            .replaceAll("'", "\"")
            .replace('method', '"method"')
            .replace('maxTimeMs', '"maxTimeMs"')
            .replace('maxIterations', '"maxIterations"'));
    }
    catch (e) {
        return {
            method: DEFAULT_METHOD,
            maxTimeMs: DEFAULT_MAX_TIME,
        };
    }
}
function getIvp2WebWorker(ivp) {
    var _a, _b;
    var arg = ivp.arg;
    var names = [];
    var vals = [];
    var exprs = [];
    // equations
    names = [];
    exprs = [];
    ivp.deqs.equations.forEach(function (val, key) {
        names.push(key);
        exprs.push(val);
    });
    var equations = { rightHandParts: exprs, names: names };
    // expressions
    names = [];
    exprs = [];
    if (ivp.exprs !== null) {
        ivp.exprs.forEach(function (val, key) {
            names.push(key);
            exprs.push(val);
        });
    }
    var expressions = { rightHandParts: exprs, names: names };
    // initial values
    names = [];
    vals = [];
    ivp.deqs.solutionNames.forEach(function (name) {
        var _a;
        names.push(name);
        vals.push(((_a = ivp.inits.get(name)) === null || _a === void 0 ? void 0 : _a.value) || 0);
    });
    var initVals = { names: names, vals: vals };
    // costants
    names = [];
    vals = [];
    if (ivp.consts !== null) {
        ivp.consts.forEach(function (val, key) {
            names.push(key);
            vals.push(val.value);
        });
    }
    var consts = { names: names, vals: vals };
    // parameters
    names = [];
    vals = [];
    if (ivp.params !== null) {
        ivp.params.forEach(function (val, key) {
            names.push(key);
            vals.push(val.value);
        });
    }
    var params = { names: names, vals: vals };
    // solver settings
    var solverOpts = getMethodOpts(ivp.solverSettings);
    return {
        name: ivp.name,
        equations: equations,
        expressions: expressions,
        initVals: initVals,
        arg: { name: arg.name, vals: [arg.initial.value, arg.final.value, arg.step.value] },
        consts: consts,
        params: params,
        tolerance: Number(ivp.tolerance),
        method: (_a = solverOpts.method) !== null && _a !== void 0 ? _a : DEFAULT_METHOD,
        maxTimeMs: (_b = solverOpts.maxTimeMs) !== null && _b !== void 0 ? _b : DEFAULT_MAX_TIME,
        usedMathConsts: ivp.usedMathConsts,
        usedMathFuncs: ivp.usedMathFuncs,
    };
}
console.log('=========================================================================================');
var ivpWW = getIvp2WebWorker(ivp);
console.log('=========================================================================================');
var MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh'];
var POW_IDX = MATH_FUNCS.indexOf('pow');
var MATH_CONSTS = ['PI', 'E'];
function getFunc4worker(ivp, paramVals) {
    var lines = [];
    // Used math funcs
    ivp.usedMathFuncs.forEach(function (idx) {
        if (idx !== POW_IDX)
            lines.push("const ".concat(MATH_FUNCS[idx], " = (x) => Math.").concat(MATH_FUNCS[idx], "(x);"));
        else
            lines.push("const pow = (x, y) => Math.pow(x, y);");
    });
    // Used math consts
    ivp.usedMathConsts.forEach(function (idx) { return lines.push("const ".concat(MATH_CONSTS[idx], " = Math.").concat(MATH_CONSTS[idx], ";")); });
    // Model constants
    ivp.consts.names.forEach(function (name, idx) { return lines.push("const ".concat(name, " = ").concat(ivp.consts.vals[idx], ";")); });
    // Model parameters
    ivp.params.names.forEach(function (name, idx) { return lines.push("const ".concat(name, " = ").concat(paramVals[idx], ";")); });
    // extract arg & function values
    lines.push("const ".concat(ivp.arg.name, " = arguments[0];"));
    ivp.equations.names.forEach(function (name, idx) { return lines.push("const ".concat(name, " = arguments[1][").concat(idx, "];")); });
    // evaluate expressions
    ivp.expressions.names.forEach((function (name, idx) { return lines.push("const ".concat(name, " = ").concat(ivp.expressions.rightHandParts[idx], ";")); }));
    // compute output
    ivp.equations.rightHandParts.forEach((function (expr, idx) { return lines.push("arguments[2][".concat(idx, "] = ").concat(expr, ";")); }));
    return lines.join('\n');
}
function getMethod(name) {
    switch (name) {
        case 'mrt':
            return index_1.mrt;
        case 'ros3prw':
            return index_1.ros3prw;
        default:
            return index_1.ros34prw;
    }
}
var funcCode = getFunc4worker(ivpWW, new Float64Array([1, -1]));
var func = new Function(funcCode);
console.log(funcCode);
var task = {
    name: 'Extended',
    arg: { name: 't', start: 0, finish: 10, step: 0.01 },
    initial: [2, 0],
    func: func,
    tolerance: 5e-5,
    solutionColNames: ['x', 'y']
};
try {
    // Solve the problem
    var solution = getMethod(ivpWW.method)(task);
    // Output results
    console.log(task.arg.name, '    ', task.solutionColNames[0], '  ', task.solutionColNames[1]);
    var length_1 = solution[0].length;
    for (var i = 0; i < length_1; ++i)
        console.log(solution[0][i], '    ', solution[1][i], '  ', solution[2][i]);
}
catch (err) {
    console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}
