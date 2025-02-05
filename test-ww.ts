import {getIVP, solveIvp} from './index';
import {getIvp2WebWorker} from './index';

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

const ivpWW = getIvp2WebWorker(ivp);

console.log(ivpWW);

const inputs = new Float64Array([0, 1, 0.1, 2, 0, 10, -10]);
const paramsCount = 2;

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
