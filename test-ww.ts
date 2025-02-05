import {getIVP, solveIvp} from './index';
import {getIvp2WebWorker} from './index';

/*
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

const inputs = new Float64Array([0, 10, 0.1, 2, 0, 10, -10]);
const ivp = getIVP(model);
const ivpWW = getIvp2WebWorker(ivp);

try {
    // Solve the problem  
    const solution = solveIvp(ivpWW, inputs, ivpWW.paramNames.length);
  
    // Output results
    console.log(ivp.arg.name, '    ', ivp.deqs.solutionNames[0], '  ', ivp.deqs.solutionNames[1]);
  
    const length = solution[0].length;
  
    for (let i = 0; i < length; ++i)
      console.log(solution[0][i], '    ', solution[1][i], '  ', solution[2][i]);
  } catch (err) {
    console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
  }
*/

/*
const model = `#name: Template 
#equations:
  dy/dt = -y + sin(C * t) / t

#argument: t
  initial = 0.01 {min: 0.01; max: 10}
    final = 15   {min: 15; max: 150}
    step = 0.01  {min: 0.001; max: 0.1}

#constants:
  C = 2

#inits:  
  y = 0 {min: 0; max: 9}`;

const ivp = getIVP(model);

const inputs = new Float64Array([0.001, 15, 0.01, 0]);

const ivpWW = getIvp2WebWorker(ivp);

console.log(ivpWW);



try {
    // Solve the problem  
    const solution = solveIvp(ivpWW, inputs, ivpWW.paramNames.length);
  
    // Output results
    console.log(ivp.arg.name, '    ', ivp.deqs.solutionNames[0]);
  
    const length = solution[0].length;
  
    for (let i = 0; i < length; ++i)
      console.log(solution[0][i], '    ', solution[1][i]);
  } catch (err) {
    console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
  }*/

const model =  `#name: Chem react
#tags: model
#description: Mass-action kinetics illustration
#comment: 
  Source: https://doi.org/10.1002/ijch.201800003.
#equations:
  dx1/dt = -k1 * x1 + k2 * (x2)**2 + k3 * x1 * x3 
           - k4 * (x1)**2 - 2 * k5 * (x1)**2 + k6 * x2 * x4

  dx2/dt = 2 * k1 * x1 - 2 * k2 * (x2)**2 
           + k5 * (x1)**2 - k6 * x2 * x4

  dx3/dt = -k3 * x1 * x3 + k4 * (x1)**2 + k6 * x2 * x4

  dx4/dt = k5 * (x1)**2 - k6 * x2 * x4

#parameters:
  k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}
  k2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}
  k3 = 1.2 {category: Reaction parameters; min: 0.1; max: 5}
  k4 = 3.4 {category: Reaction parameters; min: 0.1; max: 5}
  k5 = 2.3 {category: Reaction parameters; min: 0.1; max: 5}
  k6 = 4.5 {category: Reaction parameters; min: 0.1; max: 5}

#inits:
  x1 = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x2 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x3 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x4 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}

#argument: t
  initial = 0 {units: min; caption: Initial; category: Time; min: 0; max: 5} [Initial time of simulation]
  final = 6 {units: min; caption: Final; category: Time; min: 6; max: 10} [Final time of simulation]
  step = 0.01 {units: min; caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]

#tolerance: 5e-5

#meta.solver: {maxTimeMs: 10}`;

const ivp = getIVP(model);

const inputs = new Float64Array([0, 6, 0.01, 1, 0, 0, 0, 0.7, 0.9, 1.2, 3.4, 2.3, 4.5]);

const ivpWW = getIvp2WebWorker(ivp);

console.log(ivpWW);

try {
    // Solve the problem  
    const solution = solveIvp(ivpWW, inputs);
  
    // Output results
    console.log(ivp.arg.name, '    ', ivp.deqs.solutionNames[0], '    ', ivp.deqs.solutionNames[1], '    ', ivp.deqs.solutionNames[2], '    ', ivp.deqs.solutionNames[3]);
  
    const length = solution[0].length;
  
    for (let i = 0; i < length; ++i)
      console.log(solution[0][i], '    ', solution[1][i], '    ', solution[2][i], '    ', solution[3][i], '    ', solution[4][i]);
  } catch (err) {
    console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
  }