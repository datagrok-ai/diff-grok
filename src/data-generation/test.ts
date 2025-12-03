// Here, we consider modeling queues: https://en.wikipedia.org/wiki/M/M/c_queue
import * as DGL from '../../index';

// 1. Model specification: the M|M|2|2 model
const model = `#name: Chem react
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

#output:
  t {caption: Time}
  x1
  x4 {caption: Conc}

#parameters:
  k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5; caption: param1}
  k2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}
  k3 = 1.2 {category: Reaction parameters; min: 0.1; max: 5; caption: fixed param1}
  k4 = 3.4 {category: Reaction parameters; min: 0.1; max: 5; caption: fixed param2}
  k5 = 2.3 {category: Reaction parameters; min: 0.1; max: 5}
  k6 = 4.5

#inits:
  x1 = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x2 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x3 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x4 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}

#argument: t
  initial = 0 {units: min; caption: Initial; category: Time} [Initial time of simulation]
  final = 6 {units: min; caption: Final; category: Time; min: 6; max: 10} [Final time of simulation]
  step = 0.01 {units: min; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simulation]

#tolerance: 5e-5`;

const gen = new DGL.SyntheticDataGenerator(model);

//console.log(gen);

// 3. Perform computations
// try {
//   // 3.1) Extract names of outputs
//   const outputNames = DGL.getOutputNames(ivp);
//   const outSize = outputNames.length;

//   // 3.2) Set model inputs
//   const inputs = {
//     _t0: 0, // Initial time of simulation
//     _t1: 60, // Final time of simulation
//     _h: 1, // Time step of simulation
//     p0: 1, // Probability that initially there are NO customers
//     p1: 0, // Probability that initially there is ONE customer
//     p2: 0, // Probability that initially there are TWO customers
//     p3: 0, // Probability that initially there are THREE customers
//     arrival: 10, // Mean arrival time
//     service: 100, // Mean service time
//   };
//   const inputVector = DGL.getInputVector(inputs, ivp);

//   // 3.3) Create a pipeline
//   const creator = DGL.getPipelineCreator(ivp);
//   const pipeline = creator.getPipeline(inputVector);

//   // 3.4) Apply pipeline to perform computations
//   const solution = DGL.applyPipeline(pipeline, ivpWW, inputVector);
//   // 3.5) Print results

//   // 3.5.1) Table header
//   let line = '';
//   outputNames.forEach((name) => line += name + '      ');
//   console.log(line);

//   // 3.5.2) Table with solution
//   const length = solution[0].length;
//   for (let i = 0; i < length; ++i) {
//     line = '';

//     for (let j = 0; j < outSize; ++j)
//       line += solution[j][i].toFixed(8) + '     ';
//     console.log(line);
//   }
// } catch (err) {
//   console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
// }
