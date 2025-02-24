/* eslint-disable max-len */
/** This example shows how to apply pipelines and models with customized outputs.
    This approach can be used for in-webworkers analysis of models.

    Here, we consider modeling queues: https://en.wikipedia.org/wiki/M/M/c_queue
*/

import * as DGL from '../../index';

/** 1. Model specification: the M|M|2|2 model  */
const model = `#name: M|M|2|2
#description: Modelling the system with two servers & two waiting places (EXPLAINED)
#equations:  
  dp0/dt = -la * p0 + mu * p1
  dp1/dt = la * p0 - (la + mu) * p1 + 2 * mu * p2
  dp2/dt = la * p1 - (la + 2 * mu) * p2 + 2 * mu * p3
  dp3/dt = la * p2 - (la + 2 * mu) * p3 + 2 * mu * p4

#expressions:
  la = 1 / arrival
  mu = 1 / service
  p4 = 1 - p0 - p1 - p2 - p3
  busy = 1 - p0 - p1
  queue = p3 + 2 * p4

#argument: t
  initial = 0     {min: 0;      max: 10;   caption: start;   category: Time; units: min}  [Initial time of simulation]
    final = 60    {min: 20;     max: 100;  caption: finish;  category: Time; units: min}  [Final time of simulation]
     step = 1     {min: 0.01;   max: 0.1;  caption: step;    category: Time; units: min}  [Time step of simulation] 

#inits:
  p0 = 1   {min: 0; max: 1; category: Initial state; caption: empty}   [Probability that initially there are NO customers]
  p1 = 0   {min: 0; max: 1; category: Initial state; caption: single}  [Probability that initially there is ONE customer]
  p2 = 0   {min: 0; max: 1; category: Initial state; caption: two}     [Probability that initially there are TWO customers]
  p3 = 0   {min: 0; max: 1; category: Initial state; caption: three}   [Probability that initially there are THREE customers]

#output:
  t     {caption: Time, min}
  p0    {caption: P(Empty)}
  busy  {caption: P(Full load)}
  queue {caption: E(Queue)}

#parameters:
  arrival = 10   {min: 1; max: 100; category: Means; caption: arrival;  units: min} [Mean arrival time]
  service = 100  {min: 1; max: 100; category: Means; caption: service; units: min} [Mean service time]`;

/** 2. Generate IVP-objects: for the main thread & for computations in webworkers */
const ivp = DGL.getIVP(model);
const ivpWW = DGL.getIvp2WebWorker(ivp);

/** 3. Create pipeline */
const pipeline: DGL.Pipeline = {
  wrappers: [ // define as many wrappers as you need
    {
      preproc: null, // define a code that preprocess inputs
      out: null, // specify a code that extracts customized outputs
      postproc: null, // define a code that performs postprocessing of inputs
    },
  ],
  out: DGL.getOutputCode(ivp), // final output of computations
};

/** 4. Perform computations */
try {
  // 4.1) Extract names of outputs
  const outputNames = DGL.getOutputNames(ivp);
  const outSize = outputNames.length;

  // 4.2) Set model inputs
  const inputs = new Float64Array([
    0, // Initial time of simulation
    60, // Final time of simulation
    1, // Time step of simulation
    1, // Probability that initially there are NO customers
    0, // Probability that initially there is ONE customer
    0, // Probability that initially there are TWO customers
    0, // Probability that initially there are THREE customers
    10, // Mean arrival time
    100, // Mean service time
  ]);

  // 4.3) Apply pipeline to perform computations
  const solution = DGL.applyPipeline(pipeline, ivpWW, inputs);

  // 4.4) Print results

  // 4.4.1) Table header
  let line = '';
  outputNames.forEach((name) => line += name + '      ');
  console.log(line);

  // 4.4.2) Table with solution
  const length = solution[0].length;
  for (let i = 0; i < length; ++i) {
    line = '';

    for (let j = 0; j < outSize; ++j)
      line += solution[j][i].toFixed(8) + '     ';

    console.log(line);
  }
} catch (err) {
  console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}

