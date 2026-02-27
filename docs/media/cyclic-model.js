/* eslint-disable max-len */
/** This example shows how to apply pipelines and cyclic models.
    This approach can be used for in-webworkers analysis of models.

    Here, we consider pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model.
*/
import * as DGL from '../../index';
/** 1. Model specification */
const model = `#name: PK-PD
#tags: model
#description: Pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model
#equations:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
  d(peri)/dt  = Q * C2 - Q * C3
  d(eff)/dt  = Kin - Kout * (1 - C2/(EC50 + C2)) * eff

#expressions:
  C2 = centr / V2
  C3 = peri / V3

#loop:
  _count = 10 {caption: count; category: Dosing; min: 1; max: 20} [Number of doses]
  depot += dose

#argument: t
  _t0 = 0 {units: h; caption: begin; category: Dosing; min: 0; max: 1} [Begin of dosing interval]
  _t1 = 12 {units: h; caption: end; category: Dosing; min: 5; max: 15} [End of dosing interval]
  _h = 1 {units: h; caption: step; category: Dosing; min: 0.01; max: 0.1} [Time step of simulation]  

#inits:  
  depot = 0 {category: Initial values}
  centr = 0 {category: Initial values} [Central]
  peri = 0 {category: Initial values} [Peripheral]
  eff = 0.2 {category: Initial values} [Effective compartment rate]

#parameters:  
  dose = 1e4 {category: Dosing; min: 1e3; max: 2e4; step: 1e3} [Dosage]
  KA = 0.3 {caption: rate constant; category: Parameters; min: 0.1; max: 1}
  CL = 2 {caption: clearance; category: Parameters; min: 1; max: 5}
  V2 = 4 {caption: central volume; category: Parameters; min: 1; max: 10} [Central compartment volume]
  Q = 1 {caption: inter rate; category: Parameters; min: 0.1; max: 1} [Intercompartmental rate]
  V3 = 30 {caption: peri volume; category: Parameters; min: 20; max: 40} [Peripheral compartment volume]
  EC50 = 8 {caption: effect; category: Parameters; min: 1; max: 10}
  Kin = 0.2 {caption: Kin; category: Parameters; min: 0.1; max: 0.5} [The first-order production constant]
  Kout = 0.2 {caption: Kout; category: Parameters; min: 0.1; max: 0.5} [The first-order dissipation rate constant]
  
#tolerance: 1e-9

#meta.solver: {method: 'rkdp'; maxTimeMs: 5000}`;
/** 2. Generate IVP-objects: for the main thread & for computations in webworkers */
const ivp = DGL.getIVP(model);
const ivpWW = DGL.getIvp2WebWorker(ivp);
/** 3. Perform computations */
try {
    // 3.1) Extract names of outputs
    const outputNames = DGL.getOutputNames(ivp);
    const outSize = outputNames.length;
    // 3.2) Set model inputs
    const inputs = {
        _count: 10,
        _t0: 0,
        _t1: 16,
        _h: 1,
        depot: 0,
        centr: 0,
        peri: 0,
        eff: 0.2,
        dose: 10000,
        KA: 0.63,
        CL: 3.2,
        V2: 6.58,
        Q: 0.622,
        V3: 35.6,
        EC50: 5.41,
        Kin: 0.272,
        Kout: 0.276,
    };
    const inputVector = DGL.getInputVector(inputs, ivp);
    // 3.3) Create a pipeline
    const creator = DGL.getPipelineCreator(ivp);
    const pipeline = creator.getPipeline(inputVector);
    // 3.4) Apply pipeline to perform computations
    const solution = DGL.applyPipeline(pipeline, ivpWW, inputVector);
    // 3.5) Print results
    // 3.5.1) Table header
    let line = '         ';
    outputNames.forEach((name) => line += name + '           ');
    console.log(line);
    // 3.5.2) Table with solution
    const length = solution[0].length;
    for (let i = 0; i < length; ++i) {
        line = '';
        for (let j = 0; j < outSize; ++j)
            line += solution[j][i].toFixed(8) + '     ';
        console.log(line);
    }
}
catch (err) {
    console.log('Simulation failed: ', err instanceof Error ? err.message : 'Unknown problem!');
}
//# sourceMappingURL=cyclic-model.js.map