/* eslint-disable max-len */
import {mrt, ros3prw, ros34prw} from '../../index';

export const methods = new Map([
  ['MRT', mrt],
  ['ROS3PRw', ros3prw],
  ['ROS34PRw', ros34prw],
]);

export const MAX_MAD = 0.1;
export const TIMEOUT_MS = 10000;

const BASIC_MODEL = `#name: M|M|2|2
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
  _t0 = 0     {min: 0;      max: 10;   caption: start;   category: Time; units: min}  [Initial time of simulation]
  _t1 = 60    {min: 20;     max: 100;  caption: finish;  category: Time; units: min}  [Final time of simulation]
  _h  = 1     {min: 0.01;   max: 0.1;  caption: step;    category: Time; units: min}  [Time step of simulation] 

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

const BASIC_MODEL_INPUTS = {
  _t0: 0,
  _t1: 60,
  _h: 1,
  p0: 1,
  p1: 0,
  p2: 0,
  p3: 0,
  arrival: 10,
  service: 100,
};

const BASIC_MODEL_OUTPUTS = 4;

const CYCLIC_MODEL = `#name: PK-PD
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
  
#tolerance: 1e-9`;

const CYCLIC_MODEL_INPUTS = {
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

const CYCLIC_MODEL_OUTPUTS = 5;

const MULTISTAGE_MODEL = `#name: GA-production
#tags: model
#description: Gluconic acid (GA) production by Aspergillus niger modeling
#equations:
  dX/dt = rX
  dS/dt = -gamma * rX - lambda * X
  dO/dt = Kla * (Cod - O) - delta * rX - phi * X
  dP/dt = alpha * rX + beta * X

#expressions:
  mu = muM * S / (Ks + S) * O / (Ko + O)
  rX = mu * X

#argument: t, 1-st stage
  _t0 = 0 {units: h; caption: initial; category: Misc} [Start of the process]
  _t1 = 60 {units: h; caption: 1-st stage; category: Durations; min: 20; max: 80} [Duration of the 1-st stage]
  step = 0.1 {units: h; caption: step; category: Misc; min: 0.01; max: 1} [Time step of simulation]

#update: 2-nd stage
  duration = overall - _t1
  S += 70

#inits:  
  X = 5 {units: kg/m³; caption: biomass; category: Initial concentrations; min: 1; max: 10} [Aspergillus niger biomass]
  S = 150 {units: kg/m³; caption: glucose; category: Initial concentrations; min: 50; max: 200} [Glucose]
  O = 7 {units: kg/m³; caption: oxygen; category: Initial concentrations; min: 1; max: 10} [Dissolved oxygen]
  P = 0 {units: kg/m³; caption: acid; category: Initial concentrations; min: 0; max: 0.1} [Gluconic acid]

#output:
  t {caption: time}
  X {caption: biomass}
  S {caption: glucose}
  O {caption: oxygen}
  P {caption: acid}

#parameters:
  overall = 100 {units: h; category: Durations; min: 100; max: 140} [Overall duration]
  muM = 0.668 {units: 1/h; category: Parameters} [Monod type model parameter]
  alpha = 2.92 {category: Parameters} [Monod type model parameter]
  beta = 0.131 {units: 1/h; category: Parameters} [Monod type model parameter]
  gamma = 2.12 {category: Parameters} [Monod type model parameter]
  lambda = 0.232 {units: 1/h; category: Parameters} [Monod type model parameter]
  delta = 0.278 {category: Parameters} [Monod type model parameter]
  phi = 4.87e-3 {units: 1/h; category: Parameters} [Monod type model parameter]
  Ks = 1.309e2 {units: g/L; category: Parameters} [Monod type model parameter]
  Ko = 3.63e-4 {units: g/L; category: Parameters} [Monod type model parameter]
  Kla = 1.7e-2 {units: 1/s; category: Parameters} [Volumetric mass transfer coefficient]
  Cod = 15 {units: kg/m³; category: Parameters} [Liquid phase dissolved oxygen saturation concentration]
  
#tolerance: 1e-9`;

const MULTISTAGE_MODEL_INPUTS = {
  _t0: 0,
  _t1: 60,
  _h: 1,
  X: 5,
  S: 150,
  O: 7,
  P: 0,
  overall: 100,
  muM: 0.668,
  alpha: 2.92,
  beta: 0.131,
  gamma: 2.12,
  lambda: 0.232,
  delta: 0.278,
  phi: 4.87e-3,
  Ks: 1.309e2,
  Ko: 3.63e-4,
  Kla: 1.7e-2,
  Cod: 15,
};

const MULTISTAGE_MODEL_OUTPUTS = 5;

type TestProblem = {
  model: string,
  inputs: Record<string, number>,
  outputsCount: number,
};

export const problems = new Map<string, TestProblem>([
  ['basic', {model: BASIC_MODEL, inputs: BASIC_MODEL_INPUTS, outputsCount: BASIC_MODEL_OUTPUTS}],
  ['cyclic', {model: CYCLIC_MODEL, inputs: CYCLIC_MODEL_INPUTS, outputsCount: CYCLIC_MODEL_OUTPUTS}],
  ['multistage', {model: MULTISTAGE_MODEL, inputs: MULTISTAGE_MODEL_INPUTS, outputsCount: MULTISTAGE_MODEL_OUTPUTS}],
]);
