/** Example: GA production model (multi-stage) → Markdown, no cdot. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Acid Production
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
  start = 0 {units: h; caption: initial; category: Misc}
  stage1 = 60 {units: h; caption: 1-st stage; category: Durations; min: 20; max: 80}
  step = 0.1 {units: h; caption: step; category: Misc; min: 0.01; max: 1}

#update: 2-nd stage
  duration = overall - _t1
  S += 70

#inits:
  X = 5 {units: kg/m³; caption: biomass; category: Initial concentrations; min: 1; max: 10}
  S = 150 {units: kg/m³; caption: glucose; category: Initial concentrations; min: 50; max: 200}
  O = 7 {units: kg/m³; caption: oxygen; category: Initial concentrations; min: 1; max: 10}
  P = 0 {units: kg/m³; caption: acid; category: Initial concentrations; min: 0; max: 0.1}

#output:
  t {caption: time}
  X {caption: biomass}
  S {caption: glucose}
  O {caption: oxygen}
  P {caption: acid}

#parameters:
  overall = 100 {units: h; category: Durations; min: 100; max: 140}
  muM = 0.668 {units: 1/h; category: Parameters}
  alpha = 2.92 {category: Parameters}
  beta = 0.131 {units: 1/h; category: Parameters}
  gamma = 2.12 {category: Parameters}
  lambda = 0.232 {units: 1/h; category: Parameters}
  delta = 0.278 {category: Parameters}
  phi = 4.87e-3 {units: 1/h; category: Parameters}
  Ks = 1.309e2 {units: g/L; category: Parameters}
  Ko = 3.63e-4 {units: g/L; category: Parameters}
  Kla = 1.7e-2 {units: 1/s; category: Parameters}
  Cod = 15 {units: kg/m³; category: Parameters}

#tolerance: 1e-9`;

// Markdown with juxtaposition — good for multi-stage models
const markdown = convertIvpToLatex(MODEL, {
  format: 'markdown',
  useCdot: false,
});
console.log(markdown);
