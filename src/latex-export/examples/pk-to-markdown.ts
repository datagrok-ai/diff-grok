/** Example: PK model (with loop) → Markdown with default options. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: PK
#description: Pharmacokinetic (PK) simulation: one-compartment model
#equations:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2

#expressions:
  C2 = centr / V2

#loop:
  count = 1 {category: Dosing; min: 1; max: 10}
  depot += dose

#argument: t
  start = 0 {units: h; caption: begin; category: Dosing; min: 0; max: 1}
  final = 12 {units: h; caption: end; category: Dosing; min: 5; max: 15}
  step = 0.01 {units: h; caption: step; category: Dosing; min: 0.01; max: 0.1}

#inits:
  depot = 0 {category: Initial values}
  centr = 0 {category: Initial values}

#parameters:
  dose = 1E4 {category: Dosing; min: 1E3; max: 2E4; step: 1E3}
  KA = 0.3 {caption: rate constant; category: Parameters; min: 0.1; max: 1}
  CL = 2 {caption: clearance; category: Parameters; min: 1; max: 5}
  V2 = 4 {caption: central volume; category: Parameters; min: 1; max: 10}

#tolerance: 1E-9`;

// Compact Markdown — loop model
const markdown = convertIvpToLatex(MODEL, {format: 'markdown', compact: true});
console.log(markdown);
