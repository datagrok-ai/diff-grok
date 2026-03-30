/** Example: Extended model → Markdown with all sections. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Extended
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
  P1 = 1 {category: Parameters; min: 1; max: 10}
  P2 = -1 {category: Parameters; min: -10; max: -1}

#inits:
  x = 2 {category: Initial values; min: 0; max: 5}
  y = 0 {category: Initial values; min: -2; max: 2}

#argument: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10}
  finish = 10 {caption: Final; category: Time; min: 10; max: 20}
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001}

#tolerance: 5e-5`;

// Compact Markdown output
const markdown = convertIvpToLatex(MODEL, {format: 'markdown', compact: true});
console.log(markdown);
