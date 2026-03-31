/** Example: energy-n-control model → Markdown, no constants. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Energy-n-Control
#description: The energy-n-control demo model
#comment:
  This is an artificial example. It illustrates the use of the Math module
  and creating custom functions using JavaScript.

#equations:
  dx/dt = E1 * y * weight
  dy/dt = E2 * x * weight

#expressions:
  E1 = P2 * exp(-t)
  E2 = P2 * cos(t)
  weight = sin(PI * P1)

  // this makes further code shorter
  ceil = Math.ceil

  // simple if-then-else custom function
  func1 = (p, t) => (p < 8) ? ceil(t) : ceil(t**2 / 20)

  // another custom function using the previous one
  func2 = (p, t) => (p < 5) ? ceil(exp(-t) * 10) : func1(p, t);

  // usage of custom function
  control = (P1 < 3) ? ceil(cos(10 * t)) : func2(P1, t)

  energy = x**2 + y**2

#parameters:
  P1 = 1.1 {caption: frequency; category: Parameters; min: 1; max: 10}
  P2 = -1.1 {caption: shift; category: Parameters; min: -10; max: -1}

#inits:
  x = 2 {category: Initial state; min: 0; max: 5}
  y = 0 {category: Initial state; min: -2; max: 2}

#output:
  t {caption: time}
  energy
  control

#argument: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10}
  finish = 10 {caption: Final; category: Time; min: 10; max: 20}
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001}

#tolerance: 5e-5`;

// Compact Markdown, suppress constants table
const markdown = convertIvpToLatex(MODEL, {
  format: 'markdown',
  compact: true,
  includeConstants: false,
});
console.log(markdown);
