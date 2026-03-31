/** Example: Robertson model → LaTeX without metadata (name, description, comment). */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Robertson
#description: Robertson's chemical reaction model
#comment: This is classic example of stiff ODEs.
#equations:
  dA/dt = -0.04 * A + 1e4 * B * C
  dB/dt = 0.04 * A - 1e4 * B * C - 3e7 * B**2
  dC/dt = 3e7 * B**2

#inits:
  A = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  B = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  C = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}

#argument: t
  start = 0 {units: sec; caption: Initial; category: Time; min: 0; max: 1}
  finish = 40 {units: sec; caption: Final; category: Time; min: 2; max: 50}
  step = 0.01 {units: sec; caption: Step; category: Time; min: 0.01; max: 0.1}

#tolerance: 1e-7`;

// Compact, no metadata
const latex = convertIvpToLatex(MODEL, {compact: true, includeMetadata: false});
console.log(latex);
