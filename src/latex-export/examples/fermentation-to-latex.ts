/** Example: fermentation model → LaTeX, no parameters table. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Fermentation
#description: Simulation of fermentation process in the ethanol production
#comment:
  Source: https://core.ac.uk/download/pdf/11737483.pdf.
#equations:
  dP/dt = r * X
  dS/dt = -q * X
  dX/dt = V * S / (K + S) * X

#argument: t
  initial = 0 {units: d; caption: Initial; category: Time; min: 0; max: 10}
  final = 70 {units: d; caption: Final; category: Time; min: 15; max: 100}
  step = 0.01 {units: d; caption: Step; category: Time; min: 0.01; max: 1}

#inits:
  P = 4.276 {units: ml/ml; caption: ethanol, P; category: Initials; min: 3; max: 6}
  S = 0.3185 {units: mg/ml; caption: glucose, S; category: Initials; min: 0.1; max: 0.5}
  X = 9.2E-2 {units: mg/ml; caption: saccharomyces, X; category: Initials; min: 0.01; max: 0.2}

#parameters:
  r = 1.3455 {units: mg/ml; category: Parameters; min: 1; max: 2}
  q = 1.1129E-2 {units: mg/ml; category: Parameters; min: 0.01; max: 0.2}
  V = 8.7E-2 {units: mg/ml; category: Parameters; min: 0.01; max: 0.2}
  K = 6.28E-2 {category: Parameters}

#tolerance: 1E-7`;

// Compact LaTeX, suppress parameters
const latex = convertIvpToLatex(MODEL, {compact: true, includeParameters: false});
console.log(latex);
