/** Example: nimotuzumab model → LaTeX, no inits, no cdot. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Nimotuzumab
#description: Nimotuzumab disposition model
#comment:
  Source: https://www.mdpi.com/1999-4923/12/12/1147
#equations:
  dA1/dt = (-(CL * A3 / V1 + Q / V1) * A1 + Q / V2 * A2 - kint * Rtot * A1 / (Kss + A1 / V1))
           / (1 + Rtot * Kss / (Kss + A1 / V1)**2)

  dA2/dt = (Q / V1 * A1 - Q / V2 * A2) / (1 + Rtotp * Kss / (Kss + A2 / V2)**2)

  dA3/dt = Kin * (1 + Smax * A1**gamma / (S50**gamma + A1**gamma)) - Kout * A3

#argument: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10}
  finish = 70 {caption: Final; category: Time; min: 50; max: 100}
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 1}

#inits:
  A1 = 0.43 {category: Initials; min: 0.2; max: 0.6}
  A2 = 0.55 {category: Initials; min: 0.3; max: 0.8}
  A3 = 10.32 {category: Initials; min: 5; max: 15}

#output:
  t {caption: Time, h}
  A1 {caption: Central}
  A2 {caption: Periferal}

#parameters:
  CL = 9.64E-3 {category: Parameters; min: 0.005; max: 0.015}
  V1 = 2.63 {category: Parameters; min: 0.7; max: 3}
  Q = 2.88E-2 {category: Parameters; min: 0.015; max: 0.035}
  V2 = 9.92E-3 {category: Parameters; min: 0.0067; max: 0.0265}
  Kss = 15.5 {category: Parameters; min: 10; max: 45}
  kint = 1E-2 {category: Parameters; min: 0.01; max: 0.05}
  Rtot = 1.05E-2 {category: Parameters; min: 0.007; max: 0.032}
  Rtotp = 956 {category: Parameters; min: 120; max: 3500}
  Kin = 1.33E-2 {category: Parameters; min: 0.01; max: 0.05}
  Kout = 1.33E-2 {category: Parameters; min: 0.001; max: 0.02}
  S50 = 8.57 {category: Parameters; min: 5; max: 12}
  Smax = 3.18 {category: Parameters; min: 1; max: 5}
  gamma = 0.5 {category: Parameters; min: 0.1; max: 1}

#tolerance: 1E-7`;

// Compact LaTeX, suppress initial conditions, juxtaposition
const latex = convertIvpToLatex(MODEL, {
  compact: true,
  includeInits: false,
  useCdot: false,
});
console.log(latex);
