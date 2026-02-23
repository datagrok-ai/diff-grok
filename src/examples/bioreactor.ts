/** Example. Solve the following initial value problem

      dx / dt = x + y - t
      dy / dt = x * y + t
      x(0) = 1
      y(0) = -1

    on [0, 2] with the step 0.01.
 */

import {ODEs, lsoda, mrt} from '../../index';

// constants
const VLinit = 7.2;
const VTV = 10;
const speed = 400;
const diam = 6;
const power = 2.1;
const pH = 7.4;
const k1red = 0.05604;
const k1ox = 0.0108;
const k2Fd = 1.35;
const k2Fa = 110400000;
const k2Kd = 0.04038;
const k2Ka = 120000000;
const k3FKa = 181200000;
const k3FKd = 0.01188;
const k4ox = 0.0108;
const k4red = 0.05604;
const ktox = 0.005;
const krcyst = 0;
const pO2sat = 100;
const pKa2MEA = 8.18;
const H = 1.072069378;
const R = 0.082;

const _t0 = 0;
const _t1 = 1000;
const _h = 1;
const FFox = 0.2;
const KKox = 0.2;
const FFred = 0.1;
const KKred = 0.1;
const Ffree = 0;
const Kfree = 0;
const FKred = 0;
const FKox = 0;
const MEAthiol = 15;
const CO2 = 0.12;
const yO2P = 0.209;
const CYST = 0;
const VL = 7.2;
const qin = 1;
const yO2in = 0.21;
const T = 300;
const P = 1;
const switchTime = 135;

// the problem definition
const task = {
  name: 'Bioreactor',
  arg: {name: 't', start: _t0, finish: _t1, step: _h},
  initial: [FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, FKox, MEAthiol, CO2, yO2P, CYST, VL],
  func: (t: number, _y: Float64Array, _output: Float64Array) => {
    // extract function values
    const FFox = _y[0];
    const KKox = _y[1];
    const FFred = _y[2];
    const KKred = _y[3];
    const Ffree = _y[4];
    const Kfree = _y[5];
    const FKred = _y[6];
    const FKox = _y[7];
    const MEAthiol = _y[8];
    const CO2 = _y[9];
    const yO2P = _y[10];
    const CYST = _y[11];
    const VL = _y[12];

    // evaluate expressions
    const KF = pow(VL, -0.65) * 0.065 * pow(speed**3 * diam**5 * power / 2.16e12, 0.361);
    const MA = MEAthiol * pow(10, pH - pKa2MEA);
    const qout = qin - KF * (yO2P * H - CO2) * VL * R * T / (P * 1000);
    const OTR = KF * (yO2P * H - CO2);
    const Fin = t < switchTime ? 0 : 0.025;
    const Fper = t < switchTime ? 0.025 : Fin;
    const E0 = VLinit / VL;
    const E1 = (MA * E0)**2;
    const E11 = k1red * FFox * E0 * E1;
    const E12 = k1ox * FFred * E0;
    const E31 = k2Fd * FFred * E0;
    const E32 = k2Fa * (Ffree * E0)**2 * E1;
    const E21 = k1red * KKox * E0 * E1;
    const E22 = k1ox * KKred * E0;
    const E41 = k2Kd * KKred * E0;
    const E42 = k2Ka * (Kfree * E0)**2 * E1;
    const E51 = k3FKa * Ffree * E0 * Kfree * E0;
    const E52 = k3FKd * FKred * E0;
    const E61 = k4ox * FKred * E0 * (CYST * E0)**2;
    const E62 = k4red * FKox * E0 * E1;
    const E70 = E0 * CO2;
    const E71 = (MEAthiol * E0)**2;
    const E72 = (E70 >= 0) ? sqrt(E70) : 0;

    // compute output
    _output[0] = -E11 + E12;
    _output[1] = -E21 + E22;
    _output[2] = E11 - E12 - E31 + E32;
    _output[3] = E21 - E22 - E41 + E42;
    _output[4] = 2 * (E31 - E32) - E51 + E52;
    _output[5] = 2 * (E41 - E42) - E51 + E52;
    _output[6] = E51 - E52 - E61 + E62;
    _output[7] = E61 - E62;
    _output[8] = 2 * (-E11 + E12 - E21 + E22 + E31 + E41 - E32 - E42 - E62 - ktox * E71 * E72) -
      (MEAthiol + MA) * (Fin + Fper) / VL;
    _output[9] = (Fin * pO2sat * 0.0022 - 2 * Fper * CO2) / VL + OTR - 0.5 * ktox * E71 * E72;
    _output[10] = -OTR * (VL / (VTV - VL)) * R * T * P + yO2in * qin - yO2P * qout;
    _output[11] = ktox * E71 * E72 - krcyst * CYST * E0 - (Fin + Fper) * CYST / VL;
    _output[12] = Fin - Fper;
  },
  tolerance: 0.00005,
  solutionColNames: ['FFox', 'KKox', 'FFred', 'KKred', 'Ffree', 'Kfree', 'FKred',
    'FKox', 'MEAthiol', 'CO2', 'yO2P', 'CYST', 'VL'],
};

// used Math-functions
const pow = (x: number, y: number) => Math.pow(x, y);
const tan = (x: number) => Math.tan(x);
const sqrt = (x: number) => Math.sqrt(x);
const exp = (x: number) => Math.exp(x);

const solution = lsoda(task);
