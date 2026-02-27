/**
 * TypeScript port of SUNDIALS CVODE v7.5.0 (https://github.com/LLNL/sundials)
 *
 * Copyright (c) 2025, Lawrence Livermore National Security, University of
 *   Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security and
 *   Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * Copyright (c) 2026 Datagrok.
 *
 * Licensed under the BSD-3-Clause License. See THIRD_PARTY_LICENSES
 * for the full license text and provenance chain of the original code.
 */

import {CvodeMem} from './common';

/**
 * cvSetBDF - Compute BDF method coefficients l[] and tq[]
 *
 * The polynomial l(x) = l_0 + l_1*x + ... + l_q*x^q is given by:
 *   Lambda(x) = (1 + x/xi*_q) * product_{i=1}^{q-1} (1 + x/xi_i)
 * where xi_i = (t_n - t_{n-i}) / h
 *
 * The tq array holds test quantities for error and convergence tests:
 *   tq[1] = coefficient for error at order q-1
 *   tq[2] = coefficient for error at order q
 *   tq[3] = coefficient for error at order q+1
 *   tq[4] = convergence test coefficient
 *   tq[5] = derivative for order q+2 vector
 */
export function cvSetBDF(mem: CvodeMem): void {
  let alpha0: number; let alpha0_hat: number; let xi_inv: number; let xistar_inv: number; let hsum: number;

  mem.cv_l[0] = 1.0;
  mem.cv_l[1] = 1.0;
  xi_inv = 1.0;
  xistar_inv = 1.0;
  for (let i = 2; i <= mem.cv_q; i++) mem.cv_l[i] = 0.0;

  alpha0 = -1.0;
  alpha0_hat = -1.0;
  hsum = mem.cv_h;

  if (mem.cv_q > 1) {
    for (let j = 2; j < mem.cv_q; j++) {
      hsum += mem.cv_tau[j - 1];
      xi_inv = mem.cv_h / hsum;
      alpha0 -= 1.0 / j;
      for (let i = j; i >= 1; i--)
        mem.cv_l[i] += mem.cv_l[i - 1] * xi_inv;
    }

    // j = q
    alpha0 -= 1.0 / mem.cv_q;
    xistar_inv = -mem.cv_l[1] - alpha0;
    hsum += mem.cv_tau[mem.cv_q - 1];
    xi_inv = mem.cv_h / hsum;
    alpha0_hat = -mem.cv_l[1] - xi_inv;

    for (let i = mem.cv_q; i >= 1; i--)
      mem.cv_l[i] += mem.cv_l[i - 1] * xistar_inv;
  }

  cvSetTqBDF(mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/**
 * cvSetTqBDF - Set test quantity array for BDF method
 */
function cvSetTqBDF(mem: CvodeMem, hsum: number, alpha0: number,
  alpha0_hat: number, xi_inv: number, xistar_inv: number): void {
  const A1 = 1.0 - alpha0_hat + alpha0;
  const A2 = 1.0 + mem.cv_q * A1;
  mem.cv_tq[2] = Math.abs(A1 / (alpha0 * A2));
  mem.cv_tq[5] = Math.abs(A2 * xistar_inv / (mem.cv_l[mem.cv_q] * xi_inv));

  if (mem.cv_qwait === 1) {
    if (mem.cv_q > 1) {
      const C = xistar_inv / mem.cv_l[mem.cv_q];
      const A3 = alpha0 + 1.0 / mem.cv_q;
      const A4 = alpha0_hat + xi_inv;
      const Cpinv = (1.0 - A4 + A3) / A3;
      mem.cv_tq[1] = Math.abs(C * Cpinv);
    } else
      mem.cv_tq[1] = 1.0;

    hsum += mem.cv_tau[mem.cv_q];
    const xi_inv2 = mem.cv_h / hsum;
    const A5 = alpha0 - 1.0 / (mem.cv_q + 1);
    const A6 = alpha0_hat - xi_inv2;
    const Cppinv = (1.0 - A6 + A5) / A2;
    mem.cv_tq[3] = Math.abs(Cppinv / (xi_inv2 * (mem.cv_q + 2) * A5));
  }
  mem.cv_tq[4] = mem.cv_nlscoef / mem.cv_tq[2];
}

/**
 * cvIncreaseBDF - Adjust history array for BDF order increase
 */
export function cvIncreaseBDF(mem: CvodeMem): void {
  const N = mem.cv_N;

  for (let i = 0; i <= mem.cv_qmax; i++) mem.cv_l[i] = 0.0;
  mem.cv_l[2] = 1.0;
  let alpha1 = 1.0;
  let prod = 1.0;
  let xiold = 1.0;
  let alpha0 = -1.0;
  let hsum = mem.cv_hscale;

  if (mem.cv_q > 1) {
    for (let j = 1; j < mem.cv_q; j++) {
      hsum += mem.cv_tau[j + 1];
      const xi = hsum / mem.cv_hscale;
      prod *= xi;
      alpha0 -= 1.0 / (j + 1);
      alpha1 += 1.0 / xi;
      for (let i = j + 2; i >= 2; i--)
        mem.cv_l[i] = mem.cv_l[i] * xiold + mem.cv_l[i - 1];

      xiold = xi;
    }
  }

  const A1 = (-alpha0 - alpha1) / prod;
  const znL = mem.cv_zn[mem.cv_L];
  const znSaved = mem.cv_zn[mem.cv_indx_acor];
  for (let i = 0; i < N; i++) znL[i] = A1 * znSaved[i];

  if (mem.cv_q > 1) {
    for (let j = 2; j <= mem.cv_q; j++) {
      const c = mem.cv_l[j];
      const znj = mem.cv_zn[j];
      for (let i = 0; i < N; i++) znj[i] += c * znL[i];
    }
  }
}

/**
 * cvDecreaseBDF - Adjust history array for BDF order decrease
 */
export function cvDecreaseBDF(mem: CvodeMem): void {
  const N = mem.cv_N;

  for (let i = 0; i <= mem.cv_qmax; i++) mem.cv_l[i] = 0.0;
  mem.cv_l[2] = 1.0;
  let hsum = 0.0;

  for (let j = 1; j <= mem.cv_q - 2; j++) {
    hsum += mem.cv_tau[j];
    const xi = hsum / mem.cv_hscale;
    for (let i = j + 2; i >= 2; i--)
      mem.cv_l[i] = mem.cv_l[i] * xi + mem.cv_l[i - 1];
  }

  if (mem.cv_q > 2) {
    const znq = mem.cv_zn[mem.cv_q];
    for (let j = 2; j < mem.cv_q; j++) {
      const c = -mem.cv_l[j];
      const znj = mem.cv_zn[j];
      for (let i = 0; i < N; i++) znj[i] += c * znq[i];
    }
  }
}
