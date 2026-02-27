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

import {CvodeMem, L_MAX} from './common';

/**
 * cvAltSum - Alternating sum: sum_{i=0}^{iend} (-1)^i * a[i] / (i + k)
 */
function cvAltSum(iend: number, a: Float64Array, k: number): number {
  if (iend < 0) return 0.0;
  let sum = 0.0;
  let sign = 1;
  for (let i = 0; i <= iend; i++) {
    sum += sign * (a[i] / (i + k));
    sign = -sign;
  }
  return sum;
}

/**
 * cvAdamsStart - Generate product polynomial coefficients for Adams method
 * Returns hsum.
 */
function cvAdamsStart(mem: CvodeMem, m: Float64Array): number {
  let hsum = mem.cv_h;
  m[0] = 1.0;
  for (let i = 1; i <= mem.cv_q; i++) m[i] = 0.0;

  for (let j = 1; j < mem.cv_q; j++) {
    if (j === mem.cv_q - 1 && mem.cv_qwait === 1) {
      const sum = cvAltSum(mem.cv_q - 2, m, 2);
      mem.cv_tq[1] = mem.cv_q * sum / m[mem.cv_q - 2];
    }
    const xi_inv = mem.cv_h / hsum;
    for (let i = j; i >= 1; i--)
      m[i] += m[i - 1] * xi_inv;

    hsum += mem.cv_tau[j];
  }
  return hsum;
}

/**
 * cvAdamsFinish - Complete Adams l and tq computation
 */
function cvAdamsFinish(mem: CvodeMem, m: Float64Array, M: Float64Array, hsum: number): void {
  const M0_inv = 1.0 / M[0];

  mem.cv_l[0] = 1.0;
  for (let i = 1; i <= mem.cv_q; i++)
    mem.cv_l[i] = M0_inv * (m[i - 1] / i);


  const xi = hsum / mem.cv_h;
  const xi_inv = 1.0 / xi;

  mem.cv_tq[2] = M[1] * M0_inv / xi;
  mem.cv_tq[5] = xi / mem.cv_l[mem.cv_q];

  if (mem.cv_qwait === 1) {
    for (let i = mem.cv_q; i >= 1; i--)
      m[i] += m[i - 1] * xi_inv;

    M[2] = cvAltSum(mem.cv_q, m, 2);
    mem.cv_tq[3] = M[2] * M0_inv / mem.cv_L;
  }

  mem.cv_tq[4] = mem.cv_nlscoef / mem.cv_tq[2];
}

/**
 * cvSetAdams - Compute Adams method coefficients l[] and tq[]
 */
export function cvSetAdams(mem: CvodeMem): void {
  const m = new Float64Array(L_MAX);
  const M = new Float64Array(3);

  if (mem.cv_q === 1) {
    mem.cv_l[0] = 1.0;
    mem.cv_l[1] = 1.0;
    mem.cv_tq[1] = 1.0;
    mem.cv_tq[5] = 1.0;
    mem.cv_tq[2] = 0.5;
    mem.cv_tq[3] = 1.0 / 12.0;
    mem.cv_tq[4] = mem.cv_nlscoef / mem.cv_tq[2];
    return;
  }

  const hsum = cvAdamsStart(mem, m);
  M[0] = cvAltSum(mem.cv_q - 1, m, 1);
  M[1] = cvAltSum(mem.cv_q - 1, m, 2);
  cvAdamsFinish(mem, m, M, hsum);
}

/**
 * cvAdjustAdams - Adjust history array for Adams order change
 */
export function cvAdjustAdams(mem: CvodeMem, deltaq: number): void {
  const N = mem.cv_N;

  if (deltaq === 1) {
    // On order increase, set new column of zn to zero
    mem.cv_zn[mem.cv_L].fill(0.0, 0, N);
    return;
  }

  // Order decrease: adjust zn[j] by multiples of zn[q]
  for (let i = 0; i <= mem.cv_qmax; i++) mem.cv_l[i] = 0.0;
  mem.cv_l[1] = 1.0;
  let hsum = 0.0;

  for (let j = 1; j <= mem.cv_q - 2; j++) {
    hsum += mem.cv_tau[j];
    const xi = hsum / mem.cv_hscale;
    for (let i = j + 1; i >= 1; i--)
      mem.cv_l[i] = mem.cv_l[i] * xi + mem.cv_l[i - 1];
  }

  for (let j = 1; j <= mem.cv_q - 2; j++)
    mem.cv_l[j + 1] = mem.cv_q * (mem.cv_l[j] / (j + 1));


  if (mem.cv_q > 2) {
    const znq = mem.cv_zn[mem.cv_q];
    for (let j = 2; j < mem.cv_q; j++) {
      const c = -mem.cv_l[j];
      const znj = mem.cv_zn[j];
      for (let i = 0; i < N; i++) znj[i] += c * znq[i];
    }
  }
}
