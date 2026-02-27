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

// CVODE TypeScript port - Initial step size estimation
// Ported from SUNDIALS CVODE cvode.c (cvHin, cvUpperBoundH0, cvYddNorm)

import { CvodeMem, wrmsNorm, vLinearSum, vScale, CV_SUCCESS,
         CV_RHSFUNC_FAIL, RHSFUNC_RECVR, CV_TOO_CLOSE,
         HLB_FACTOR, HUB_FACTOR, H_BIAS, MAX_ITERS, TINY } from './common';

/**
 * cvHin - Compute initial step size h0.
 *
 * The algorithm:
 * 1. Compute a bound hub on |h| from ydot and tdist
 * 2. Set hlb = HLB_FACTOR * tround  (lower bound)
 * 3. Iterate to estimate ||y''|| and compute h candidates
 * 4. h0 = H_BIAS * geometric mean of hlb and hub candidates
 */
export function cvHin(mem: CvodeMem, tout: number): number {
  const N = mem.cv_N;
  const sign = (tout - mem.cv_tn >= 0.0) ? 1.0 : -1.0;
  const tdist = Math.abs(tout - mem.cv_tn);
  const tround = Math.max(Math.abs(mem.cv_tn), Math.abs(tout)) * 2.2204460492503131e-16;

  if (tdist < 2.0 * tround) return CV_TOO_CLOSE;

  // Set lower and upper bounds on h0, initial guess
  const hlb = HLB_FACTOR * tround;
  let hub = cvUpperBoundH0(mem, tdist);

  // Initial guess - geometric mean
  let hg = Math.sqrt(hlb * hub);
  if (hub < 100.0 * hlb) {
    hg = 0.5 * hub;
  }

  // Iterate to estimate ydd norm
  let hnew = hg;
  let retval: number;

  for (let count = 0; count < MAX_ITERS; count++) {
    hg = hnew;

    // Compute y''norm via finite difference
    const result = cvYddNorm(mem, hg * sign);
    retval = result.retval;
    if (retval !== CV_SUCCESS) return retval;

    const yddnrm = result.yddnrm;

    // If yddnrm is zero, use current hg
    if (yddnrm * hub * hub > 2.0) {
      hnew = Math.sqrt(2.0 / yddnrm);
    } else {
      hnew = Math.sqrt(hg * hub);
    }

    // Refine: make hnew be a geometric mean biased toward smaller
    hnew = H_BIAS * hnew + (1.0 - H_BIAS) * hg;

    // Check convergence
    if (hnew / hg > 2.0 || hnew / hg < 0.5) continue;
    // If converged, break
    break;
  }

  // Bound hnew
  hnew = Math.min(hnew, hub);
  hnew = Math.max(hnew, hlb);

  // Apply bounds on step size
  if (mem.cv_hmax_inv > 0.0) {
    hnew = Math.min(hnew, 1.0 / mem.cv_hmax_inv);
  }
  if (mem.cv_hmin > 0.0) {
    hnew = Math.max(hnew, mem.cv_hmin);
  }

  // Bound by tdist
  hnew = Math.min(hnew, tdist);

  mem.cv_h = sign * hnew;
  mem.cv_next_h = mem.cv_h;
  mem.cv_hscale = mem.cv_h;

  // Scale zn[1] = h * y'(t0)
  vScale(mem.cv_h, mem.cv_zn[1], mem.cv_zn[1], N);

  return CV_SUCCESS;
}

/**
 * cvUpperBoundH0 - Compute upper bound on initial step size
 * Based on ||y'|| weighted norm
 */
function cvUpperBoundH0(mem: CvodeMem, tdist: number): number {
  const N = mem.cv_N;

  // Get WRMS norm of y' (zn[1] = f(t0, y0) before scaling)
  const hub_inv = wrmsNorm(N, mem.cv_zn[1], mem.cv_ewt);

  // Bound based on tdist
  let hub = HUB_FACTOR * tdist;

  // Use smaller of two bounds
  if (hub * hub_inv > 1.0) {
    hub = 1.0 / hub_inv;
  }

  return hub;
}

/**
 * cvYddNorm - Estimate ||y''|| via finite difference
 * y(t0 + hg) ~ y0 + hg*y'0
 * y'(t0 + hg) = f(t0 + hg, y(t0 + hg))
 * y'' ~ (y'(t0+hg) - y'(t0)) / hg
 */
function cvYddNorm(mem: CvodeMem, hg: number): { retval: number; yddnrm: number } {
  const N = mem.cv_N;

  // y = zn[0] + hg * zn[1] (zn[1] still holds unscaled f at this point)
  vLinearSum(hg, mem.cv_zn[1], 1.0, mem.cv_zn[0], mem.cv_y, N);

  // Evaluate f(t0 + hg, y)
  const retval = mem.cv_f!(mem.cv_tn + hg, mem.cv_y, mem.cv_tempv, mem.cv_user_data);
  mem.cv_nfe++;
  if (retval < 0) return { retval: CV_RHSFUNC_FAIL, yddnrm: 0 };
  if (retval > 0) return { retval: RHSFUNC_RECVR, yddnrm: 0 };

  // tempv = (f(t0+hg, y) - zn[1]) / hg  (estimate of y'')
  const hg_inv = 1.0 / hg;
  vLinearSum(hg_inv, mem.cv_tempv, -hg_inv, mem.cv_zn[1], mem.cv_tempv, N);

  const yddnrm = wrmsNorm(N, mem.cv_tempv, mem.cv_ewt);

  return { retval: CV_SUCCESS, yddnrm };
}
