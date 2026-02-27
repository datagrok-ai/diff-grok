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

// CVODE rootfinding routines
// Ported from SUNDIALS cvode.c v7.5.0 (rootfinding section, lines ~821-4737)

import {
  CvodeMem,
  CvodeRootFn,
  cvodeGetDky,
  vLinearSum,
  vScale,
  CV_SUCCESS,
  CV_RTFUNC_FAIL,
  RTFOUND,
  CLOSERT,
  UROUND,
  FUZZ_FACTOR,
  CV_NORMAL,
  CV_ONE_STEP,
} from './common';

// ============================================================================
// Local constants (matching SUNDIALS cvode.c)
// ============================================================================

const ZERO = 0.0;
const HALF = 0.5;
const ONE = 1.0;
const TWO = 2.0;
const FIVE = 5.0;
const HUNDRED = 100.0;
const PT1 = 0.1;

// ============================================================================
// Exported functions
// ============================================================================

/**
 * CVodeRootInit - Initialize rootfinding for CVODE.
 *
 * C reference: CVodeRootInit (cvode.c, lines ~821-1012)
 *
 * If nrtfn <= 0, disables rootfinding and returns CV_SUCCESS.
 * Otherwise, allocates arrays for root functions and sets defaults.
 *
 * @param mem   - CVODE memory block
 * @param nrtfn - number of root functions
 * @param g     - root function callback
 * @returns CV_SUCCESS on success
 */
export function cvodeRootInit(mem: CvodeMem, nrtfn: number, g: CvodeRootFn): number {
  // If nrtfn <= 0, disable rootfinding
  if (nrtfn <= 0) {
    mem.cv_nrtfn = 0;
    mem.cv_gfun = null;
    return CV_SUCCESS;
  }

  mem.cv_nrtfn = nrtfn;
  mem.cv_gfun = g;

  // Allocate Float64Arrays for g values
  mem.cv_glo = new Float64Array(nrtfn);
  mem.cv_ghi = new Float64Array(nrtfn);
  mem.cv_grout = new Float64Array(nrtfn);

  // Allocate Int32Arrays for root info and direction
  mem.cv_iroots = new Int32Array(nrtfn);
  mem.cv_rootdir = new Int32Array(nrtfn);

  // Allocate Uint8Array for active flags
  mem.cv_gactive = new Uint8Array(nrtfn);

  // Set default values: rootdir = 0 (both directions), gactive = 1 (all active)
  mem.cv_rootdir.fill(0);
  mem.cv_gactive.fill(1);

  return CV_SUCCESS;
}

/**
 * cvRcheck1 - First-step root check.
 *
 * C reference: cvRcheck1 (cvode.c, lines ~4251-4300)
 *
 * Completes initialization of rootfinding memory and checks whether g
 * has a zero both at and very near the initial point of the IVP.
 *
 * @param mem - CVODE memory block
 * @returns CV_SUCCESS or CV_RTFUNC_FAIL
 */
export function cvRcheck1(mem: CvodeMem): number {
  const nrtfn = mem.cv_nrtfn;
  const N = mem.cv_N;

  for (let i = 0; i < nrtfn; i++) {
    mem.cv_iroots[i] = 0;
  }

  mem.cv_tlo = mem.cv_tn;
  mem.cv_ttol = (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h)) * UROUND * HUNDRED;

  // Evaluate g at initial t and check for zero values
  const retval = mem.cv_gfun!(mem.cv_tlo, mem.cv_zn[0], mem.cv_glo, mem.cv_user_data);
  mem.cv_nge = 1;
  if (retval !== 0) {
    return CV_RTFUNC_FAIL;
  }

  let zroot = false;
  for (let i = 0; i < nrtfn; i++) {
    if (Math.abs(mem.cv_glo[i]) === ZERO) {
      zroot = true;
      mem.cv_gactive[i] = 0; // SUNFALSE
    }
  }
  if (!zroot) {
    return CV_SUCCESS;
  }

  // Some g_i is zero at t0; look at g at t0+(small increment)
  const hratio = Math.max(mem.cv_ttol / Math.abs(mem.cv_h), PT1);
  const smallh = hratio * mem.cv_h;
  const tplus = mem.cv_tlo + smallh;

  // y = zn[0] + hratio * zn[1]
  vLinearSum(ONE, mem.cv_zn[0], hratio, mem.cv_zn[1], mem.cv_y, N);

  const retval2 = mem.cv_gfun!(tplus, mem.cv_y, mem.cv_ghi, mem.cv_user_data);
  mem.cv_nge++;
  if (retval2 !== 0) {
    return CV_RTFUNC_FAIL;
  }

  // Re-activate components where g moved away from zero
  for (let i = 0; i < nrtfn; i++) {
    if (mem.cv_gactive[i] === 0 && Math.abs(mem.cv_ghi[i]) !== ZERO) {
      mem.cv_gactive[i] = 1; // SUNTRUE
      mem.cv_glo[i] = mem.cv_ghi[i];
    }
  }

  return CV_SUCCESS;
}

/**
 * cvRcheck2 - Subsequent-step root check for exact zeros at last root.
 *
 * C reference: cvRcheck2 (cvode.c, lines ~4323-4385)
 *
 * Checks for exact zeros of g at the last root found, for close pairs
 * of zeros (error condition), and for a new root at a nearby point.
 *
 * @param mem - CVODE memory block
 * @returns CV_SUCCESS, CLOSERT, RTFOUND, or CV_RTFUNC_FAIL
 */
export function cvRcheck2(mem: CvodeMem): number {
  const nrtfn = mem.cv_nrtfn;
  const N = mem.cv_N;

  if (mem.cv_irfnd === 0) {
    return CV_SUCCESS;
  }

  // Get y at tlo via interpolation
  cvodeGetDky(mem, mem.cv_tlo, 0, mem.cv_y);

  // Evaluate g(tlo, y)
  let retval = mem.cv_gfun!(mem.cv_tlo, mem.cv_y, mem.cv_glo, mem.cv_user_data);
  mem.cv_nge++;
  if (retval !== 0) {
    return CV_RTFUNC_FAIL;
  }

  let zroot = false;
  for (let i = 0; i < nrtfn; i++) {
    mem.cv_iroots[i] = 0;
  }
  for (let i = 0; i < nrtfn; i++) {
    if (mem.cv_gactive[i] === 0) continue;
    if (Math.abs(mem.cv_glo[i]) === ZERO) {
      zroot = true;
      mem.cv_iroots[i] = 1;
    }
  }
  if (!zroot) {
    return CV_SUCCESS;
  }

  // One or more g_i has a zero at tlo. Check g at tlo+smallh.
  mem.cv_ttol = (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h)) * UROUND * HUNDRED;
  const smallh = (mem.cv_h > ZERO) ? mem.cv_ttol : -mem.cv_ttol;
  const tplus = mem.cv_tlo + smallh;

  if ((tplus - mem.cv_tn) * mem.cv_h >= ZERO) {
    const hratio = smallh / mem.cv_h;
    // y = y + hratio * zn[1]
    vLinearSum(ONE, mem.cv_y, hratio, mem.cv_zn[1], mem.cv_y, N);
  } else {
    cvodeGetDky(mem, tplus, 0, mem.cv_y);
  }

  retval = mem.cv_gfun!(tplus, mem.cv_y, mem.cv_ghi, mem.cv_user_data);
  mem.cv_nge++;
  if (retval !== 0) {
    return CV_RTFUNC_FAIL;
  }

  // Check for close roots (error return), new zero at tlo+smallh,
  // and g_i that changed from zero to nonzero.
  zroot = false;
  for (let i = 0; i < nrtfn; i++) {
    if (mem.cv_gactive[i] === 0) continue;
    if (Math.abs(mem.cv_ghi[i]) === ZERO) {
      if (mem.cv_iroots[i] === 1) {
        return CLOSERT;
      }
      zroot = true;
      mem.cv_iroots[i] = 1;
    } else {
      if (mem.cv_iroots[i] === 1) {
        mem.cv_glo[i] = mem.cv_ghi[i];
      }
    }
  }
  if (zroot) {
    return RTFOUND;
  }

  return CV_SUCCESS;
}

/**
 * cvRcheck3 - Root search between tlo and tn (or tout).
 *
 * C reference: cvRcheck3 (cvode.c, lines ~4400-4453)
 *
 * Interfaces to cvRootfind to look for a root of g between tlo and
 * either tn or tout, whichever comes first.
 *
 * @param mem - CVODE memory block
 * @returns CV_SUCCESS, RTFOUND, or CV_RTFUNC_FAIL
 */
export function cvRcheck3(mem: CvodeMem): number {
  const nrtfn = mem.cv_nrtfn;
  const N = mem.cv_N;

  // Set thi = tn or tout, whichever comes first; set y = y(thi)
  if (mem.cv_taskc === CV_ONE_STEP) {
    mem.cv_thi = mem.cv_tn;
    vScale(ONE, mem.cv_zn[0], mem.cv_y, N);
  }
  if (mem.cv_taskc === CV_NORMAL) {
    if ((mem.cv_toutc - mem.cv_tn) * mem.cv_h >= ZERO) {
      mem.cv_thi = mem.cv_tn;
      vScale(ONE, mem.cv_zn[0], mem.cv_y, N);
    } else {
      mem.cv_thi = mem.cv_toutc;
      cvodeGetDky(mem, mem.cv_thi, 0, mem.cv_y);
    }
  }

  // Evaluate g(thi, y) and call cvRootfind
  const retval = mem.cv_gfun!(mem.cv_thi, mem.cv_y, mem.cv_ghi, mem.cv_user_data);
  mem.cv_nge++;
  if (retval !== 0) {
    return CV_RTFUNC_FAIL;
  }

  mem.cv_ttol = (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h)) * UROUND * HUNDRED;
  const ier = cvRootfind(mem);
  if (ier === CV_RTFUNC_FAIL) {
    return CV_RTFUNC_FAIL;
  }

  // Re-activate g components that moved away from zero
  for (let i = 0; i < nrtfn; i++) {
    if (mem.cv_gactive[i] === 0 && mem.cv_grout[i] !== ZERO) {
      mem.cv_gactive[i] = 1; // SUNTRUE
    }
  }

  mem.cv_tlo = mem.cv_trout;
  for (let i = 0; i < nrtfn; i++) {
    mem.cv_glo[i] = mem.cv_grout[i];
  }

  // If no root found, return CV_SUCCESS
  if (ier === CV_SUCCESS) {
    return CV_SUCCESS;
  }

  // If a root was found, interpolate to get y(trout) and return
  cvodeGetDky(mem, mem.cv_trout, 0, mem.cv_y);
  return RTFOUND;
}

// ============================================================================
// Internal (non-exported) functions
// ============================================================================

/**
 * cvRootfind - Illinois algorithm for root finding.
 *
 * C reference: cvRootfind (cvode.c, lines ~4532-4737)
 *
 * Solves for a root of g(t) between tlo and thi using a modified secant
 * method (Illinois algorithm). Only roots of odd multiplicity (sign changes)
 * or exact zeros are found. When multiple roots exist, the one closest to
 * tlo is returned.
 *
 * @param mem - CVODE memory block
 * @returns CV_SUCCESS, RTFOUND, or CV_RTFUNC_FAIL
 */
function cvRootfind(mem: CvodeMem): number {
  const nrtfn = mem.cv_nrtfn;

  let imax = 0;
  let alph: number;
  let tmid: number;
  let gfrac: number;
  let maxfrac: number;
  let fracint: number;
  let fracsub: number;

  // First check for change in sign in ghi or for a zero in ghi.
  maxfrac = ZERO;
  let zroot = false;
  let sgnchg = false;

  for (let i = 0; i < nrtfn; i++) {
    if (mem.cv_gactive[i] === 0) continue;
    if (Math.abs(mem.cv_ghi[i]) === ZERO) {
      if (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO) {
        zroot = true;
      }
    } else {
      if ((mem.cv_glo[i] * mem.cv_ghi[i] < 0) &&
          (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO)) {
        gfrac = Math.abs(mem.cv_ghi[i] / (mem.cv_ghi[i] - mem.cv_glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = true;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  // If no sign change was found, reset trout and grout.
  // Return CV_SUCCESS if no zero was found, or set iroots and return RTFOUND.
  if (!sgnchg) {
    mem.cv_trout = mem.cv_thi;
    for (let i = 0; i < nrtfn; i++) {
      mem.cv_grout[i] = mem.cv_ghi[i];
    }
    if (!zroot) {
      return CV_SUCCESS;
    }
    for (let i = 0; i < nrtfn; i++) {
      mem.cv_iroots[i] = 0;
      if (mem.cv_gactive[i] === 0) continue;
      if ((Math.abs(mem.cv_ghi[i]) === ZERO) &&
          (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO)) {
        mem.cv_iroots[i] = mem.cv_glo[i] > 0 ? -1 : 1;
      }
    }
    return RTFOUND;
  }

  // Initialize alph
  alph = ONE;

  // A sign change was found. Loop to locate nearest root.
  let side = 0;
  let sideprev = -1;

  for (;;) {
    // If interval size is already less than tolerance ttol, break.
    if (Math.abs(mem.cv_thi - mem.cv_tlo) <= mem.cv_ttol) {
      break;
    }

    // Set weight alph.
    // On the first two passes, set alph = 1. Thereafter, reset alph
    // according to the side of the subinterval where the sign change was found.
    if (sideprev === side) {
      alph = (side === 2) ? alph * TWO : alph * HALF;
    } else {
      alph = ONE;
    }

    // Set next root approximation tmid and get g(tmid).
    // If tmid is too close to tlo or thi, adjust it inward.
    tmid = mem.cv_thi -
      (mem.cv_thi - mem.cv_tlo) * mem.cv_ghi[imax] /
      (mem.cv_ghi[imax] - alph * mem.cv_glo[imax]);

    if (Math.abs(tmid - mem.cv_tlo) < HALF * mem.cv_ttol) {
      fracint = Math.abs(mem.cv_thi - mem.cv_tlo) / mem.cv_ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF / fracint;
      tmid = mem.cv_tlo + fracsub * (mem.cv_thi - mem.cv_tlo);
    }
    if (Math.abs(mem.cv_thi - tmid) < HALF * mem.cv_ttol) {
      fracint = Math.abs(mem.cv_thi - mem.cv_tlo) / mem.cv_ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF / fracint;
      tmid = mem.cv_thi - fracsub * (mem.cv_thi - mem.cv_tlo);
    }

    // Interpolate to get y(tmid) and evaluate g(tmid)
    cvodeGetDky(mem, tmid, 0, mem.cv_y);
    const retval = mem.cv_gfun!(tmid, mem.cv_y, mem.cv_grout, mem.cv_user_data);
    mem.cv_nge++;
    if (retval !== 0) {
      return CV_RTFUNC_FAIL;
    }

    // Check in which subinterval g changes sign, and reset imax.
    // Set side = 1 if sign change is on low side, or 2 if on high side.
    maxfrac = ZERO;
    zroot = false;
    sgnchg = false;
    sideprev = side;

    for (let i = 0; i < nrtfn; i++) {
      if (mem.cv_gactive[i] === 0) continue;
      if (Math.abs(mem.cv_grout[i]) === ZERO) {
        if (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO) {
          zroot = true;
        }
      } else {
        if ((mem.cv_glo[i] * mem.cv_grout[i] < 0) &&
            (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO)) {
          gfrac = Math.abs(mem.cv_grout[i] /
            (mem.cv_grout[i] - mem.cv_glo[i]));
          if (gfrac > maxfrac) {
            sgnchg = true;
            maxfrac = gfrac;
            imax = i;
          }
        }
      }
    }

    if (sgnchg) {
      // Sign change found in (tlo, tmid); replace thi with tmid.
      mem.cv_thi = tmid;
      for (let i = 0; i < nrtfn; i++) {
        mem.cv_ghi[i] = mem.cv_grout[i];
      }
      side = 1;
      // Stop at root thi if converged; otherwise loop.
      if (Math.abs(mem.cv_thi - mem.cv_tlo) <= mem.cv_ttol) {
        break;
      }
      continue;
    }

    if (zroot) {
      // No sign change in (tlo, tmid), but g = 0 at tmid; return root tmid.
      mem.cv_thi = tmid;
      for (let i = 0; i < nrtfn; i++) {
        mem.cv_ghi[i] = mem.cv_grout[i];
      }
      break;
    }

    // No sign change in (tlo, tmid), and no zero at tmid.
    // Sign change must be in (tmid, thi). Replace tlo with tmid.
    mem.cv_tlo = tmid;
    for (let i = 0; i < nrtfn; i++) {
      mem.cv_glo[i] = mem.cv_grout[i];
    }
    side = 2;
    // Stop at root thi if converged; otherwise loop back.
    if (Math.abs(mem.cv_thi - mem.cv_tlo) <= mem.cv_ttol) {
      break;
    }
  } // End of root-search loop

  // Reset trout and grout, set iroots, and return RTFOUND.
  mem.cv_trout = mem.cv_thi;
  for (let i = 0; i < nrtfn; i++) {
    mem.cv_grout[i] = mem.cv_ghi[i];
    mem.cv_iroots[i] = 0;
    if (mem.cv_gactive[i] === 0) continue;
    if ((Math.abs(mem.cv_ghi[i]) === ZERO) &&
        (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO)) {
      mem.cv_iroots[i] = mem.cv_glo[i] > 0 ? -1 : 1;
    }
    if ((mem.cv_glo[i] * mem.cv_ghi[i] < 0) &&
        (mem.cv_rootdir[i] * mem.cv_glo[i] <= ZERO)) {
      mem.cv_iroots[i] = mem.cv_glo[i] > 0 ? -1 : 1;
    }
  }

  return RTFOUND;
}
