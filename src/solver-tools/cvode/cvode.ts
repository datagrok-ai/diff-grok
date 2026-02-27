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

// CVODE TypeScript port - Main integration loop
// Ported from SUNDIALS CVODE v7.5.0 cvode.c

import {
  CvodeMem, CvodeRhsFn,
  CV_ADAMS, CV_BDF,
  CV_NORMAL, CV_ONE_STEP,
  CV_SUCCESS, CV_TSTOP_RETURN, CV_ROOT_RETURN,
  CV_TOO_MUCH_WORK, CV_TOO_MUCH_ACC, CV_ERR_FAILURE, CV_CONV_FAILURE,
  CV_RHSFUNC_FAIL, CV_FIRST_RHSFUNC_ERR, CV_RTFUNC_FAIL, CV_TOO_CLOSE,
  CV_LSETUP_FAIL, CV_LSOLVE_FAIL, CV_UNREC_RHSFUNC_ERR, CV_REPTD_RHSFUNC_ERR,
  CV_LINIT_FAIL, CV_NLS_FAIL,
  CV_NN, CV_SS, CV_SV,
  ADAMS_Q_MAX, BDF_Q_MAX,
  HMIN_DEFAULT, HMAX_INV_DEFAULT, MXSTEP_DEFAULT, MXHNIL_DEFAULT,
  MSBP_DEFAULT, DGMAX_LSETUP_DEFAULT,
  ETA_MIN_FX_DEFAULT, ETA_MAX_FX_DEFAULT, ETA_MAX_FS_DEFAULT,
  ETA_MAX_ES_DEFAULT, ETA_MAX_GS_DEFAULT,
  ETA_MIN_DEFAULT, ETA_MIN_EF_DEFAULT, ETA_MAX_EF_DEFAULT, ETA_CF_DEFAULT,
  SMALL_NST_DEFAULT, SMALL_NEF_DEFAULT,
  MXNEF, MXNCF, MXNEF1,
  DO_ERROR_TEST, PREDICT_AGAIN, TRY_AGAIN, FIRST_CALL,
  PREV_CONV_FAIL, PREV_ERR_FAIL,
  RHSFUNC_RECVR, SUN_NLS_CONV_RECVR,
  CV_NO_FAILURES,
  FUZZ_FACTOR, UROUND, CORTES, ONEPSM, TINY,
  ADDON, BIAS1, BIAS2, BIAS3, LONG_WAIT,
  NUM_TESTS, L_MAX,
  RTFOUND, CLOSERT,
  wrmsNorm, vLinearSum, vScale, cvodeGetDky,
} from './common';

import { cvSetBDF, cvIncreaseBDF, cvDecreaseBDF } from './cvode_bdf';
import { cvSetAdams, cvAdjustAdams } from './cvode_adams';
import { cvHin } from './cvode_hin';
import { cvNls } from './cvode_nls';
import { cvRcheck1, cvRcheck2, cvRcheck3 } from './cvode_root';

// ============================================================================
// PUBLIC EXPORTED FUNCTIONS
// ============================================================================

/** cvodeCreate - Create and initialize CvodeMem */
export function cvodeCreate(lmm: number): CvodeMem {
  if (lmm !== CV_ADAMS && lmm !== CV_BDF) {
    throw new Error('cvodeCreate: bad lmm, must be CV_ADAMS or CV_BDF');
  }

  const maxord = (lmm === CV_ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;
  const mem = new CvodeMem();

  mem.cv_lmm = lmm;
  mem.cv_qmax = maxord;
  mem.cv_qmax_alloc = maxord;

  // Initialize ssdat as 6x4 array of zeros (1-based indexing: i=1..5, k=1..3)
  mem.cv_ssdat = [];
  for (let i = 0; i < 6; i++) {
    mem.cv_ssdat[i] = [0, 0, 0, 0];
  }

  mem.cv_nlscoef = CORTES;

  return mem;
}

/** cvodeInit - Initialize problem */
export function cvodeInit(mem: CvodeMem, f: CvodeRhsFn, t0: number, y0: Float64Array): number {
  mem.cv_f = f;
  mem.cv_tn = t0;
  mem.cv_N = y0.length;
  const N = mem.cv_N;

  // Allocate Nordsieck history array: zn[0..qmax_alloc+1]
  const allocQ = mem.cv_qmax_alloc;
  mem.cv_zn = [];
  for (let i = 0; i <= allocQ + 1; i++) {
    mem.cv_zn[i] = new Float64Array(N);
  }

  // Allocate work vectors
  mem.cv_ewt = new Float64Array(N);
  mem.cv_acor = new Float64Array(N);
  mem.cv_tempv = new Float64Array(N);
  mem.cv_ftemp = new Float64Array(N);
  mem.cv_vtemp1 = new Float64Array(N);
  mem.cv_vtemp2 = new Float64Array(N);
  mem.cv_vtemp3 = new Float64Array(N);
  mem.cv_y = new Float64Array(N);

  // Copy y0 into zn[0]
  mem.cv_zn[0].set(y0);

  // Set step parameters
  mem.cv_q = 1;
  mem.cv_L = 2;
  mem.cv_qwait = mem.cv_L;
  mem.cv_etamax = mem.cv_eta_max_fs;

  mem.cv_qu = 0;
  mem.cv_hu = 0.0;
  mem.cv_tolsf = 1.0;

  // Reset all counters
  mem.cv_nst = 0;
  mem.cv_nfe = 0;
  mem.cv_ncfn = 0;
  mem.cv_netf = 0;
  mem.cv_nni = 0;
  mem.cv_nnf = 0;
  mem.cv_nsetups = 0;
  mem.cv_nhnil = 0;
  mem.cv_nstlp = 0;
  mem.cv_nscon = 0;
  mem.cv_nge = 0;
  mem.cv_irfnd = 0;

  mem.cv_h0u = 0.0;
  mem.cv_next_h = 0.0;
  mem.cv_next_q = 0;

  // Initialize stability limit detection data
  mem.cv_nor = 0;
  mem.cv_ssdat = [];
  for (let i = 0; i < 6; i++) {
    mem.cv_ssdat[i] = [0, 0, 0, 0];
  }

  // Set indx_acor to qmax_alloc (the zn slot used to save acor)
  mem.cv_indx_acor = mem.cv_qmax_alloc;

  mem.cv_MallocDone = true;

  return CV_SUCCESS;
}

/** cvode - Main integration loop */
export function cvode(
  mem: CvodeMem, tout: number, yout: Float64Array, itask: number,
): { flag: number; t: number } {
  let istate: number;
  let troundoff: number;
  let retval: number;
  let ier: number;
  let tret = mem.cv_tn;

  // Check itask
  if (itask !== CV_NORMAL && itask !== CV_ONE_STEP) {
    return { flag: -99, t: mem.cv_tn };
  }

  if (itask === CV_NORMAL) { mem.cv_toutc = tout; }
  mem.cv_taskc = itask;

  // --- nst == 0: first step initialization ---
  if (mem.cv_nst === 0) {
    mem.cv_tretlast = tret = mem.cv_tn;

    // Check inputs
    ier = cvInitialSetup(mem);
    if (ier !== CV_SUCCESS) { return { flag: ier, t: tret }; }

    // Call f at (t0, y0), set zn[1] = y'(t0)
    retval = mem.cv_f!(mem.cv_tn, mem.cv_zn[0], mem.cv_zn[1], mem.cv_user_data);
    mem.cv_nfe++;
    if (retval < 0) { return { flag: CV_RHSFUNC_FAIL, t: tret }; }
    if (retval > 0) { return { flag: CV_FIRST_RHSFUNC_ERR, t: tret }; }

    // Test tstop legality
    if (mem.cv_tstopset) {
      if ((mem.cv_tstop - mem.cv_tn) * (tout - mem.cv_tn) <= 0.0) {
        return { flag: -99, t: tret }; // CV_ILL_INPUT
      }
    }

    // Set initial h
    mem.cv_h = mem.cv_hin;
    if (mem.cv_h !== 0.0 && (tout - mem.cv_tn) * mem.cv_h < 0.0) {
      return { flag: -99, t: tret };
    }
    if (mem.cv_h === 0.0) {
      let tout_hin = tout;
      if (mem.cv_tstopset && (tout - mem.cv_tn) * (tout - mem.cv_tstop) > 0.0) {
        tout_hin = mem.cv_tstop;
      }
      const hflag = cvHin(mem, tout_hin);
      if (hflag !== CV_SUCCESS) {
        istate = cvHandleFailure(mem, hflag);
        return { flag: istate, t: tret };
      }
    }

    // Enforce hmax and hmin
    let rh = Math.abs(mem.cv_h) * mem.cv_hmax_inv;
    if (rh > 1.0) { mem.cv_h /= rh; }
    if (Math.abs(mem.cv_h) < mem.cv_hmin) {
      mem.cv_h *= mem.cv_hmin / Math.abs(mem.cv_h);
    }

    // Check for approach to tstop
    if (mem.cv_tstopset) {
      if ((mem.cv_tn + mem.cv_h - mem.cv_tstop) * mem.cv_h > 0.0) {
        mem.cv_h = (mem.cv_tstop - mem.cv_tn) * (1.0 - 4.0 * UROUND);
      }
    }

    // Scale zn[1] by h
    mem.cv_hscale = mem.cv_h;
    mem.cv_h0u = mem.cv_h;
    mem.cv_hprime = mem.cv_h;
    vScale(mem.cv_h, mem.cv_zn[1], mem.cv_zn[1], mem.cv_N);

    // Check for roots at/near t0
    if (mem.cv_nrtfn > 0) {
      retval = cvRcheck1(mem);
      if (retval === CV_RTFUNC_FAIL) {
        return { flag: CV_RTFUNC_FAIL, t: tret };
      }
    }
  }

  // --- nst > 0: stop tests ---
  if (mem.cv_nst > 0) {
    troundoff = FUZZ_FACTOR * UROUND *
      (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h));

    // Root check from last step
    if (mem.cv_nrtfn > 0) {
      const irfndp = mem.cv_irfnd;

      retval = cvRcheck2(mem);

      if (retval === CLOSERT) {
        return { flag: -99, t: mem.cv_tlo };
      } else if (retval === CV_RTFUNC_FAIL) {
        return { flag: CV_RTFUNC_FAIL, t: mem.cv_tlo };
      } else if (retval === RTFOUND) {
        mem.cv_tretlast = tret = mem.cv_tlo;
        return { flag: CV_ROOT_RETURN, t: tret };
      }

      // If tn is distinct from tretlast, check remaining interval
      if (Math.abs(mem.cv_tn - mem.cv_tretlast) > troundoff) {
        retval = cvRcheck3(mem);

        if (retval === CV_SUCCESS) {
          mem.cv_irfnd = 0;
          if (irfndp === 1 && itask === CV_ONE_STEP) {
            mem.cv_tretlast = tret = mem.cv_tn;
            yout.set(mem.cv_zn[0]);
            return { flag: CV_SUCCESS, t: tret };
          }
        } else if (retval === RTFOUND) {
          mem.cv_irfnd = 1;
          mem.cv_tretlast = tret = mem.cv_tlo;
          return { flag: CV_ROOT_RETURN, t: tret };
        } else if (retval === CV_RTFUNC_FAIL) {
          return { flag: CV_RTFUNC_FAIL, t: mem.cv_tlo };
        }
      }
    }

    // Test for tn at tstop
    if (mem.cv_tstopset) {
      if (Math.abs(mem.cv_tn - mem.cv_tstop) <= troundoff) {
        if ((tout - mem.cv_tstop) * mem.cv_h >= 0.0 ||
            Math.abs(tout - mem.cv_tstop) <= troundoff) {
          if (mem.cv_tstopinterp) {
            cvodeGetDky(mem, mem.cv_tstop, 0, yout);
          } else {
            yout.set(mem.cv_zn[0]);
          }
          mem.cv_tretlast = tret = mem.cv_tstop;
          mem.cv_tstopset = false;
          return { flag: CV_TSTOP_RETURN, t: tret };
        }
      } else if ((mem.cv_tn + mem.cv_hprime - mem.cv_tstop) * mem.cv_h > 0.0) {
        mem.cv_hprime = (mem.cv_tstop - mem.cv_tn) * (1.0 - 4.0 * UROUND);
        mem.cv_eta = mem.cv_hprime / mem.cv_h;
      }
    }

    // In NORMAL mode, test if tout was reached
    if (itask === CV_NORMAL && (mem.cv_tn - tout) * mem.cv_h >= 0.0) {
      mem.cv_tretlast = tret = tout;
      cvodeGetDky(mem, tout, 0, yout);
      return { flag: CV_SUCCESS, t: tret };
    }

    // In ONE_STEP mode, test if tn was returned
    if (itask === CV_ONE_STEP &&
        Math.abs(mem.cv_tn - mem.cv_tretlast) > troundoff) {
      mem.cv_tretlast = tret = mem.cv_tn;
      yout.set(mem.cv_zn[0]);
      return { flag: CV_SUCCESS, t: tret };
    }
  }

  // --- Main integration loop ---
  let nstloc = 0;
  istate = CV_SUCCESS; // will be overwritten

  for (;;) {
    mem.cv_next_h = mem.cv_h;
    mem.cv_next_q = mem.cv_q;

    // Reset and check ewt
    if (mem.cv_nst > 0) {
      const ewtsetOK = cvEwtSet(mem);
      if (ewtsetOK !== 0) {
        istate = -99; // CV_ILL_INPUT
        mem.cv_tretlast = tret = mem.cv_tn;
        yout.set(mem.cv_zn[0]);
        break;
      }
    }

    // Check for too many steps
    if (mem.cv_mxstep > 0 && nstloc >= mem.cv_mxstep) {
      istate = CV_TOO_MUCH_WORK;
      mem.cv_tretlast = tret = mem.cv_tn;
      yout.set(mem.cv_zn[0]);
      break;
    }

    // Check for too much accuracy requested
    const nrm = wrmsNorm(mem.cv_N, mem.cv_zn[0], mem.cv_ewt);
    mem.cv_tolsf = UROUND * nrm;
    if (mem.cv_tolsf > 1.0) {
      istate = CV_TOO_MUCH_ACC;
      mem.cv_tretlast = tret = mem.cv_tn;
      yout.set(mem.cv_zn[0]);
      mem.cv_tolsf *= 2.0;
      break;
    } else {
      mem.cv_tolsf = 1.0;
    }

    // Check for h below roundoff
    if (mem.cv_tn + mem.cv_h === mem.cv_tn) {
      mem.cv_nhnil++;
    }

    // Take a step
    const kflag = cvStep(mem);

    // Process failed step
    if (kflag !== CV_SUCCESS) {
      istate = cvHandleFailure(mem, kflag);
      mem.cv_tretlast = tret = mem.cv_tn;
      yout.set(mem.cv_zn[0]);
      break;
    }

    nstloc++;

    // If tstop was reached, reset tn
    if (mem.cv_tstopset) {
      troundoff = FUZZ_FACTOR * UROUND *
        (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h));
      if (Math.abs(mem.cv_tn - mem.cv_tstop) <= troundoff) {
        mem.cv_tn = mem.cv_tstop;
      }
    }

    // Check for root in last step
    if (mem.cv_nrtfn > 0) {
      retval = cvRcheck3(mem);

      if (retval === RTFOUND) {
        mem.cv_irfnd = 1;
        istate = CV_ROOT_RETURN;
        mem.cv_tretlast = tret = mem.cv_tlo;
        break;
      } else if (retval === CV_RTFUNC_FAIL) {
        istate = CV_RTFUNC_FAIL;
        break;
      }
    }

    // Check tstop in loop
    if (mem.cv_tstopset) {
      troundoff = FUZZ_FACTOR * UROUND *
        (Math.abs(mem.cv_tn) + Math.abs(mem.cv_h));

      if (Math.abs(mem.cv_tn - mem.cv_tstop) <= troundoff) {
        if ((tout - mem.cv_tstop) * mem.cv_h >= 0.0 ||
            Math.abs(tout - mem.cv_tstop) <= troundoff) {
          if (mem.cv_tstopinterp) {
            cvodeGetDky(mem, mem.cv_tstop, 0, yout);
          } else {
            yout.set(mem.cv_zn[0]);
          }
          mem.cv_tretlast = tret = mem.cv_tstop;
          mem.cv_tstopset = false;
          istate = CV_TSTOP_RETURN;
          break;
        }
      } else if ((mem.cv_tn + mem.cv_hprime - mem.cv_tstop) * mem.cv_h > 0.0) {
        mem.cv_hprime = (mem.cv_tstop - mem.cv_tn) * (1.0 - 4.0 * UROUND);
        mem.cv_eta = mem.cv_hprime / mem.cv_h;
      }
    }

    // In NORMAL mode, check if tout reached
    if (itask === CV_NORMAL && (mem.cv_tn - tout) * mem.cv_h >= 0.0) {
      istate = CV_SUCCESS;
      mem.cv_tretlast = tret = tout;
      cvodeGetDky(mem, tout, 0, yout);
      mem.cv_next_q = mem.cv_qprime;
      mem.cv_next_h = mem.cv_hprime;
      break;
    }

    // In ONE_STEP mode, copy y and exit
    if (itask === CV_ONE_STEP) {
      istate = CV_SUCCESS;
      mem.cv_tretlast = tret = mem.cv_tn;
      yout.set(mem.cv_zn[0]);
      mem.cv_next_q = mem.cv_qprime;
      mem.cv_next_h = mem.cv_hprime;
      break;
    }
  }

  return { flag: istate, t: tret };
}

// ============================================================================
// INTERNAL FUNCTIONS
// ============================================================================

/** cvInitialSetup - First-step checks */
function cvInitialSetup(mem: CvodeMem): number {
  // Did the user specify tolerances?
  if (mem.cv_itol === CV_NN) {
    mem.cv_error = 'cvInitialSetup: no tolerances set';
    return -99; // CV_ILL_INPUT
  }

  // Load initial error weights
  const ier = cvEwtSet(mem);
  if (ier !== 0) {
    mem.cv_error = 'cvInitialSetup: bad ewt';
    return -99; // CV_ILL_INPUT
  }

  // Call linit function if it exists
  if (mem.cv_linit !== null) {
    const linit_ret = mem.cv_linit(mem);
    if (linit_ret !== 0) {
      mem.cv_error = 'cvInitialSetup: linear solver init failed';
      return CV_LINIT_FAIL;
    }
  }

  return CV_SUCCESS;
}

/** cvEwtSet - Set error weights based on tolerance type */
function cvEwtSet(mem: CvodeMem): number {
  switch (mem.cv_itol) {
    case CV_SS: return cvEwtSetSS(mem);
    case CV_SV: return cvEwtSetSV(mem);
    default: return -1;
  }
}

/** cvEwtSetSS - ewt[i] = 1/(reltol*|y[i]| + Sabstol) */
function cvEwtSetSS(mem: CvodeMem): number {
  const N = mem.cv_N;
  const rtol = mem.cv_reltol;
  const atol = mem.cv_Sabstol;
  const y = mem.cv_zn[0];
  const ewt = mem.cv_ewt;

  for (let i = 0; i < N; i++) {
    const tol = rtol * Math.abs(y[i]) + atol;
    if (mem.cv_atolmin0 && tol <= 0.0) { return -1; }
    ewt[i] = 1.0 / tol;
  }
  return 0;
}

/** cvEwtSetSV - ewt[i] = 1/(reltol*|y[i]| + Vabstol[i]) */
function cvEwtSetSV(mem: CvodeMem): number {
  const N = mem.cv_N;
  const rtol = mem.cv_reltol;
  const vatol = mem.cv_Vabstol;
  const y = mem.cv_zn[0];
  const ewt = mem.cv_ewt;

  for (let i = 0; i < N; i++) {
    const tol = rtol * Math.abs(y[i]) + vatol[i];
    if (mem.cv_atolmin0 && tol <= 0.0) { return -1; }
    ewt[i] = 1.0 / tol;
  }
  return 0;
}

/** cvStep - One internal step */
function cvStep(mem: CvodeMem): number {
  let dsm = 0.0;
  let ncf = 0;
  let nef = 0;
  let nflag: number;
  let kflag: number;
  let eflag: number;

  // If step size changed, update history array
  if (mem.cv_nst > 0 && mem.cv_hprime !== mem.cv_h) {
    cvAdjustParams(mem);
  }

  const saved_t = mem.cv_tn;
  nflag = FIRST_CALL;

  // Looping point for step attempts
  for (;;) {
    cvPredict(mem);
    cvSet(mem);

    nflag = cvNls(mem, nflag);

    const nflagRef = { value: nflag };
    const ncfRef = { value: ncf };
    kflag = cvHandleNFlag(mem, nflagRef, saved_t, ncfRef);
    nflag = nflagRef.value;
    ncf = ncfRef.value;

    // Go back if we need to predict again
    if (kflag === PREDICT_AGAIN) { continue; }

    // Return if nonlinear solve failed and recovery is not possible
    if (kflag !== DO_ERROR_TEST) { return kflag; }

    // Perform error test
    const nflagRef2 = { value: nflag };
    const nefRef = { value: nef };
    const dsmRef = { value: dsm };
    eflag = cvDoErrorTest(mem, nflagRef2, saved_t, nefRef, dsmRef);
    nflag = nflagRef2.value;
    nef = nefRef.value;
    dsm = dsmRef.value;

    // Go back if we need to try again
    if (eflag === TRY_AGAIN) { continue; }

    // Return if error test failed and recovery is not possible
    if (eflag !== CV_SUCCESS) { return eflag; }

    // Error test passed, break from loop
    break;
  }

  // Step successful: update data and consider order/step changes
  cvCompleteStep(mem);
  cvPrepareNextStep(mem, dsm);

  // BDF stability limit detection
  if (mem.cv_sldeton) { cvBDFStab(mem); }

  mem.cv_etamax = (mem.cv_nst <= mem.cv_small_nst)
    ? mem.cv_eta_max_es
    : mem.cv_eta_max_gs;

  // Rescale acor to be the estimated local error vector
  vScale(mem.cv_tq[2], mem.cv_acor, mem.cv_acor, mem.cv_N);

  return CV_SUCCESS;
}

/** cvAdjustParams - adjust order and rescale when step size changes */
function cvAdjustParams(mem: CvodeMem): void {
  if (mem.cv_qprime !== mem.cv_q) {
    if (!mem.first_step_after_resize) {
      cvAdjustOrder(mem, mem.cv_qprime - mem.cv_q);
    }
    mem.cv_q = mem.cv_qprime;
    mem.cv_L = mem.cv_q + 1;
    mem.cv_qwait = mem.cv_L;
  }
  cvRescale(mem);
}

/** cvAdjustOrder - handle order change by deltaq */
function cvAdjustOrder(mem: CvodeMem, deltaq: number): void {
  if (mem.cv_q === 2 && deltaq !== 1) { return; }

  if (mem.cv_lmm === CV_ADAMS) {
    cvAdjustAdams(mem, deltaq);
  } else if (mem.cv_lmm === CV_BDF) {
    if (deltaq === 1) { cvIncreaseBDF(mem); }
    else if (deltaq === -1) { cvDecreaseBDF(mem); }
  }
}

/** cvRescale - multiply zn[j] by eta^j for j=1..q */
function cvRescale(mem: CvodeMem): void {
  const N = mem.cv_N;
  let factor = mem.cv_eta;
  for (let j = 1; j <= mem.cv_q; j++) {
    const znj = mem.cv_zn[j];
    for (let i = 0; i < N; i++) {
      znj[i] *= factor;
    }
    factor *= mem.cv_eta;
  }
  mem.cv_h = mem.cv_hscale * mem.cv_eta;
  mem.cv_next_h = mem.cv_h;
  mem.cv_hscale = mem.cv_h;
  mem.cv_nscon = 0;
}

/** cvPredict - advance tn and compute predicted Nordsieck array */
function cvPredict(mem: CvodeMem): void {
  const N = mem.cv_N;

  mem.cv_tn += mem.cv_h;
  if (mem.cv_tstopset) {
    if ((mem.cv_tn - mem.cv_tstop) * mem.cv_h > 0.0) {
      mem.cv_tn = mem.cv_tstop;
    }
  }

  for (let k = 1; k <= mem.cv_q; k++) {
    for (let j = mem.cv_q; j >= k; j--) {
      const znjm1 = mem.cv_zn[j - 1];
      const znj = mem.cv_zn[j];
      for (let i = 0; i < N; i++) {
        znjm1[i] += znj[i];
      }
    }
  }
}

/** cvSet - compute method coefficients l, tq, rl1, gamma, gamrat */
function cvSet(mem: CvodeMem): void {
  if (mem.cv_lmm === CV_ADAMS) { cvSetAdams(mem); }
  else { cvSetBDF(mem); }

  mem.cv_rl1 = 1.0 / mem.cv_l[1];
  mem.cv_gamma = mem.cv_h * mem.cv_rl1;
  if (mem.cv_nst === 0) { mem.cv_gammap = mem.cv_gamma; }
  mem.cv_gamrat = (mem.cv_nst > 0)
    ? mem.cv_gamma / mem.cv_gammap
    : 1.0;
}

/** cvHandleNFlag - take action on nonlinear solver return */
function cvHandleNFlag(
  mem: CvodeMem, nflagRef: { value: number }, saved_t: number, ncfRef: { value: number },
): number {
  const nflag = nflagRef.value;

  if (nflag === CV_SUCCESS) { return DO_ERROR_TEST; }

  // Nonlinear solve failed; increment ncfn and restore zn
  mem.cv_ncfn++;
  cvRestore(mem, saved_t);

  // Return if failed unrecoverably
  if (nflag < 0) {
    if (nflag === CV_LSETUP_FAIL) { return CV_LSETUP_FAIL; }
    else if (nflag === CV_LSOLVE_FAIL) { return CV_LSOLVE_FAIL; }
    else if (nflag === CV_RHSFUNC_FAIL) { return CV_RHSFUNC_FAIL; }
    else { return CV_NLS_FAIL; }
  }

  // Recoverable error
  ncfRef.value++;
  mem.cv_etamax = 1.0;

  // If maxncf failures or |h| = hmin, return failure
  if (Math.abs(mem.cv_h) <= mem.cv_hmin * ONEPSM ||
      ncfRef.value === mem.cv_maxncf) {
    if (nflag === SUN_NLS_CONV_RECVR) { return CV_CONV_FAILURE; }
    if (nflag === RHSFUNC_RECVR) { return CV_REPTD_RHSFUNC_ERR; }
  }

  // Reduce step size; return to reattempt
  mem.cv_eta = Math.max(mem.cv_eta_cf, mem.cv_hmin / Math.abs(mem.cv_h));
  nflagRef.value = PREV_CONV_FAIL;
  cvRescale(mem);

  return PREDICT_AGAIN;
}

/** cvRestore - undo prediction, restore tn and zn */
function cvRestore(mem: CvodeMem, saved_t: number): void {
  const N = mem.cv_N;

  mem.cv_tn = saved_t;
  for (let k = 1; k <= mem.cv_q; k++) {
    for (let j = mem.cv_q; j >= k; j--) {
      const znjm1 = mem.cv_zn[j - 1];
      const znj = mem.cv_zn[j];
      for (let i = 0; i < N; i++) {
        znjm1[i] -= znj[i];
      }
    }
  }
}

/** cvDoErrorTest - local error test */
function cvDoErrorTest(
  mem: CvodeMem, nflagRef: { value: number }, saved_t: number,
  nefRef: { value: number }, dsmRef: { value: number },
): number {
  const dsm = mem.cv_acnrm * mem.cv_tq[2];
  dsmRef.value = dsm;

  // If dsm <= 1, test passed
  if (dsm <= 1.0) { return CV_SUCCESS; }

  // Test failed
  nefRef.value++;
  mem.cv_netf++;
  nflagRef.value = PREV_ERR_FAIL;
  cvRestore(mem, saved_t);

  // At maxnef failures or |h| = hmin, return failure
  if (Math.abs(mem.cv_h) <= mem.cv_hmin * ONEPSM ||
      nefRef.value === mem.cv_maxnef) {
    return CV_ERR_FAILURE;
  }

  // Set etamax = 1 to prevent step size increase
  mem.cv_etamax = 1.0;

  // Compute eta from dsm and rescale
  if (nefRef.value <= MXNEF1) {
    mem.cv_eta = 1.0 / (Math.pow(BIAS2 * dsm, 1.0 / mem.cv_L) + ADDON);
    mem.cv_eta = Math.max(mem.cv_eta_min_ef,
      Math.max(mem.cv_eta, mem.cv_hmin / Math.abs(mem.cv_h)));
    if (nefRef.value >= mem.cv_small_nef) {
      mem.cv_eta = Math.min(mem.cv_eta, mem.cv_eta_max_ef);
    }
    cvRescale(mem);
    return TRY_AGAIN;
  }

  // After MXNEF1 failures, force order reduction
  if (mem.cv_q > 1) {
    mem.cv_eta = Math.max(mem.cv_eta_min_ef,
      mem.cv_hmin / Math.abs(mem.cv_h));
    cvAdjustOrder(mem, -1);
    mem.cv_L = mem.cv_q;
    mem.cv_q--;
    mem.cv_qwait = mem.cv_L;
    cvRescale(mem);
    return TRY_AGAIN;
  }

  // Already at order 1: restart by reloading zn
  mem.cv_eta = Math.max(mem.cv_eta_min_ef,
    mem.cv_hmin / Math.abs(mem.cv_h));
  mem.cv_h *= mem.cv_eta;
  mem.cv_next_h = mem.cv_h;
  mem.cv_hscale = mem.cv_h;
  mem.cv_qwait = LONG_WAIT;
  mem.cv_nscon = 0;

  const retval = mem.cv_f!(mem.cv_tn, mem.cv_zn[0], mem.cv_tempv, mem.cv_user_data);
  mem.cv_nfe++;
  if (retval < 0) { return CV_RHSFUNC_FAIL; }
  if (retval > 0) { return CV_UNREC_RHSFUNC_ERR; }

  vScale(mem.cv_h, mem.cv_tempv, mem.cv_zn[1], mem.cv_N);

  return TRY_AGAIN;
}

/** cvCompleteStep - update state after successful step */
function cvCompleteStep(mem: CvodeMem): void {
  const N = mem.cv_N;

  mem.cv_nst++;
  mem.cv_nscon++;
  mem.cv_hu = mem.cv_h;
  mem.cv_qu = mem.cv_q;

  mem.first_step_after_resize = false;

  // Shift tau array
  for (let i = mem.cv_q; i >= 2; i--) {
    mem.cv_tau[i] = mem.cv_tau[i - 1];
  }
  if (mem.cv_q === 1 && mem.cv_nst > 1) {
    mem.cv_tau[2] = mem.cv_tau[1];
  }
  mem.cv_tau[1] = mem.cv_h;

  // Apply corrections to zn: zn[j] += l[j] * acor
  for (let j = 0; j <= mem.cv_q; j++) {
    const lj = mem.cv_l[j];
    const znj = mem.cv_zn[j];
    const acor = mem.cv_acor;
    for (let i = 0; i < N; i++) {
      znj[i] += lj * acor[i];
    }
  }

  // Decrement qwait, save acor if qwait==1 and q < qmax
  mem.cv_qwait--;
  if (mem.cv_qwait === 1 && mem.cv_q !== mem.cv_qmax) {
    mem.cv_zn[mem.cv_qmax].set(mem.cv_acor);
    mem.cv_saved_tq5 = mem.cv_tq[5];
    mem.cv_indx_acor = mem.cv_qmax;
  }
}

/** cvPrepareNextStep - set stepsize and order for next step */
function cvPrepareNextStep(mem: CvodeMem, dsm: number): void {
  if (mem.cv_etamax === 1.0) {
    // Defer step size or order changes
    mem.cv_qwait = Math.max(mem.cv_qwait, 2);
    mem.cv_qprime = mem.cv_q;
    mem.cv_hprime = mem.cv_h;
    mem.cv_eta = 1.0;
  } else {
    // Compute etaq
    mem.cv_etaq = 1.0 / (Math.pow(BIAS2 * dsm, 1.0 / mem.cv_L) + ADDON);

    if (mem.cv_qwait !== 0) {
      // No order change, just adjust eta
      mem.cv_eta = mem.cv_etaq;
      mem.cv_qprime = mem.cv_q;
      cvSetEta(mem);
    } else {
      // Consider order change
      mem.cv_qwait = 2;
      mem.cv_etaqm1 = cvComputeEtaqm1(mem);
      mem.cv_etaqp1 = cvComputeEtaqp1(mem);
      cvChooseEta(mem);
      cvSetEta(mem);
    }
  }
}

/** cvSetEta - adjust eta by heuristic limits and hmax */
function cvSetEta(mem: CvodeMem): void {
  if (mem.cv_eta > mem.cv_eta_min_fx && mem.cv_eta < mem.cv_eta_max_fx) {
    // Eta within fixed-step bounds, retain step size
    mem.cv_eta = 1.0;
    mem.cv_hprime = mem.cv_h;
  } else {
    if (mem.cv_eta >= mem.cv_eta_max_fx) {
      // Increase step size
      mem.cv_eta = Math.min(mem.cv_eta, mem.cv_etamax);
      mem.cv_eta /= Math.max(1.0,
        Math.abs(mem.cv_h) * mem.cv_hmax_inv * mem.cv_eta);
    } else {
      // Reduce step size
      mem.cv_eta = Math.max(mem.cv_eta, mem.cv_eta_min);
      mem.cv_eta = Math.max(mem.cv_eta,
        mem.cv_hmin / Math.abs(mem.cv_h));
    }
    mem.cv_hprime = mem.cv_h * mem.cv_eta;
    if (mem.cv_qprime < mem.cv_q) { mem.cv_nscon = 0; }
  }
}

/** cvComputeEtaqm1 - eta for possible order decrease */
function cvComputeEtaqm1(mem: CvodeMem): number {
  mem.cv_etaqm1 = 0.0;
  if (mem.cv_q > 1) {
    const ddn = wrmsNorm(mem.cv_N, mem.cv_zn[mem.cv_q], mem.cv_ewt) * mem.cv_tq[1];
    mem.cv_etaqm1 = 1.0 / (Math.pow(BIAS1 * ddn, 1.0 / mem.cv_q) + ADDON);
  }
  return mem.cv_etaqm1;
}

/** cvComputeEtaqp1 - eta for possible order increase */
function cvComputeEtaqp1(mem: CvodeMem): number {
  mem.cv_etaqp1 = 0.0;
  if (mem.cv_q !== mem.cv_qmax) {
    if (mem.cv_saved_tq5 === 0.0) { return mem.cv_etaqp1; }
    const cquot = (mem.cv_tq[5] / mem.cv_saved_tq5) *
      Math.pow(mem.cv_h / mem.cv_tau[2], mem.cv_L);

    // tempv = acor - cquot * zn[qmax]
    const N = mem.cv_N;
    const tempv = mem.cv_tempv;
    const znqmax = mem.cv_zn[mem.cv_qmax];
    const acor = mem.cv_acor;
    for (let i = 0; i < N; i++) {
      tempv[i] = acor[i] - cquot * znqmax[i];
    }

    const dup = wrmsNorm(N, tempv, mem.cv_ewt) * mem.cv_tq[3];
    mem.cv_etaqp1 = 1.0 / (Math.pow(BIAS3 * dup, 1.0 / (mem.cv_L + 1)) + ADDON);
  }
  return mem.cv_etaqp1;
}

/** cvChooseEta - choose max of etaqm1, etaq, etaqp1 */
function cvChooseEta(mem: CvodeMem): void {
  const etam = Math.max(mem.cv_etaqm1, Math.max(mem.cv_etaq, mem.cv_etaqp1));

  if (etam > mem.cv_eta_min_fx && etam < mem.cv_eta_max_fx) {
    mem.cv_eta = 1.0;
    mem.cv_qprime = mem.cv_q;
  } else {
    // Preference: keep > decrease > increase
    if (etam === mem.cv_etaq) {
      mem.cv_eta = mem.cv_etaq;
      mem.cv_qprime = mem.cv_q;
    } else if (etam === mem.cv_etaqm1) {
      mem.cv_eta = mem.cv_etaqm1;
      mem.cv_qprime = mem.cv_q - 1;
    } else {
      mem.cv_eta = mem.cv_etaqp1;
      mem.cv_qprime = mem.cv_q + 1;

      if (mem.cv_lmm === CV_BDF) {
        // Store Delta_n in zn[qmax] for order increase
        mem.cv_zn[mem.cv_qmax].set(mem.cv_acor);
      }
    }
  }
}

/** cvHandleFailure - set error and return flag */
function cvHandleFailure(mem: CvodeMem, flag: number): number {
  switch (flag) {
    case CV_ERR_FAILURE:
      mem.cv_error = `Error test failures at t=${mem.cv_tn}, h=${mem.cv_h}`;
      break;
    case CV_CONV_FAILURE:
      mem.cv_error = `Convergence failures at t=${mem.cv_tn}, h=${mem.cv_h}`;
      break;
    case CV_LSETUP_FAIL:
      mem.cv_error = `Linear solver setup failed at t=${mem.cv_tn}`;
      break;
    case CV_LSOLVE_FAIL:
      mem.cv_error = `Linear solver solve failed at t=${mem.cv_tn}`;
      break;
    case CV_RHSFUNC_FAIL:
      mem.cv_error = `RHS function failed at t=${mem.cv_tn}`;
      break;
    case CV_UNREC_RHSFUNC_ERR:
      mem.cv_error = `Unrecoverable RHS function error at t=${mem.cv_tn}`;
      break;
    case CV_REPTD_RHSFUNC_ERR:
      mem.cv_error = `Repeated RHS function errors at t=${mem.cv_tn}`;
      break;
    case CV_RTFUNC_FAIL:
      mem.cv_error = `Root function failed at t=${mem.cv_tn}`;
      break;
    case CV_TOO_CLOSE:
      mem.cv_error = 'tout too close to t0';
      break;
    default:
      mem.cv_error = `Unrecognized failure flag: ${flag}`;
      break;
  }
  return flag;
}

/** cvBDFStab - BDF stability limit detection (STALD) */
function cvBDFStab(mem: CvodeMem): void {
  // If order >= 3, save scaled derivative data
  if (mem.cv_q >= 3) {
    for (let k = 1; k <= 3; k++) {
      for (let i = 5; i >= 2; i--) {
        mem.cv_ssdat[i][k] = mem.cv_ssdat[i - 1][k];
      }
    }

    let factorial = 1;
    for (let i = 1; i <= mem.cv_q - 1; i++) { factorial *= i; }

    const sq = factorial * mem.cv_q * (mem.cv_q + 1) * mem.cv_acnrm /
      Math.max(mem.cv_tq[5], TINY);
    const sqm1 = factorial * mem.cv_q *
      wrmsNorm(mem.cv_N, mem.cv_zn[mem.cv_q], mem.cv_ewt);
    const sqm2 = factorial *
      wrmsNorm(mem.cv_N, mem.cv_zn[mem.cv_q - 1], mem.cv_ewt);

    mem.cv_ssdat[1][1] = sqm2 * sqm2;
    mem.cv_ssdat[1][2] = sqm1 * sqm1;
    mem.cv_ssdat[1][3] = sq * sq;
  }

  if (mem.cv_qprime >= mem.cv_q) {
    // If enough data, call stability limit detection
    if (mem.cv_q >= 3 && mem.cv_nscon >= mem.cv_q + 5) {
      const ldflag = cvSLdet(mem);
      if (ldflag > 3) {
        // Stability limit violation: reduce order
        mem.cv_qprime = mem.cv_q - 1;
        mem.cv_eta = mem.cv_etaqm1;
        mem.cv_eta = Math.min(mem.cv_eta, mem.cv_etamax);
        mem.cv_eta = mem.cv_eta /
          Math.max(1.0, Math.abs(mem.cv_h) * mem.cv_hmax_inv * mem.cv_eta);
        mem.cv_hprime = mem.cv_h * mem.cv_eta;
        mem.cv_nor = mem.cv_nor + 1;
      }
    }
  } else {
    // Let order increase happen, reset nscon
    mem.cv_nscon = 0;
  }
}

/** cvSLdet - Stability limit detection algorithm (~300 lines from C) */
function cvSLdet(mem: CvodeMem): number {
  let kmin = 0;
  let kflag = 0;

  // Allocate local arrays with 1-based indexing
  const rat: number[][] = [];
  for (let i = 0; i < 5; i++) rat[i] = [0, 0, 0, 0];
  const rav = [0, 0, 0, 0];
  const qkr = [0, 0, 0, 0];
  const sigsq = [0, 0, 0, 0];
  const smax = [0, 0, 0, 0];
  const ssmax = [0, 0, 0, 0];
  const drr = [0, 0, 0, 0];
  const rrc = [0, 0, 0, 0];
  const sqmx = [0, 0, 0, 0];
  const vrat = [0, 0, 0, 0, 0];

  // qjk[j][k], qc[i][k], qco[i][k] - all 1-based
  const qjk: number[][] = [];
  for (let j = 0; j < 4; j++) qjk[j] = [0, 0, 0, 0];
  const qc: number[][] = [];
  for (let i = 0; i < 6; i++) qc[i] = [0, 0, 0, 0];
  const qco: number[][] = [];
  for (let i = 0; i < 6; i++) qco[i] = [0, 0, 0, 0];

  // Cutoffs and tolerances
  const rrcut = 0.98;
  const vrrtol = 1.0e-4;
  const vrrt2 = 5.0e-4;
  const sqtol = 1.0e-3;
  const rrtol = 1.0e-2;

  let rr = 0.0;

  // Get maxima, minima, and variances, form quartic coefficients
  for (let k = 1; k <= 3; k++) {
    let smink = mem.cv_ssdat[1][k];
    let smaxk = 0.0;

    for (let i = 1; i <= 5; i++) {
      smink = Math.min(smink, mem.cv_ssdat[i][k]);
      smaxk = Math.max(smaxk, mem.cv_ssdat[i][k]);
    }

    if (smink < TINY * smaxk) { return -1; }

    smax[k] = smaxk;
    ssmax[k] = smaxk * smaxk;

    let sumrat = 0.0;
    let sumrsq = 0.0;
    for (let i = 1; i <= 4; i++) {
      rat[i][k] = mem.cv_ssdat[i][k] / mem.cv_ssdat[i + 1][k];
      sumrat += rat[i][k];
      sumrsq += rat[i][k] * rat[i][k];
    }
    rav[k] = 0.25 * sumrat;
    vrat[k] = Math.abs(0.25 * sumrsq - rav[k] * rav[k]);

    qc[5][k] = mem.cv_ssdat[1][k] * mem.cv_ssdat[3][k] -
               mem.cv_ssdat[2][k] * mem.cv_ssdat[2][k];
    qc[4][k] = mem.cv_ssdat[2][k] * mem.cv_ssdat[3][k] -
               mem.cv_ssdat[1][k] * mem.cv_ssdat[4][k];
    qc[3][k] = 0.0;
    qc[2][k] = mem.cv_ssdat[2][k] * mem.cv_ssdat[5][k] -
               mem.cv_ssdat[3][k] * mem.cv_ssdat[4][k];
    qc[1][k] = mem.cv_ssdat[4][k] * mem.cv_ssdat[4][k] -
               mem.cv_ssdat[3][k] * mem.cv_ssdat[5][k];

    for (let i = 1; i <= 5; i++) { qco[i][k] = qc[i][k]; }
  }

  // Isolate normal or nearly-normal matrix case
  const vmin = Math.min(vrat[1], Math.min(vrat[2], vrat[3]));
  const vmax = Math.max(vrat[1], Math.max(vrat[2], vrat[3]));

  if (vmin < vrrtol * vrrtol) {
    if (vmax > vrrt2 * vrrt2) { return -2; }
    else {
      rr = (rav[1] + rav[2] + rav[3]) / 3.0;
      let drrmax = 0.0;
      for (let k = 1; k <= 3; k++) {
        const adrr = Math.abs(rav[k] - rr);
        drrmax = Math.max(drrmax, adrr);
      }
      if (drrmax > vrrt2) { return -3; }
      kflag = 1;
    }
  } else {
    // Use quartics to get rr
    if (Math.abs(qco[1][1]) < TINY * ssmax[1]) { return -4; }

    let tem = qco[1][2] / qco[1][1];
    for (let i = 2; i <= 5; i++) { qco[i][2] = qco[i][2] - tem * qco[i][1]; }
    qco[1][2] = 0.0;

    tem = qco[1][3] / qco[1][1];
    for (let i = 2; i <= 5; i++) { qco[i][3] = qco[i][3] - tem * qco[i][1]; }
    qco[1][3] = 0.0;

    if (Math.abs(qco[2][2]) < TINY * ssmax[2]) { return -4; }

    tem = qco[2][3] / qco[2][2];
    for (let i = 3; i <= 5; i++) { qco[i][3] = qco[i][3] - tem * qco[i][2]; }

    if (Math.abs(qco[4][3]) < TINY * ssmax[3]) { return -4; }

    rr = -qco[5][3] / qco[4][3];

    if (rr < TINY || rr > 100.0) { return -5; }

    for (let k = 1; k <= 3; k++) {
      qkr[k] = qc[5][k] + rr * (qc[4][k] + rr * rr * (qc[2][k] + rr * qc[1][k]));
    }

    let sqmax = 0.0;
    for (let k = 1; k <= 3; k++) {
      const saqk = Math.abs(qkr[k]) / ssmax[k];
      if (saqk > sqmax) { sqmax = saqk; }
    }

    if (sqmax < sqtol) {
      kflag = 2;
    } else {
      // Newton corrections to improve rr
      for (let it = 1; it <= 3; it++) {
        for (let k = 1; k <= 3; k++) {
          const qp = qc[4][k] + rr * rr * (3.0 * qc[2][k] + rr * 4.0 * qc[1][k]);
          drr[k] = 0.0;
          if (Math.abs(qp) > TINY * ssmax[k]) { drr[k] = -qkr[k] / qp; }
          rrc[k] = rr + drr[k];
        }

        for (let k = 1; k <= 3; k++) {
          const s = rrc[k];
          let sqmaxk = 0.0;
          for (let j = 1; j <= 3; j++) {
            qjk[j][k] = qc[5][j] + s * (qc[4][j] + s * s * (qc[2][j] + s * qc[1][j]));
            const saqj = Math.abs(qjk[j][k]) / ssmax[j];
            if (saqj > sqmaxk) { sqmaxk = saqj; }
          }
          sqmx[k] = sqmaxk;
        }

        let sqmin = sqmx[1] + 1.0;
        for (let k = 1; k <= 3; k++) {
          if (sqmx[k] < sqmin) {
            kmin = k;
            sqmin = sqmx[k];
          }
        }
        rr = rrc[kmin];

        if (sqmin < sqtol) {
          kflag = 3;
          break;
        } else {
          for (let j = 1; j <= 3; j++) { qkr[j] = qjk[j][kmin]; }
        }
      }

      if (kflag !== 3) {
        // sqmin > sqtol after all Newton iterations
        return -6;
      }
    }
  }

  // Given rr, find sigsq[k] and verify rr
  for (let k = 1; k <= 3; k++) {
    const rsa = mem.cv_ssdat[1][k];
    const rsb = mem.cv_ssdat[2][k] * rr;
    const rsc = mem.cv_ssdat[3][k] * rr * rr;
    const rsd = mem.cv_ssdat[4][k] * rr * rr * rr;
    const rd1a = rsa - rsb;
    const rd1b = rsb - rsc;
    const rd1c = rsc - rsd;
    const rd2a = rd1a - rd1b;
    const rd2b = rd1b - rd1c;
    const rd3a = rd2a - rd2b;

    if (Math.abs(rd1b) < TINY * smax[k]) { return -7; }

    const cest1 = -rd3a / rd1b;
    if (cest1 < TINY || cest1 > 4.0) { return -7; }

    const corr1 = (rd2b / cest1) / (rr * rr);
    sigsq[k] = mem.cv_ssdat[3][k] + corr1;
  }

  if (sigsq[2] < TINY) { return -8; }

  const ratp = sigsq[3] / sigsq[2];
  const ratm = sigsq[1] / sigsq[2];
  const qfac1 = 0.25 * (mem.cv_q * mem.cv_q - 1);
  const qfac2 = 2.0 / (mem.cv_q - 1);
  const bb = ratp * ratm - 1.0 - qfac1 * ratp;
  const tem2 = 1.0 - qfac2 * bb;

  if (Math.abs(tem2) < TINY) { return -8; }

  const rrb = 1.0 / tem2;

  if (Math.abs(rrb - rr) > rrtol) { return -9; }

  // Check if rr is above cutoff
  if (rr > rrcut) {
    if (kflag === 1) { kflag = 4; }
    if (kflag === 2) { kflag = 5; }
    if (kflag === 3) { kflag = 6; }
  }

  return kflag;
}

// Rootfinding functions imported from cvode_root.ts (cvRcheck1, cvRcheck2, cvRcheck3)
