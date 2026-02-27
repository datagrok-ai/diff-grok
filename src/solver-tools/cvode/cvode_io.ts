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

// CVODE TypeScript port - Configuration and statistics I/O functions
// Ported from SUNDIALS CVODE v7.5.0 (cvode_io.c)

import {
  CvodeMem,
  CV_SUCCESS,
  CV_SS,
  CV_SV,
  ADAMS_Q_MAX,
  BDF_Q_MAX,
  CV_ADAMS,
} from './common';

// ============================================================================
// Statistics interface
// ============================================================================

export interface CvodeStats {
  nSteps: number;
  nRhsEvals: number;
  nLinSolvSetups: number;
  nErrTestFails: number;
  nNonlinSolvIters: number;
  nNonlinSolvConvFails: number;
  nJacEvals: number;
  nGEvals: number;
  lastOrder: number;
  lastStep: number;
  currentTime: number;
}

// ============================================================================
// Tolerance functions
// ============================================================================

/**
 * Set scalar relative and absolute tolerances.
 * Requires rtol >= 0 and atol >= 0.
 */
export function cvodeSStolerances(mem: CvodeMem, rtol: number, atol: number): number {
  if (rtol < 0) {
    mem.cv_error = 'cvodeSStolerances: rtol must be >= 0';
    return -1;
  }
  if (atol < 0) {
    mem.cv_error = 'cvodeSStolerances: atol must be >= 0';
    return -1;
  }

  mem.cv_reltol = rtol;
  mem.cv_Sabstol = atol;
  mem.cv_itol = CV_SS;
  mem.cv_atolmin0 = (atol === 0);

  return CV_SUCCESS;
}

/**
 * Set scalar relative tolerance and vector absolute tolerance.
 * Requires rtol >= 0 and all atol[i] >= 0.
 */
export function cvodeSVtolerances(mem: CvodeMem, rtol: number, atol: Float64Array): number {
  if (rtol < 0) {
    mem.cv_error = 'cvodeSVtolerances: rtol must be >= 0';
    return -1;
  }
  for (let i = 0; i < atol.length; i++) {
    if (atol[i] < 0) {
      mem.cv_error = `cvodeSVtolerances: atol[${i}] must be >= 0`;
      return -1;
    }
  }

  mem.cv_reltol = rtol;
  mem.cv_Vabstol = new Float64Array(atol);
  mem.cv_itol = CV_SV;

  let minAtol = Infinity;
  for (let i = 0; i < atol.length; i++) {
    if (atol[i] < minAtol) minAtol = atol[i];
  }
  mem.cv_atolmin0 = (minAtol === 0);
  mem.cv_VabstolMallocDone = true;

  return CV_SUCCESS;
}

// ============================================================================
// Optional input functions
// ============================================================================

/**
 * Set the maximum number of internal steps per call to CVode.
 */
export function cvodeSetMaxNumSteps(mem: CvodeMem, mxstep: number): number {
  mem.cv_mxstep = mxstep;
  return CV_SUCCESS;
}

/**
 * Set the maximum method order. Clamped to [1, qmax] where qmax is
 * ADAMS_Q_MAX or BDF_Q_MAX depending on the linear multistep method.
 */
export function cvodeSetMaxOrd(mem: CvodeMem, maxord: number): number {
  const qmax = (mem.cv_lmm === CV_ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;
  if (maxord < 1) maxord = 1;
  if (maxord > qmax) maxord = qmax;
  mem.cv_qmax = maxord;
  return CV_SUCCESS;
}

/**
 * Set the maximum step size. If hmax <= 0, no upper limit is imposed
 * (cv_hmax_inv = 0). Otherwise cv_hmax_inv = 1 / hmax.
 */
export function cvodeSetMaxStep(mem: CvodeMem, hmax: number): number {
  if (hmax <= 0) {
    mem.cv_hmax_inv = 0;
  } else {
    mem.cv_hmax_inv = 1.0 / hmax;
  }
  return CV_SUCCESS;
}

/**
 * Set the minimum step size (stored as absolute value).
 */
export function cvodeSetMinStep(mem: CvodeMem, hmin: number): number {
  mem.cv_hmin = Math.abs(hmin);
  return CV_SUCCESS;
}

/**
 * Set the initial step size hint.
 */
export function cvodeSetInitStep(mem: CvodeMem, hin: number): number {
  mem.cv_hin = hin;
  return CV_SUCCESS;
}

/**
 * Set a stop time. Integration will not proceed past tstop.
 */
export function cvodeSetStopTime(mem: CvodeMem, tstop: number): number {
  mem.cv_tstop = tstop;
  mem.cv_tstopset = true;
  mem.cv_tstopinterp = true;
  return CV_SUCCESS;
}

/**
 * Attach user data that will be passed to the RHS, Jacobian, and root
 * functions.
 */
export function cvodeSetUserData(mem: CvodeMem, data: any): number {
  mem.cv_user_data = data;
  return CV_SUCCESS;
}

/**
 * Enable or disable BDF stability limit detection.
 */
export function cvodeSetStabLimDet(mem: CvodeMem, on: boolean): number {
  mem.cv_sldeton = on;
  return CV_SUCCESS;
}

// ============================================================================
// Statistics getters
// ============================================================================

/** Return the total number of internal steps taken. */
export function cvodeGetNumSteps(mem: CvodeMem): number {
  return mem.cv_nst;
}

/** Return the total number of RHS evaluations. */
export function cvodeGetNumRhsEvals(mem: CvodeMem): number {
  return mem.cv_nfe;
}

/** Return the number of linear solver setup calls. */
export function cvodeGetNumLinSolvSetups(mem: CvodeMem): number {
  return mem.cv_nsetups;
}

/** Return the number of error test failures. */
export function cvodeGetNumErrTestFails(mem: CvodeMem): number {
  return mem.cv_netf;
}

/** Return the number of nonlinear solver iterations. */
export function cvodeGetNumNonlinSolvIters(mem: CvodeMem): number {
  return mem.cv_nni;
}

/** Return the number of nonlinear solver convergence failures. */
export function cvodeGetNumNonlinSolvConvFails(mem: CvodeMem): number {
  return mem.cv_ncfn;
}

/** Return the number of Jacobian evaluations (from the linear solver). */
export function cvodeGetNumJacEvals(mem: CvodeMem): number {
  return mem.cv_lmem?.nje ?? 0;
}

/** Return the total number of root-function evaluations. */
export function cvodeGetNumGEvals(mem: CvodeMem): number {
  return mem.cv_nge;
}

/** Return the method order used on the last successful step. */
export function cvodeGetLastOrder(mem: CvodeMem): number {
  return mem.cv_qu;
}

/** Return the step size used on the last successful step. */
export function cvodeGetLastStep(mem: CvodeMem): number {
  return mem.cv_hu;
}

/** Return the current internal time reached by the solver. */
export function cvodeGetCurrentTime(mem: CvodeMem): number {
  return mem.cv_tn;
}

/**
 * Return a snapshot of all integrator statistics as a plain object.
 */
export function cvodeGetIntegratorStats(mem: CvodeMem): CvodeStats {
  return {
    nSteps: cvodeGetNumSteps(mem),
    nRhsEvals: cvodeGetNumRhsEvals(mem),
    nLinSolvSetups: cvodeGetNumLinSolvSetups(mem),
    nErrTestFails: cvodeGetNumErrTestFails(mem),
    nNonlinSolvIters: cvodeGetNumNonlinSolvIters(mem),
    nNonlinSolvConvFails: cvodeGetNumNonlinSolvConvFails(mem),
    nJacEvals: cvodeGetNumJacEvals(mem),
    nGEvals: cvodeGetNumGEvals(mem),
    lastOrder: cvodeGetLastOrder(mem),
    lastStep: cvodeGetLastStep(mem),
    currentTime: cvodeGetCurrentTime(mem),
  };
}
