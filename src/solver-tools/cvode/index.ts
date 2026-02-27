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

// CVODE TypeScript port - Barrel exports
// Public API for the CVODE solver module

// High-level class API
export { Cvode } from './cvode_class';
export type { OdeFunction, CvodeOptions, CvodeSolveResult } from './cvode_class';

// Low-level procedural API
export { cvodeCreate, cvodeInit, cvode } from './cvode';
export { cvodeGetDky } from './common';
export { cvodeSetLinearSolver, cvodeSetJacFn } from './cvode_ls';
export { cvodeRootInit } from './cvode_root';

// Configuration and statistics
export type { CvodeStats } from './cvode_io';
export {
  cvodeSStolerances,
  cvodeSVtolerances,
  cvodeSetMaxNumSteps,
  cvodeSetMaxOrd,
  cvodeSetMaxStep,
  cvodeSetMinStep,
  cvodeSetInitStep,
  cvodeSetStopTime,
  cvodeSetUserData,
  cvodeSetStabLimDet,
  cvodeGetNumSteps,
  cvodeGetNumRhsEvals,
  cvodeGetNumLinSolvSetups,
  cvodeGetNumErrTestFails,
  cvodeGetNumNonlinSolvIters,
  cvodeGetNumNonlinSolvConvFails,
  cvodeGetNumJacEvals,
  cvodeGetNumGEvals,
  cvodeGetLastOrder,
  cvodeGetLastStep,
  cvodeGetCurrentTime,
  cvodeGetIntegratorStats,
} from './cvode_io';

// Types and constants
export type { CvodeRhsFn, CvodeJacFn, CvodeRootFn } from './common';
export { CvodeMem, CvLsMem } from './common';
export {
  CV_ADAMS, CV_BDF,
  CV_NORMAL, CV_ONE_STEP,
  CV_SUCCESS, CV_TSTOP_RETURN, CV_ROOT_RETURN,
  CV_TOO_MUCH_WORK, CV_TOO_MUCH_ACC, CV_ERR_FAILURE,
  CV_CONV_FAILURE, CV_RHSFUNC_FAIL,
} from './common';
