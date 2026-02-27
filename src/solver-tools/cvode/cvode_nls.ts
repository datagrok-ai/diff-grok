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

import {CvodeMem, wrmsNorm, vLinearSum, vScale, vConst,
  CV_SUCCESS, CV_BDF, CV_ADAMS, CV_RHSFUNC_FAIL, CV_LSETUP_FAIL, CV_LSOLVE_FAIL,
  RHSFUNC_RECVR, SUN_NLS_CONV_RECVR, SUN_NLS_CONTINUE,
  CV_NO_FAILURES, CV_FAIL_BAD_J, CV_FAIL_OTHER,
  FIRST_CALL, PREV_CONV_FAIL, PREV_ERR_FAIL, DO_ERROR_TEST,
  NLS_MAXCOR, CRDOWN, RDIV, ONEPSM} from './common';

/**
 * cvNls - Solve the nonlinear system for one step.
 *
 * For BDF: Newton iteration solving G(ycor) = 0 where
 *   G(ycor) = rl1*zn[1] + ycor - gamma*f(tn, zn[0]+ycor)
 * The Newton update is: (I - gamma*J) * delta = -G(ycor)
 *
 * For Adams: Fixed-point iteration
 *   ycor_new = h * rl1 * f(tn, zn[0]+ycor) - rl1*zn[1]
 */
export function cvNls(mem: CvodeMem, nflag: number): number {
  const N = mem.cv_N;
  let callSetup: boolean;

  // Decide whether to call linear solver setup
  if (mem.cv_lsetup) {
    mem.convfail = (nflag === FIRST_CALL || nflag === PREV_ERR_FAIL) ?
      CV_NO_FAILURES :
      CV_FAIL_OTHER;

    callSetup = (nflag === PREV_CONV_FAIL) || (nflag === PREV_ERR_FAIL) ||
                (mem.cv_nst === 0) || (mem.first_step_after_resize) ||
                (mem.cv_nst >= mem.cv_nstlp + mem.cv_msbp) ||
                (Math.abs(mem.cv_gamrat - 1.0) > mem.cv_dgmax_lsetup);
  } else {
    mem.cv_crate = 1.0;
    callSetup = false;
  }

  // Initial guess for correction: acor = 0
  vConst(0.0, mem.cv_acor, N);

  // Set acnrmcur to false
  mem.cv_acnrmcur = false;

  // Dispatch to Newton or fixed-point
  let flag: number;
  if (mem.cv_lmm === CV_BDF)
    flag = cvNewtonIteration(mem, callSetup);
  else
    flag = cvFixedPointIteration(mem);


  if (flag !== CV_SUCCESS) return flag;

  // Solve successful - update y = zn[0] + acor
  vLinearSum(1.0, mem.cv_zn[0], 1.0, mem.cv_acor, mem.cv_y, N);

  // Compute acnrm if not already done
  if (!mem.cv_acnrmcur)
    mem.cv_acnrm = wrmsNorm(N, mem.cv_acor, mem.cv_ewt);


  // Update Jacobian status
  mem.cv_jcur = false;

  return CV_SUCCESS;
}

/**
 * Newton iteration for BDF method.
 * Solves: (I - gamma*J) * delta = -G(ycor)
 * where G(ycor) = rl1*zn[1] + ycor - gamma*f(tn, zn[0]+ycor)
 */
function cvNewtonIteration(mem: CvodeMem, callSetup: boolean): number {
  const N = mem.cv_N;
  let retval: number;

  // Call linear solver setup if needed
  if (callSetup && mem.cv_lsetup) {
    // y = zn[0] + acor (current approximation)
    vLinearSum(1.0, mem.cv_zn[0], 1.0, mem.cv_acor, mem.cv_y, N);

    // Evaluate f at current approximation
    retval = mem.cv_f!(mem.cv_tn, mem.cv_y, mem.cv_ftemp, mem.cv_user_data);
    mem.cv_nfe++;
    if (retval < 0) return CV_RHSFUNC_FAIL;
    if (retval > 0) return RHSFUNC_RECVR;

    // Call lsetup
    const jcurPtr = {value: false};
    retval = mem.cv_lsetup!(mem, mem.convfail, mem.cv_y, mem.cv_ftemp, jcurPtr,
      mem.cv_vtemp1, mem.cv_vtemp2, mem.cv_vtemp3);
    mem.cv_nsetups++;

    mem.cv_jcur = jcurPtr.value;
    mem.cv_gamrat = 1.0;
    mem.cv_gammap = mem.cv_gamma;
    mem.cv_crate = 1.0;
    mem.cv_nstlp = mem.cv_nst;

    if (retval < 0) return CV_LSETUP_FAIL;
    if (retval > 0) return SUN_NLS_CONV_RECVR;
  }

  // Newton iteration loop
  for (let m = 0; m < NLS_MAXCOR; m++) {
    mem.cv_nni++;

    // y = zn[0] + acor
    vLinearSum(1.0, mem.cv_zn[0], 1.0, mem.cv_acor, mem.cv_y, N);

    // Evaluate f(tn, y)
    retval = mem.cv_f!(mem.cv_tn, mem.cv_y, mem.cv_ftemp, mem.cv_user_data);
    mem.cv_nfe++;
    if (retval < 0) return CV_RHSFUNC_FAIL;
    if (retval > 0) return RHSFUNC_RECVR;

    // Compute residual: b = rl1*zn[1] + acor - gamma*ftemp
    // Store in tempv as the RHS for the linear solve
    const b = mem.cv_tempv;
    for (let i = 0; i < N; i++)
      b[i] = mem.cv_rl1 * mem.cv_zn[1][i] + mem.cv_acor[i] - mem.cv_gamma * mem.cv_ftemp[i];


    // Solve: (I - gamma*J) * delta = -b, so pass b (which is G, so we negate)
    // Actually the convention is lsolve solves (I-gamma*J)*x = b
    // and b here is -G so we negate
    for (let i = 0; i < N; i++) b[i] = -b[i];

    retval = mem.cv_lsolve!(mem, b, mem.cv_ewt, mem.cv_y, mem.cv_ftemp);
    if (retval < 0) return CV_LSOLVE_FAIL;
    if (retval > 0) return SUN_NLS_CONV_RECVR;

    // Update acor += delta (delta is in b after solve)
    for (let i = 0; i < N; i++) mem.cv_acor[i] += b[i];

    // Convergence test
    const del = wrmsNorm(N, b, mem.cv_ewt);
    const convTestResult = cvNlsConvTest(mem, del, m, mem.cv_tq[4]);

    if (convTestResult === CV_SUCCESS) return CV_SUCCESS;
    if (convTestResult === SUN_NLS_CONV_RECVR) {
      mem.cv_nnf++;
      return SUN_NLS_CONV_RECVR;
    }
    // SUN_NLS_CONTINUE: continue iterating
  }

  // Max iterations reached without convergence
  mem.cv_nnf++;
  return SUN_NLS_CONV_RECVR;
}

/**
 * Fixed-point iteration for Adams method.
 * ycor = rl1 * (h * f(tn, zn[0]+ycor) - zn[1])
 */
function cvFixedPointIteration(mem: CvodeMem): number {
  const N = mem.cv_N;
  let retval: number;

  for (let m = 0; m < NLS_MAXCOR; m++) {
    mem.cv_nni++;

    // y = zn[0] + acor
    vLinearSum(1.0, mem.cv_zn[0], 1.0, mem.cv_acor, mem.cv_y, N);

    // Evaluate f(tn, y) -> store in tempv
    retval = mem.cv_f!(mem.cv_tn, mem.cv_y, mem.cv_tempv, mem.cv_user_data);
    mem.cv_nfe++;
    if (retval < 0) return CV_RHSFUNC_FAIL;
    if (retval > 0) return RHSFUNC_RECVR;

    // Compute new ycor: rl1 * (h * f - zn[1])
    // First: tempv = h * tempv - zn[1]
    for (let i = 0; i < N; i++)
      mem.cv_tempv[i] = mem.cv_h * mem.cv_tempv[i] - mem.cv_zn[1][i];

    // Then scale by rl1
    const newAcor = mem.cv_vtemp1;
    vScale(mem.cv_rl1, mem.cv_tempv, newAcor, N);

    // Compute delta = newAcor - acor
    const delta = mem.cv_tempv;
    for (let i = 0; i < N; i++)
      delta[i] = newAcor[i] - mem.cv_acor[i];


    // Update acor
    for (let i = 0; i < N; i++) mem.cv_acor[i] = newAcor[i];

    // Convergence test
    const del = wrmsNorm(N, delta, mem.cv_ewt);
    const convTestResult = cvNlsConvTest(mem, del, m, mem.cv_tq[4]);

    if (convTestResult === CV_SUCCESS) return CV_SUCCESS;
    if (convTestResult === SUN_NLS_CONV_RECVR) {
      mem.cv_nnf++;
      return SUN_NLS_CONV_RECVR;
    }
    // SUN_NLS_CONTINUE: continue iterating
  }

  // Max iterations reached
  mem.cv_nnf++;
  return SUN_NLS_CONV_RECVR;
}

/**
 * cvNlsConvTest - Test for convergence of the nonlinear iteration.
 *
 * del: WRMS norm of the current update
 * m: iteration count (0-based)
 * tol: tolerance (tq[4])
 *
 * Returns: CV_SUCCESS if converged, SUN_NLS_CONV_RECVR if diverging,
 *          SUN_NLS_CONTINUE if not yet converged
 */
function cvNlsConvTest(mem: CvodeMem, del: number, m: number, tol: number): number {
  // Update convergence rate estimate
  if (m > 0)
    mem.cv_crate = Math.max(CRDOWN * mem.cv_crate, del / mem.cv_delp);


  // Test for convergence
  const dcon = del * Math.min(1.0, mem.cv_crate) / tol;

  if (dcon <= 1.0) {
    // Converged
    if (m === 0)
      mem.cv_acnrm = del;
    else
      mem.cv_acnrm = wrmsNorm(mem.cv_N, mem.cv_acor, mem.cv_ewt);

    mem.cv_acnrmcur = true;
    return CV_SUCCESS;
  }

  // Check for divergence
  if (m >= 1 && del > RDIV * mem.cv_delp)
    return SUN_NLS_CONV_RECVR;


  // Save del and continue
  mem.cv_delp = del;
  return SUN_NLS_CONTINUE;
}
