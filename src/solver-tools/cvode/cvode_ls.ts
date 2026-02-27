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

// CVODE TypeScript port - Dense linear solver interface
// Translated from SUNDIALS cvode_ls.c
//
// Instead of SUNDIALS' SUNLinearSolver / SUNMatrix abstraction, we use
// direct Float64Array[] matrices and dgefa / dgesl from dense_linalg.ts.

import { CvodeMem, CvLsMem, CvodeJacFn, wrmsNorm,
         CV_SUCCESS, CV_RHSFUNC_FAIL, RHSFUNC_RECVR,
         CV_NO_FAILURES, CV_FAIL_BAD_J, CV_FAIL_OTHER,
         CVLS_MSBJ, CVLS_DGMAX, ONEPSM } from './common';
import { dgefa, dgesl } from './dense_linalg';

/**
 * cvodeSetLinearSolver - Initialize the dense linear solver.
 * Sets up CvLsMem, allocates matrix, sets function pointers on CvodeMem.
 */
export function cvodeSetLinearSolver(mem: CvodeMem): number {
  const N = mem.cv_N;
  const lmem = new CvLsMem();

  // Allocate dense matrix A (N x N, row-major)
  lmem.A = new Array(N);
  lmem.savedJ = new Array(N);
  for (let i = 0; i < N; i++) {
    lmem.A[i] = new Float64Array(N);
    lmem.savedJ[i] = new Float64Array(N);
  }
  lmem.pivots = new Int32Array(N);

  // Set defaults
  lmem.jacDQ = true;
  lmem.jacFn = null;
  lmem.jbad = true;
  lmem.nje = 0;
  lmem.nfeDQ = 0;
  lmem.nstlj = 0;
  lmem.tnlj = 0;
  lmem.msbj = CVLS_MSBJ;
  lmem.dgmax_jbad = CVLS_DGMAX;

  // Attach to CvodeMem
  mem.cv_lmem = lmem;
  mem.cv_linit = cvLsInitialize;
  mem.cv_lsetup = cvLsSetup;
  mem.cv_lsolve = cvLsSolve;
  mem.cv_lfree = cvLsFree;

  return CV_SUCCESS;
}

/**
 * cvodeSetJacFn - Set a user-supplied Jacobian function.
 */
export function cvodeSetJacFn(mem: CvodeMem, jac: CvodeJacFn | null): number {
  const lmem = mem.cv_lmem;
  if (!lmem) return -1;

  if (jac) {
    lmem.jacDQ = false;
    lmem.jacFn = jac;
  } else {
    lmem.jacDQ = true;
    lmem.jacFn = null;
  }
  return CV_SUCCESS;
}

/**
 * cvLsInitialize - Initialize linear solver counters.
 */
function cvLsInitialize(mem: CvodeMem): number {
  const lmem = mem.cv_lmem!;
  lmem.nje = 0;
  lmem.nfeDQ = 0;
  lmem.nstlj = 0;
  lmem.jbad = true;
  return CV_SUCCESS;
}

/**
 * cvLsSetup - Compute Jacobian and form/factor the linear system matrix A = I - gamma*J.
 *
 * Heuristic for when to recompute Jacobian:
 * - On first call (nstlj == 0)
 * - If convfail is CV_FAIL_BAD_J (NLS flagged Jacobian as bad)
 * - If enough steps since last eval (nst >= nstlj + msbj)
 * - If gamma changed significantly (|gamrat - 1| > dgmax_jbad)
 */
function cvLsSetup(mem: CvodeMem, convfail: number, ypred: Float64Array,
                   fpred: Float64Array, jcurPtr: { value: boolean },
                   tmp1: Float64Array, tmp2: Float64Array,
                   tmp3: Float64Array): number {
  const lmem = mem.cv_lmem!;
  const N = mem.cv_N;
  let retval: number;

  // Determine if Jacobian needs updating
  const jbad = lmem.jbad;
  let jok = !jbad;

  if (!jbad) {
    // Check gamma ratio
    if (Math.abs(mem.cv_gamrat - 1.0) > lmem.dgmax_jbad) {
      jok = false;
    }
    // Check steps since last Jacobian
    if (mem.cv_nst >= lmem.nstlj + lmem.msbj) {
      jok = false;
    }
    // Check convfail
    if (convfail === CV_FAIL_BAD_J) {
      jok = false;
    }
  }

  if (jok) {
    // Reuse saved Jacobian, just form A = I - gamma*J
    jcurPtr.value = false;

    // Copy savedJ to A and form I - gamma*J
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        lmem.A[i][j] = -mem.cv_gamma * lmem.savedJ[i][j];
      }
      lmem.A[i][i] += 1.0;
    }
  } else {
    // Recompute Jacobian
    lmem.nje++;
    lmem.nstlj = mem.cv_nst;
    lmem.jbad = false;
    jcurPtr.value = true;

    // Evaluate Jacobian
    if (lmem.jacDQ) {
      retval = cvLsDQJac(mem, mem.cv_tn, ypred, fpred, lmem.savedJ, tmp1, tmp2);
    } else {
      retval = lmem.jacFn!(mem.cv_tn, ypred, fpred, lmem.savedJ, mem.cv_user_data);
    }

    if (retval < 0) {
      lmem.jbad = true;
      return -1;
    }
    if (retval > 0) {
      lmem.jbad = true;
      return 1;
    }

    // Form A = I - gamma*J
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        lmem.A[i][j] = -mem.cv_gamma * lmem.savedJ[i][j];
      }
      lmem.A[i][i] += 1.0;
    }
  }

  // Factor A using LU decomposition
  const info = dgefa(lmem.A, N, lmem.pivots);
  if (info !== 0) {
    lmem.jbad = true;
    return 1; // recoverable
  }

  return CV_SUCCESS;
}

/**
 * cvLsSolve - Solve the linear system A*x = b using LU factors.
 */
function cvLsSolve(mem: CvodeMem, b: Float64Array, weight: Float64Array,
                   ycur: Float64Array, fcur: Float64Array): number {
  const lmem = mem.cv_lmem!;
  const N = mem.cv_N;

  dgesl(lmem.A, N, lmem.pivots, b);

  // Apply gamma/gammap correction if gamma has changed since last setup
  if (mem.cv_gamrat !== 1.0) {
    const factor = 2.0 / (1.0 + mem.cv_gamrat);
    for (let i = 0; i < N; i++) {
      b[i] *= factor;
    }
  }

  return CV_SUCCESS;
}

/**
 * cvLsDQJac - Compute difference-quotient Jacobian approximation.
 *
 * For each column j:
 *   sigma_j = max(|y_j|, 1/ewt_j) * sqrt(uround)
 *   J[:,j] = (f(t, y + sigma_j*e_j) - f(t, y)) / sigma_j
 */
function cvLsDQJac(mem: CvodeMem, t: number, y: Float64Array, fy: Float64Array,
                   J: Float64Array[], tmp1: Float64Array, tmp2: Float64Array): number {
  const N = mem.cv_N;
  const lmem = mem.cv_lmem!;
  const uround = 2.2204460492503131e-16;
  const srur = Math.sqrt(uround);
  const fnorm_val = wrmsNorm(N, fy, mem.cv_ewt);
  const minInc = (fnorm_val !== 0.0) ? (1000.0 * Math.abs(mem.cv_h) * uround * N * fnorm_val) : 1.0;

  // Use tmp1 as perturbed y, tmp2 as perturbed f
  for (let j = 0; j < N; j++) {
    const yjsaved = y[j];
    let inc = Math.max(srur * Math.abs(yjsaved), minInc / mem.cv_ewt[j]);
    // Ensure inc has same sign (positive since we want to increase)

    // Perturb y[j]
    tmp1.set(y);
    tmp1[j] = yjsaved + inc;
    inc = tmp1[j] - yjsaved; // exact increment after rounding

    // Evaluate f(t, y + inc*e_j)
    const retval = mem.cv_f!(t, tmp1, tmp2, mem.cv_user_data);
    lmem.nfeDQ++;
    mem.cv_nfe++;
    if (retval !== 0) return retval;

    // Compute column j of Jacobian: (f_perturbed - f) / inc
    const inc_inv = 1.0 / inc;
    for (let i = 0; i < N; i++) {
      J[i][j] = (tmp2[i] - fy[i]) * inc_inv;
    }
  }

  return CV_SUCCESS;
}

/**
 * cvLsFree - Free linear solver memory (no-op in TS, GC handles it).
 */
function cvLsFree(mem: CvodeMem): void {
  mem.cv_lmem = null;
}
