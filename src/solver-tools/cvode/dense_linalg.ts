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

// Dense linear algebra routines for CVODE (0-based indexing).
//
// These are fresh 0-based implementations of LU factorization with partial
// pivoting and triangular solve. The existing src/blas.ts uses 1-based
// indexing inherited from the Fortran/Linpack convention; CVODE needs 0-based.
//
// Matrix layout: row-major Float64Array[], where a[i][j] is row i, column j.
// This matches the LSODA convention (see CVODE-C-to-TS.md section 2.3) and
// differs from SUNDIALS column-major storage. The algorithms are translated
// from SUNDIALS sundials_dense.c (denseGETRF / denseGETRS) with the storage
// order transposed.

/**
 * LU factorization of an n-by-n dense matrix by Gaussian elimination
 * with partial pivoting (row-major, 0-based).
 *
 * On exit, `a` is overwritten with the LU factors:
 * - The upper triangle (including diagonal) holds U.
 * - The strict lower triangle holds the multipliers L (unit diagonal implied).
 *
 * The pivot vector `ipvt` records row interchanges: at step k, row k was
 * swapped with row ipvt[k] before elimination.
 *
 * @param a   - n-by-n matrix stored as `Float64Array[]` (a[i] = row i).
 *              Overwritten with LU factors in place.
 * @param n   - matrix dimension (number of rows and columns).
 * @param ipvt - pivot index output array (`Int32Array` of length n, 0-based).
 * @returns 0 on success, or k+1 (1-based) if U[k][k] is exactly zero
 *          (the matrix is singular or nearly so). A nonzero return means
 *          the factorization completed but the factor U is singular.
 */
export function dgefa(a: Float64Array[], n: number, ipvt: Int32Array): number {
  let info = 0;

  for (let k = 0; k < n; k++) {
    // --- Find pivot: row index l in [k, n) with largest |a[l][k]| ---
    let l = k;
    let maxAbs = Math.abs(a[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(a[i][k]);
      if (v > maxAbs) {
        maxAbs = v;
        l = i;
      }
    }
    ipvt[k] = l;

    // --- Check for zero pivot ---
    if (a[l][k] === 0.0) {
      info = k + 1;
      continue;
    }

    // --- Swap rows k and l (swap row references, O(1)) ---
    if (l !== k) {
      const tmp = a[k];
      a[k] = a[l];
      a[l] = tmp;
    }

    // --- Compute multipliers: a[i][k] = a[i][k] / a[k][k], i = k+1..n-1 ---
    // These multipliers form the strict lower triangle of L.
    const pivot = a[k][k];
    const mult = 1.0 / pivot;
    for (let i = k + 1; i < n; i++) {
      a[i][k] *= mult;
    }

    // --- Eliminate: update the trailing (n-k-1)-by-(n-k-1) submatrix ---
    // a[i][j] -= a[i][k] * a[k][j],  for i = k+1..n-1, j = k+1..n-1
    for (let i = k + 1; i < n; i++) {
      const m = a[i][k]; // multiplier (already computed above)
      if (m !== 0.0) {
        const rowI = a[i];
        const rowK = a[k];
        for (let j = k + 1; j < n; j++) {
          rowI[j] -= m * rowK[j];
        }
      }
    }
  }

  return info;
}

/**
 * Solve A * x = b using the LU factorization produced by {@link dgefa}.
 *
 * On entry, `b` contains the right-hand side vector. On exit, `b` is
 * overwritten with the solution vector x.
 *
 * The solve proceeds in three stages:
 * 1. Permute b according to the pivot vector (row interchanges).
 * 2. Forward substitution: solve L * y = P * b.
 * 3. Back substitution: solve U * x = y.
 *
 * @param a    - LU-factored n-by-n matrix from dgefa (row-major, 0-based).
 * @param n    - matrix dimension.
 * @param ipvt - pivot indices from dgefa (0-based, length n).
 * @param b    - right-hand side on input (length n, 0-based), solution on output.
 */
export function dgesl(
  a: Float64Array[], n: number, ipvt: Int32Array, b: Float64Array,
): void {
  // --- Permute b according to pivot vector ---
  for (let k = 0; k < n; k++) {
    const p = ipvt[k];
    if (p !== k) {
      const tmp = b[k];
      b[k] = b[p];
      b[p] = tmp;
    }
  }

  // --- Forward substitution: solve L * y = P*b ---
  // L has unit diagonal; multipliers are in strict lower triangle of a.
  // y[i] = b[i] - sum_{k=0}^{i-1} a[i][k] * y[k],  but since we process
  // column-by-column (outer loop on k), we update all rows below k:
  for (let k = 0; k < n - 1; k++) {
    for (let i = k + 1; i < n; i++) {
      b[i] -= a[i][k] * b[k];
    }
  }

  // --- Back substitution: solve U * x = y ---
  // U is in the upper triangle (including diagonal) of a.
  for (let k = n - 1; k >= 0; k--) {
    b[k] /= a[k][k];
    for (let i = 0; i < k; i++) {
      b[i] -= a[i][k] * b[k];
    }
  }
}
