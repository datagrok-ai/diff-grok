// BLAS linear algebra routines
// All arrays are 1-indexed (element [0] unused), matching the C/Fortran convention.

/**
 * Compute dy = da * dx + dy
 * Arrays are 1-indexed with stride incx/incy.
 */
export function daxpy(
  n: number, da: number, dx: Float64Array, incx: number,
  dy: Float64Array, incy: number,
): void {
  if (n < 0 || da === 0.) return;

  if (incx !== incy || incx < 1) {
    let ix = 1;
    let iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (let i = 1; i <= n; i++) {
      dy[iy] = dy[iy] + da * dx[ix];
      ix += incx;
      iy += incy;
    }
    return;
  }

  if (incx === 1) {
    const m = n % 4;
    if (m !== 0) {
      for (let i = 1; i <= m; i++)
        dy[i] = dy[i] + da * dx[i];
      if (n < 4) return;
    }
    for (let i = m + 1; i <= n; i += 4) {
      dy[i] = dy[i] + da * dx[i];
      dy[i + 1] = dy[i + 1] + da * dx[i + 1];
      dy[i + 2] = dy[i + 2] + da * dx[i + 2];
      dy[i + 3] = dy[i + 3] + da * dx[i + 3];
    }
    return;
  }

  for (let i = 1; i <= n * incx; i += incx)
    dy[i] = da * dx[i] + dy[i];
}

/**
 * Compute dot product dx . dy
 * Arrays are 1-indexed with stride incx/incy.
 */
export function ddot(
  n: number, dx: Float64Array, incx: number,
  dy: Float64Array, incy: number,
): number {
  let dotprod = 0.;
  if (n <= 0) return dotprod;

  if (incx !== incy || incx < 1) {
    let ix = 1;
    let iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (let i = 1; i <= n; i++) {
      dotprod += dx[ix] * dy[iy];
      ix += incx;
      iy += incy;
    }
    return dotprod;
  }

  if (incx === 1) {
    for (let i = 1; i <= n; i++)
      dotprod += dx[i] * dy[i];
    return dotprod;
  }

  for (let i = 1; i <= n * incx; i += incx)
    dotprod += dx[i] * dy[i];
  return dotprod;
}

/**
 * Scale vector: dx = da * dx
 * Array is 1-indexed with stride incx.
 */
export function dscal(n: number, da: number, dx: Float64Array, incx: number): void {
  if (n <= 0) return;

  if (incx !== 1) {
    for (let i = 1; i <= n * incx; i += incx)
      dx[i] = da * dx[i];
    return;
  }

  const m = n % 5;
  if (m !== 0) {
    for (let i = 1; i <= m; i++)
      dx[i] = da * dx[i];
    if (n < 5) return;
  }
  for (let i = m + 1; i <= n; i += 5) {
    dx[i] = da * dx[i];
    dx[i + 1] = da * dx[i + 1];
    dx[i + 2] = da * dx[i + 2];
    dx[i + 3] = da * dx[i + 3];
    dx[i + 4] = da * dx[i + 4];
  }
}

/**
 * Find index of element with maximum absolute value.
 * Array is 1-indexed. Returns 1-based index (0 if n <= 0).
 */
export function idamax(n: number, dx: Float64Array, incx: number): number {
  let xindex = 0;
  if (n <= 0) return xindex;
  xindex = 1;
  if (n <= 1 || incx <= 0) return xindex;

  if (incx !== 1) {
    let dmax = Math.abs(dx[1]);
    let ii = 2;
    for (let i = 1 + incx; i <= n * incx; i += incx) {
      const xmag = Math.abs(dx[i]);
      if (xmag > dmax) {
        xindex = ii;
        dmax = xmag;
      }
      ii++;
    }
    return xindex;
  }

  let dmax = Math.abs(dx[1]);
  for (let i = 2; i <= n; i++) {
    const xmag = Math.abs(dx[i]);
    if (xmag > dmax) {
      xindex = i;
      dmax = xmag;
    }
  }
  return xindex;
}

/**
 * LU factorization by Gaussian elimination with partial pivoting.
 * a: 2D matrix (a[i][j], 1-indexed rows and columns)
 * n: matrix dimension
 * ipvt: pivot indices output (1-indexed)
 * Returns info: 0 = normal, k = U[k][k] == 0 (singular).
 */
export function dgefa(a: Float64Array[], n: number, ipvt: Int32Array): number {
  let info = 0;

  for (let k = 1; k <= n - 1; k++) {
    // Find pivot index. a[k] + k - 1 in C => subarray from index k
    const j = idamax(n - k + 1, a[k].subarray(k - 1), 1) + k - 1;
    ipvt[k] = j;

    if (a[k][j] === 0.) {
      info = k;
      continue;
    }

    // Interchange if necessary
    if (j !== k) {
      const t = a[k][j];
      a[k][j] = a[k][k];
      a[k][k] = t;
    }

    // Compute multipliers
    const t = -1. / a[k][k];
    dscal(n - k, t, a[k].subarray(k), 1);

    // Column elimination with row indexing
    for (let i = k + 1; i <= n; i++) {
      const t2 = a[i][j];
      if (j !== k) {
        a[i][j] = a[i][k];
        a[i][k] = t2;
      }
      daxpy(n - k, t2, a[k].subarray(k), 1, a[i].subarray(k), 1);
    }
  }

  ipvt[n] = n;
  if (a[n][n] === 0.) info = n;
  return info;
}

/**
 * Solve linear system using LU factors from dgefa.
 * a: LU factored matrix (1-indexed)
 * n: dimension
 * ipvt: pivot indices (1-indexed)
 * b: right-hand side on input, solution on output (1-indexed)
 * job: 0 = solve A*x=b, nonzero = solve A^T*x=b
 */
export function dgesl(
  a: Float64Array[], n: number, ipvt: Int32Array,
  b: Float64Array, job: number,
): void {
  if (job === 0) {
    // Solve L * y = b
    for (let k = 1; k <= n; k++) {
      const t = ddot(k - 1, a[k], 1, b, 1);
      b[k] = (b[k] - t) / a[k][k];
    }
    // Solve U * x = y
    for (let k = n - 1; k >= 1; k--) {
      b[k] = b[k] + ddot(n - k, a[k].subarray(k), 1, b.subarray(k), 1);
      const j = ipvt[k];
      if (j !== k) {
        const t = b[j];
        b[j] = b[k];
        b[k] = t;
      }
    }
    return;
  }

  // Solve Transpose(U) * y = b
  for (let k = 1; k <= n - 1; k++) {
    const j = ipvt[k];
    const t = b[j];
    if (j !== k) {
      b[j] = b[k];
      b[k] = t;
    }
    daxpy(n - k, t, a[k].subarray(k), 1, b.subarray(k), 1);
  }
  // Solve Transpose(L) * x = y
  for (let k = n; k >= 1; k--) {
    b[k] = b[k] / a[k][k];
    const t = -b[k];
    daxpy(k - 1, t, a[k], 1, b, 1);
  }
}

/**
 * Weighted max-norm: max(i=1..n) |v[i]| * w[i]
 * Arrays are 1-indexed.
 */
export function vmnorm(n: number, v: Float64Array, w: Float64Array): number {
  let vm = 0.;
  for (let i = 1; i <= n; i++)
    vm = Math.max(vm, Math.abs(v[i]) * w[i]);
  return vm;
}

/**
 * Weighted matrix norm consistent with vmnorm.
 * fnorm = max(i=1..n) ( w[i] * sum(j=1..n) |a[i][j]| / w[j] )
 * a is 2D (1-indexed), w is 1-indexed.
 */
export function fnorm(n: number, a: Float64Array[], w: Float64Array): number {
  let an = 0.;
  for (let i = 1; i <= n; i++) {
    let sum = 0.;
    const ap1 = a[i];
    for (let j = 1; j <= n; j++)
      sum += Math.abs(ap1[j]) / w[j];
    an = Math.max(an, sum * w[i]);
  }
  return an;
}
