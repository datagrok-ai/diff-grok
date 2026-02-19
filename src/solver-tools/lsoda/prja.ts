import { LsodaContext, ETA, SQRTETA } from './common';
import { vmnorm, fnorm, dgefa } from './blas';

/**
 * Compute and process the matrix P = I - h * el[1] * J,
 * where J is approximated by finite differences.
 * Returns 1 on success, 0 on failure (singular matrix).
 */
export function prja(ctx: LsodaContext, y: Float64Array): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  c.nje++;
  const hl0 = c.h * c.el[1];

  if (c.miter !== 2) {
    console.error('[prja] miter != 2');
    return 0;
  }

  let fac = vmnorm(neq, c.savf, c.ewt);
  let r0 = 1000. * Math.abs(c.h) * ETA * neq * fac;
  if (r0 === 0.) r0 = 1.;

  for (let j = 1; j <= neq; j++) {
    const yj = y[j];
    const r = Math.max(SQRTETA * Math.abs(yj), r0 / c.ewt[j]);
    y[j] += r;
    fac = -hl0 / r;
    ctx.func(c.tn, y.subarray(1), c.acor.subarray(1), ctx.data);
    for (let i = 1; i <= neq; i++)
      c.wm[i][j] = (c.acor[i] - c.savf[i]) * fac;
    y[j] = yj;
  }
  c.nfe += neq;

  // Compute norm of Jacobian
  c.pdnorm = fnorm(neq, c.wm, c.ewt) / Math.abs(hl0);

  // Add identity matrix
  for (let i = 1; i <= neq; i++)
    c.wm[i][i] += 1.;

  // LU decomposition on P
  const ier = dgefa(c.wm, neq, c.ipvt);
  if (ier !== 0) return 0;

  return 1;
}
