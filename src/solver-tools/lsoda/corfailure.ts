import {LsodaContext, MXNCF} from './common';

/**
 * Handle corrector convergence failure.
 * Retracts yh array and increments ncf.
 * Returns 2 if unrecoverable, 1 if step can be retried.
 */
export function corfailure(ctx: LsodaContext, told: number): number {
  const c = ctx.common!;
  const neq = ctx.neq;
  const hmin = ctx.opt!.hmin;

  c.ncf++;
  c.rmax = 2.;
  c.tn = told;

  for (let j = c.nq; j >= 1; j--) {
    for (let i1 = j; i1 <= c.nq; i1++) {
      for (let i = 1; i <= neq; i++)
        c.yh[i1][i] -= c.yh[i1 + 1][i];
    }
  }

  if (Math.abs(c.h) <= hmin * 1.00001 || c.ncf === MXNCF)
    return 2;


  c.ipup = c.miter;
  return 1;
}
