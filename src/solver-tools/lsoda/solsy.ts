import { LsodaContext } from './common';
import { dgesl } from './blas';

/**
 * Solve the linear system arising from a chord iteration.
 * Only miter=2 (dense) is supported.
 */
export function solsy(ctx: LsodaContext, y: Float64Array): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  if (c.miter !== 2) {
    throw new Error('[solsy] miter != 2 not implemented');
  }

  dgesl(c.wm, neq, c.ipvt, y, 0);
  return 1;
}
