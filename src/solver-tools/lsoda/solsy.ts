/**
 * TypeScript port of liblsoda (https://github.com/sdwfrost/liblsoda)
 * Original LSODA algorithm by Linda Petzold & Alan Hindmarsh (1983)
 *
 * Copyright (c) 2011 Yu Feng, McWilliam Cosmology Center, Carnegie Mellon University
 * Copyright (c) 2016 Simon Frost
 * Copyright (c) 2026 Datagrok
 *
 * Licensed under the MIT License. See THIRD_PARTY_LICENSES for the full
 * license text and provenance chain of the original code.
 */

import {LsodaContext} from './common';
import {dgesl} from './blas';

/**
 * Solve the linear system arising from a chord iteration.
 * Only miter=2 (dense) is supported.
 */
export function solsy(ctx: LsodaContext, y: Float64Array): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  if (c.miter !== 2)
    throw new Error('[solsy] miter != 2 not implemented');


  dgesl(c.wm, neq, c.ipvt, y, 0);
  return 1;
}
