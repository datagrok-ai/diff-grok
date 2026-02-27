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

import {LsodaContext, sm1} from './common';

/**
 * Scale step size. Rescales the yh array when h changes.
 */
export function scaleh(ctx: LsodaContext, rh: number): void {
  const c = ctx.common!;
  const neq = ctx.neq;
  const hmxi = ctx.opt!.hmxi;

  rh = Math.min(rh, c.rmax);
  rh = rh / Math.max(1., Math.abs(c.h) * hmxi * rh);

  // If meth = 1, restrict by stability region
  if (c.meth === 1) {
    c.irflag = 0;
    const pdh = Math.max(Math.abs(c.h) * c.pdlast, 0.000001);
    if ((rh * pdh * 1.00001) >= sm1[c.nq]) {
      rh = sm1[c.nq] / pdh;
      c.irflag = 1;
    }
  }

  let r = 1.;
  for (let j = 2; j <= c.nq + 1; j++) {
    r *= rh;
    for (let i = 1; i <= neq; i++)
      c.yh[j][i] *= r;
  }
  c.h *= rh;
  c.rc *= rh;
  c.ialth = c.nq + 1;
}
