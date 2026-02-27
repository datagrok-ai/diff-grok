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

import {LsodaContext, ETA} from './common';

/**
 * Interpolate solution at time t using Nordsieck history array.
 * k: derivative order (0 for solution)
 * dky: output array (1-indexed)
 * Returns 0 on success, negative on error.
 */
export function intdy(ctx: LsodaContext, t: number, k: number, dky: Float64Array): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  if (k < 0 || k > c.nq) {
    console.error(`[intdy] k = ${k} illegal`);
    return -1;
  }

  const tp = c.tn - c.hu - 100. * ETA * (c.tn + c.hu);
  if ((t - tp) * (t - c.tn) > 0.) {
    console.error(`intdy -- t = ${t} illegal. t not in interval tcur - hu to tcur`);
    return -2;
  }

  const s = (t - c.tn) / c.h;
  let ic = 1;
  for (let jj = (c.nq + 1) - k; jj <= c.nq; jj++)
    ic *= jj;
  let cc = ic;

  for (let i = 1; i <= neq; i++)
    dky[i] = cc * c.yh[c.nq + 1][i];

  for (let j = c.nq - 1; j >= k; j--) {
    const jp1 = j + 1;
    ic = 1;
    for (let jj = jp1 - k; jj <= j; jj++)
      ic *= jj;
    cc = ic;
    for (let i = 1; i <= neq; i++)
      dky[i] = cc * c.yh[jp1][i] + s * dky[i];
  }

  if (k === 0) return 0;

  const r = Math.pow(c.h, -k);
  for (let i = 1; i <= neq; i++)
    dky[i] *= r;

  return 0;
}
