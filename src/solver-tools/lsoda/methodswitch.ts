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

import {LsodaContext, ETA, RATIO, sm1, cm1, cm2} from './common';
import {vmnorm} from './blas';

/** Mutable rh output */
export interface MethodSwitchResult {
  rh: number;
}

/**
 * Consider switching between Adams (meth=1) and BDF (meth=2).
 * Modifies ctx.common.meth, miter, nq, icount, pdlast on switch.
 */
export function methodswitch(
  ctx: LsodaContext, dsm: number, pnorm: number, out: MethodSwitchResult,
): void {
  const c = ctx.common!;
  const neq = ctx.neq;
  const mxordn = ctx.opt!.mxordn;
  const mxords = ctx.opt!.mxords;

  if (c.meth === 1) {
    // Currently Adams, consider switching to BDF
    if (c.nq > 5) return;

    let rh2: number;
    let nqm2: number;

    if (dsm <= (100. * pnorm * ETA) || c.pdest === 0.) {
      if (c.irflag === 0) return;
      rh2 = 2.;
      nqm2 = Math.min(c.nq, mxords);
    } else {
      const exsm = 1. / (c.nq + 1);
      let rh1 = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);
      let rh1it = 2. * rh1;
      const pdh = c.pdlast * Math.abs(c.h);
      if ((pdh * rh1) > 0.00001)
        rh1it = sm1[c.nq] / pdh;
      rh1 = Math.min(rh1, rh1it);

      if (c.nq > mxords) {
        nqm2 = mxords;
        const lm2 = mxords + 1;
        const exm2 = 1. / lm2;
        const lm2p1 = lm2 + 1;
        const dm2 = vmnorm(neq, c.yh[lm2p1], c.ewt) / cm2[mxords];
        rh2 = 1. / (1.2 * Math.pow(dm2, exm2) + 0.0000012);
      } else {
        const dm2 = dsm * (cm1[c.nq] / cm2[c.nq]);
        rh2 = 1. / (1.2 * Math.pow(dm2, exsm) + 0.0000012);
        nqm2 = c.nq;
      }
      if (rh2 < RATIO * rh1) return;
    }

    // Switch to BDF
    out.rh = rh2;
    c.icount = 20;
    c.meth = 2;
    c.miter = 2;
    c.pdlast = 0.;
    c.nq = nqm2;
    return;
  }

  // Currently BDF (meth=2), consider switching to Adams
  const exsm = 1. / (c.nq + 1);
  let nqm1: number;
  let exm1: number;
  let rh1: number;
  let dm1: number;

  if (mxordn < c.nq) {
    nqm1 = mxordn;
    const lm1 = mxordn + 1;
    exm1 = 1. / lm1;
    const lm1p1 = lm1 + 1;
    dm1 = vmnorm(neq, c.yh[lm1p1], c.ewt) / cm1[mxordn];
    rh1 = 1. / (1.2 * Math.pow(dm1, exm1) + 0.0000012);
  } else {
    dm1 = dsm * (cm2[c.nq] / cm1[c.nq]);
    rh1 = 1. / (1.2 * Math.pow(dm1, exsm) + 0.0000012);
    nqm1 = c.nq;
    exm1 = exsm;
  }

  let rh1it = 2. * rh1;
  const pdh = c.pdnorm * Math.abs(c.h);
  if ((pdh * rh1) > 0.00001)
    rh1it = sm1[nqm1] / pdh;
  rh1 = Math.min(rh1, rh1it);

  const rh2 = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);
  if ((rh1 * RATIO) < (5. * rh2)) return;

  const alpha = Math.max(0.001, rh1);
  dm1 *= Math.pow(alpha, exm1);
  if (dm1 <= 1000. * ETA * pnorm) return;

  // Switch to Adams
  out.rh = rh1;
  c.icount = 20;
  c.meth = 1;
  c.miter = 0;
  c.pdlast = 0.;
  c.nq = nqm1;
}
