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

import {LsodaContext, ETA, MAXCOR} from './common';
import {vmnorm} from './blas';
import {prja} from './prja';
import {solsy} from './solsy';
import {corfailure} from './corfailure';

/** Mutable state passed between correction and stoda */
export interface CorrectionState {
  del: number;
  delp: number;
  m: number;
}

/**
 * Corrector iteration.
 * Returns corflag: 0 = converged, 1 = reduce step & redo prediction, 2 = failure.
 */
export function correction(
  ctx: LsodaContext, y: Float64Array, pnorm: number,
  cs: CorrectionState, told: number,
): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  cs.m = 0;
  let rate = 0.;
  cs.del = 0.;

  for (let i = 1; i <= neq; i++)
    y[i] = c.yh[1][i];

  ctx.func(c.tn, y.subarray(1), c.savf.subarray(1), ctx.data);
  c.nfe++;

  while (true) {
    if (cs.m === 0) {
      if (c.ipup > 0) {
        const ierpj = prja(ctx, y);
        c.jcur = 1;
        c.ipup = 0;
        c.rc = 1.;
        c.nslp = c.nst;
        c.crate = 0.7;
        if (!ierpj)
          return corfailure(ctx, told);
      }
      for (let i = 1; i <= neq; i++)
        c.acor[i] = 0.;
    }

    if (c.miter === 0) {
      // Functional iteration
      for (let i = 1; i <= neq; i++) {
        c.savf[i] = c.h * c.savf[i] - c.yh[2][i];
        y[i] = c.savf[i] - c.acor[i];
      }
      cs.del = vmnorm(neq, y, c.ewt);
      for (let i = 1; i <= neq; i++) {
        y[i] = c.yh[1][i] + c.el[1] * c.savf[i];
        c.acor[i] = c.savf[i];
      }
    } else {
      // Chord method
      for (let i = 1; i <= neq; i++)
        y[i] = c.h * c.savf[i] - (c.yh[2][i] + c.acor[i]);
      solsy(ctx, y);
      cs.del = vmnorm(neq, y, c.ewt);
      for (let i = 1; i <= neq; i++) {
        c.acor[i] += y[i];
        y[i] = c.yh[1][i] + c.el[1] * c.acor[i];
      }
    }

    // Test for convergence
    if (cs.del <= 100. * pnorm * ETA)
      break;

    if (cs.m !== 0 || c.meth !== 1) {
      if (cs.m !== 0) {
        let rm = 1024.0;
        if (cs.del <= 1024. * cs.delp)
          rm = cs.del / cs.delp;
        rate = Math.max(rate, rm);
        c.crate = Math.max(0.2 * c.crate, rm);
      }
      const conit = 0.5 / (c.nq + 2);
      const dcon = cs.del * Math.min(1., 1.5 * c.crate) / (c.tesco[c.nq][2] * conit);
      if (dcon <= 1.) {
        c.pdest = Math.max(c.pdest, rate / Math.abs(c.h * c.el[1]));
        if (c.pdest !== 0.)
          c.pdlast = c.pdest;
        break;
      }
    }

    // Corrector iteration failed to converge
    cs.m++;
    if (cs.m === MAXCOR || (cs.m >= 2 && cs.del > 2. * cs.delp)) {
      if (c.miter === 0 || c.jcur === 1)
        return corfailure(ctx, told);

      c.ipup = c.miter;
      // Restart corrector if Jacobian is recomputed
      cs.m = 0;
      rate = 0.;
      cs.del = 0.;
      for (let i = 1; i <= neq; i++)
        y[i] = c.yh[1][i];
      ctx.func(c.tn, y.subarray(1), c.savf.subarray(1), ctx.data);
      c.nfe++;
    } else {
      cs.delp = cs.del;
      ctx.func(c.tn, y.subarray(1), c.savf.subarray(1), ctx.data);
      c.nfe++;
    }
  }

  return 0;
}
