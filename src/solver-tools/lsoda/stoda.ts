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

import {LsodaContext, CCMAX, MSBP} from './common';
import {vmnorm} from './blas';
import {cfode} from './cfode';
import {scaleh} from './scaleh';
import {correction, CorrectionState} from './correction';
import {orderswitch, OrderSwitchResult} from './orderswitch';
import {methodswitch, MethodSwitchResult} from './methodswitch';

/**
 * Perform one step of the integration.
 *
 * jstart:  0 = first step, >0 = continue, -1 = new h/meth/miter, -2 = new h only
 * Returns kflag: 0 = success, -1 = error test failed, -2 = convergence failed.
 */
export function stoda(ctx: LsodaContext, y: Float64Array, jstart: number): number {
  const c = ctx.common!;
  const hmin = ctx.opt!.hmin;
  const mxords = ctx.opt!.mxords;
  const mxordn = ctx.opt!.mxordn;
  const neq = ctx.neq;

  let kflag = 0;
  const told = c.tn;
  c.ncf = 0;
  let delp = 0.;

  let maxord = mxordn;
  if (c.meth === 2) maxord = mxords;

  // Helper: endstoda — finalize step, scale acor
  function endstoda(): void {
    const r = 1. / c.tesco[c.nqu][2];
    for (let i = 1; i <= neq; i++)
      c.acor[i] *= r;
    c.hold = c.h;
  }

  // Helper: resetcoeff — reload el from elco for current order
  function resetcoeff(): void {
    const el0 = c.el[1];
    for (let i = 1; i <= c.nq + 1; i++)
      c.el[i] = c.elco[c.nq][i];
    c.rc = c.rc * c.el[1] / el0;
  }

  // --- Initialization based on jstart ---

  if (jstart === 0) {
    c.nq = 1;
    c.ialth = 2;
    c.rmax = 10000.;
    c.rc = 0.;
    c.crate = 0.7;
    c.hold = c.h;
    c.nslp = 0;
    c.ipup = c.miter;
    c.el[1] = 1.0;
    c.icount = 20;
    c.irflag = 0;
    c.pdest = 0.;
    c.pdlast = 0.;
    cfode(ctx, 1);
    resetcoeff();
  }

  if (jstart === -1) {
    c.ipup = c.miter;
    if (c.ialth === 1) c.ialth = 2;
    if (c.meth !== c.mused) {
      cfode(ctx, c.meth);
      c.ialth = c.nq + 1;
      resetcoeff();
    }
    if (c.h !== c.hold) {
      const rh = c.h / c.hold;
      c.h = c.hold;
      scaleh(ctx, rh);
    }
  }

  if (jstart === -2) {
    if (c.h !== c.hold) {
      const rh = c.h / c.hold;
      c.h = c.hold;
      scaleh(ctx, rh);
    }
  }

  // --- Main integration loop ---

  let dsm = 0.0;
  const cs: CorrectionState = {del: 0, delp: delp, m: 0};
  const osResult: OrderSwitchResult = {rh: 0};
  const msResult: MethodSwitchResult = {rh: 0};
  let rh: number;

  outer:
  while (true) {
    // Before the corrector starts
    c.jcur = 0;

    inner:
    while (true) {
      if (Math.abs(c.rc - 1.) > CCMAX)
        c.ipup = c.miter;
      if (c.nst >= c.nslp + MSBP)
        c.ipup = c.miter;

      c.tn += c.h;

      // Prediction: multiply yh by Pascal triangle matrix
      for (let j = c.nq; j >= 1; j--) {
        for (let i1 = j; i1 <= c.nq; i1++) {
          for (let i = 1; i <= neq; i++)
            c.yh[i1][i] += c.yh[i1 + 1][i];
        }
      }

      const pnorm = vmnorm(neq, c.yh[1], c.ewt);

      cs.del = 0;
      cs.delp = delp;
      cs.m = 0;
      const corflag = correction(ctx, y, pnorm, cs, told);
      delp = cs.delp;

      if (corflag === 0) break inner;

      if (corflag === 1) {
        rh = Math.max(0.25, hmin / Math.abs(c.h));
        scaleh(ctx, rh);
        continue inner;
      }

      if (corflag === 2) {
        kflag = -2;
        c.hold = c.h;
        jstart = 1;
        return kflag;
      }
    }

    // Local error test
    if (cs.m === 0)
      dsm = cs.del / c.tesco[c.nq][2];
    if (cs.m > 0)
      dsm = vmnorm(neq, c.acor, c.ewt) / c.tesco[c.nq][2];

    if (dsm <= 1.) {
      // Successful step
      kflag = 0;
      c.nst++;
      c.hu = c.h;
      c.nqu = c.nq;
      c.mused = c.meth;

      for (let j = 1; j <= c.nq + 1; j++) {
        const r = c.el[j];
        for (let i = 1; i <= neq; i++)
          c.yh[j][i] += r * c.acor[i];
      }

      c.icount--;
      if (c.icount < 0) {
        msResult.rh = 0;
        methodswitch(ctx, dsm, vmnorm(neq, c.yh[1], c.ewt), msResult);
        if (c.meth !== c.mused) {
          rh = Math.max(msResult.rh, hmin / Math.abs(c.h));
          scaleh(ctx, rh);
          c.rmax = 10.;
          endstoda();
          break outer;
        }
      }

      // No method switch — do step/order selection
      c.ialth--;
      if (c.ialth === 0) {
        let rhup = 0.;
        if (c.nq + 1 !== maxord + 1) {
          for (let i = 1; i <= neq; i++)
            c.savf[i] = c.acor[i] - c.yh[maxord + 1][i];
          const dup = vmnorm(neq, c.savf, c.ewt) / c.tesco[c.nq][3];
          const exup = 1. / (c.nq + 2);
          rhup = 1. / (1.4 * Math.pow(dup, exup) + 0.0000014);
        }

        osResult.rh = 0;
        const orderflag = orderswitch(ctx, rhup, dsm, osResult, kflag, maxord);

        if (orderflag === 0) {
          endstoda();
          break outer;
        }
        if (orderflag === 1) {
          rh = Math.max(osResult.rh, hmin / Math.abs(c.h));
          scaleh(ctx, rh);
          c.rmax = 10.;
          endstoda();
          break outer;
        }
        if (orderflag === 2) {
          resetcoeff();
          rh = Math.max(osResult.rh, hmin / Math.abs(c.h));
          scaleh(ctx, rh);
          c.rmax = 10.;
          endstoda();
          break outer;
        }
      }

      if (c.ialth > 1 || c.nq + 1 === maxord + 1) {
        endstoda();
        break outer;
      }

      for (let i = 1; i <= neq; i++)
        c.yh[maxord + 1][i] = c.acor[i];
      endstoda();
      break outer;
    } else {
      // Error test failed
      kflag--;
      c.tn = told;

      // Retract yh array
      for (let j = c.nq; j >= 1; j--) {
        for (let i1 = j; i1 <= c.nq; i1++) {
          for (let i = 1; i <= neq; i++)
            c.yh[i1][i] -= c.yh[i1 + 1][i];
        }
      }

      c.rmax = 2.;
      if (Math.abs(c.h) <= hmin * 1.00001) {
        kflag = -1;
        c.hold = c.h;
        jstart = 1;
        break outer;
      }

      if (kflag > -3) {
        osResult.rh = 0;
        const orderflag = orderswitch(ctx, 0., dsm, osResult, kflag, maxord);
        if (orderflag === 1 || orderflag === 0) {
          if (orderflag === 0)
            osResult.rh = Math.min(osResult.rh, 0.2);
          osResult.rh = Math.max(osResult.rh, hmin / Math.abs(c.h));
          scaleh(ctx, osResult.rh);
        }
        if (orderflag === 2) {
          resetcoeff();
          osResult.rh = Math.max(osResult.rh, hmin / Math.abs(c.h));
          scaleh(ctx, osResult.rh);
        }
        continue outer;
      } else {
        // 3 or more failures
        if (kflag === -10) {
          kflag = -1;
          c.hold = c.h;
          jstart = 1;
          break outer;
        } else {
          rh = 0.1;
          rh = Math.max(hmin / Math.abs(c.h), rh);
          c.h *= rh;
          for (let i = 1; i <= neq; i++)
            y[i] = c.yh[1][i];
          ctx.func(c.tn, y.subarray(1), c.savf.subarray(1), ctx.data);
          c.nfe++;
          for (let i = 1; i <= neq; i++)
            c.yh[2][i] = c.h * c.savf[i];
          c.ipup = c.miter;
          c.ialth = 5;
          if (c.nq === 1) continue outer;
          c.nq = 1;
          resetcoeff();
          continue outer;
        }
      }
    }
  }

  return kflag;
}
