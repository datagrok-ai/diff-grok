import { LsodaContext, sm1 } from './common';
import { vmnorm } from './blas';

/** Mutable rh output */
export interface OrderSwitchResult {
  rh: number;
}

/**
 * Decide whether to change the order and/or step size.
 * Returns orderflag: 0 = no change, 1 = change h only, 2 = change both h and nq.
 */
export function orderswitch(
  ctx: LsodaContext, rhup: number, dsm: number,
  out: OrderSwitchResult, kflag: number, maxord: number,
): number {
  const c = ctx.common!;
  const neq = ctx.neq;

  const exsm = 1. / (c.nq + 1);
  let rhsm = 1. / (1.2 * Math.pow(dsm, exsm) + 0.0000012);

  let rhdn = 0.;
  if (c.nq !== 1) {
    const ddn = vmnorm(neq, c.yh[c.nq + 1], c.ewt) / c.tesco[c.nq][1];
    const exdn = 1. / c.nq;
    rhdn = 1. / (1.3 * Math.pow(ddn, exdn) + 0.0000013);
  }

  // If meth = 1, limit rh according to stability region
  if (c.meth === 1) {
    const pdh = Math.max(Math.abs(c.h) * c.pdlast, 0.000001);
    if (c.nq + 1 < maxord + 1)
      rhup = Math.min(rhup, sm1[c.nq + 1] / pdh);
    rhsm = Math.min(rhsm, sm1[c.nq] / pdh);
    if (c.nq > 1)
      rhdn = Math.min(rhdn, sm1[c.nq - 1] / pdh);
    c.pdest = 0.;
  }

  let newq: number;
  if (rhsm >= rhup) {
    if (rhsm >= rhdn) {
      newq = c.nq;
      out.rh = rhsm;
    } else {
      newq = c.nq - 1;
      out.rh = rhdn;
      if (kflag < 0 && out.rh > 1.) out.rh = 1.;
    }
  } else {
    if (rhup <= rhdn) {
      newq = c.nq - 1;
      out.rh = rhdn;
      if (kflag < 0 && out.rh > 1.) out.rh = 1.;
    } else {
      out.rh = rhup;
      if (out.rh >= 1.1) {
        const r = c.el[c.nq + 1] / (c.nq + 1);
        c.nq = c.nq + 1;
        for (let i = 1; i <= neq; i++)
          c.yh[c.nq + 1][i] = c.acor[i] * r;
        return 2;
      } else {
        c.ialth = 3;
        return 0;
      }
    }
  }

  // If meth = 1 and h is restricted by stability, bypass 10% test
  if (c.meth === 1) {
    const pdh = Math.max(Math.abs(c.h) * c.pdlast, 0.000001);
    if (out.rh * pdh * 1.00001 < sm1[newq])
      if (kflag === 0 && out.rh < 1.1) {
        c.ialth = 3;
        return 0;
      }
  } else {
    if (kflag === 0 && out.rh < 1.1) {
      c.ialth = 3;
      return 0;
    }
  }

  if (kflag <= -2)
    out.rh = Math.min(out.rh, 0.2);

  if (newq === c.nq) {
    return 1;
  }
  c.nq = newq;
  return 2;
}
