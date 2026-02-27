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

import {
  LsodaContext, LsodaCommon, LsodaOpt, LsodaFunc, NordsieckSnapshot, ETA,
} from './common';
import {vmnorm} from './blas';
import {intdy} from './intdy';
import {stoda} from './stoda';
import {DenseOutput} from './dense';

// ---- Internal helpers ----

function captureSnapshot(c: LsodaCommon, neq: number): NordsieckSnapshot {
  const nq = c.nq;
  const yh: Float64Array[] = new Array(nq + 1);
  for (let j = 0; j <= nq; j++) {
    const row = new Float64Array(neq);
    for (let i = 0; i < neq; i++)
      row[i] = c.yh[j + 1][i + 1]; // convert 1-based to 0-based
    yh[j] = row;
  }
  return {tn: c.tn, h: c.h, hu: c.hu, nq, yh};
}

function ewset(ctx: LsodaContext, ycur: Float64Array): void {
  const c = ctx.common!;
  const neq = ctx.neq;
  const rtol = ctx.opt!.rtol;
  const atol = ctx.opt!.atol;
  for (let i = 1; i <= neq; i++)
    c.ewt[i] = rtol[i] * Math.abs(ycur[i]) + atol[i];
  for (let i = 1; i <= neq; i++)
    c.ewt[i] = 1. / c.ewt[i];
}

function checkOpt(ctx: LsodaContext, opt: LsodaOpt): boolean {
  const mxstp0 = 500;
  const mord = [0, 12, 5];

  if (ctx.state === 0) ctx.state = 1;
  if (ctx.state === 1) {
    opt.h0 = 0.;
    opt.mxordn = mord[1];
    opt.mxords = mord[2];
  }

  if (ctx.neq <= 0) {
    ctx.error = `[lsoda] neq = ${ctx.neq} is less than 1`;
    return false;
  }

  // Check rtol and atol (1-indexed)
  if (ctx.state === 1 || ctx.state === 3) {
    for (let i = 1; i <= ctx.neq; i++) {
      if (opt.rtol[i] < 0.) {
        ctx.error = `[lsoda] rtol = ${opt.rtol[i]} is less than 0.`;
        return false;
      }
      if (opt.atol[i] < 0.) {
        ctx.error = `[lsoda] atol = ${opt.atol[i]} is less than 0.`;
        return false;
      }
    }
  }

  if (opt.itask === 0) opt.itask = 1;
  if (opt.itask < 1 || opt.itask > 5) {
    ctx.error = `[lsoda] illegal itask = ${opt.itask}`;
    return false;
  }
  if (opt.ixpr < 0 || opt.ixpr > 1) {
    ctx.error = `[lsoda] ixpr = ${opt.ixpr} is illegal`;
    return false;
  }
  if (opt.mxstep < 0) {
    ctx.error = '[lsoda] mxstep < 0';
    return false;
  }
  if (opt.mxstep === 0) opt.mxstep = mxstp0;
  if (opt.mxhnil < 0) {
    ctx.error = '[lsoda] mxhnil < 0';
    return false;
  }

  if (ctx.state === 1) {
    if (opt.mxordn < 0) {
      ctx.error = `[lsoda] mxordn = ${opt.mxordn} is less than 0`;
      return false;
    }
    if (opt.mxordn === 0) opt.mxordn = 100;
    opt.mxordn = Math.min(opt.mxordn, mord[1]);
    if (opt.mxords < 0) {
      ctx.error = `[lsoda] mxords = ${opt.mxords} is less than 0`;
      return false;
    }
    if (opt.mxords === 0) opt.mxords = 100;
    opt.mxords = Math.min(opt.mxords, mord[2]);
  }

  if (opt.hmax < 0.) {
    ctx.error = '[lsoda] hmax < 0.';
    return false;
  }
  opt.hmxi = 0.;
  if (opt.hmax > 0)
    opt.hmxi = 1. / opt.hmax;
  if (opt.hmin < 0.) {
    ctx.error = '[lsoda] hmin < 0.';
    return false;
  }
  return true;
}

function allocMem(ctx: LsodaContext): boolean {
  const c = ctx.common!;
  const nyh = ctx.neq;
  const lenyh = 1 + Math.max(ctx.opt!.mxordn, ctx.opt!.mxords);

  // yh: array of (lenyh+1) Float64Arrays, each of size (nyh+1)
  c.yh = new Array(lenyh + 1);
  for (let i = 0; i <= lenyh; i++)
    c.yh[i] = new Float64Array(nyh + 1);

  // wm: array of (nyh+1) Float64Arrays, each of size (nyh+1)
  c.wm = new Array(nyh + 1);
  for (let i = 0; i <= nyh; i++)
    c.wm[i] = new Float64Array(nyh + 1);

  c.ewt = new Float64Array(nyh + 1);
  c.savf = new Float64Array(nyh + 1);
  c.acor = new Float64Array(nyh + 1);
  c.ipvt = new Int32Array(nyh + 1);

  // Allocate elco and tesco
  c.elco = new Array(13);
  for (let i = 0; i < 13; i++)
    c.elco[i] = new Float64Array(14);
  c.tesco = new Array(13);
  for (let i = 0; i < 13; i++)
    c.tesco[i] = new Float64Array(4);

  return true;
}

// ---- Public API ----

/**
 * Prepare a context with options. Must be called before lsoda().
 * Returns true on success.
 */
export function lsodaPrepare(ctx: LsodaContext, opt: LsodaOpt): boolean {
  ctx.common = new LsodaCommon();
  ctx.opt = opt;
  if (!checkOpt(ctx, opt)) return false;
  return allocMem(ctx);
}

/**
 * Reset context for reuse with a different problem.
 * All options are preserved, internal state is zeroed.
 */
export function lsodaReset(ctx: LsodaContext): void {
  const c = ctx.common!;
  // Zero all scalar state
  c.h = 0; c.hu = 0; c.rc = 0; c.tn = 0; c.tsw = 0; c.pdnorm = 0;
  c.crate = 0; c.hold = 0; c.rmax = 0; c.pdest = 0; c.pdlast = 0;
  c.ialth = 0; c.ipup = 0; c.nslp = 0; c.icount = 0; c.irflag = 0;
  c.imxer = 0; c.illin = 0; c.nhnil = 0; c.nslast = 0; c.jcur = 0;
  c.meth = 0; c.mused = 0; c.nq = 0; c.nst = 0; c.ncf = 0;
  c.nfe = 0; c.nje = 0; c.nqu = 0; c.miter = 0;
  c.el.fill(0);
  // Zero arrays
  for (let i = 0; i < c.yh.length; i++) c.yh[i].fill(0);
  for (let i = 0; i < c.wm.length; i++) c.wm[i].fill(0);
  c.ewt.fill(0); c.savf.fill(0); c.acor.fill(0); c.ipvt.fill(0);
  for (let i = 0; i < c.elco.length; i++) c.elco[i].fill(0);
  for (let i = 0; i < c.tesco.length; i++) c.tesco[i].fill(0);
  ctx.state = 1;
  ctx.error = null;
}

/**
 * Free context resources. In JS/TS, this is mostly a no-op (GC handles it).
 */
export function lsodaFree(ctx: LsodaContext): void {
  if (ctx.error)
    console.error(`unhandled error message: ${ctx.error}`);

  ctx.common = null;
}

/**
 * Main solver. Integrates from t towards tout.
 * y: 1-indexed Float64Array (y[0] unused, y[1..neq] are state values)
 * Returns { t, state } — updated time and solver state.
 */
export function lsoda(
  ctx: LsodaContext, y: Float64Array, tIn: number, tout: number,
): { t: number; state: number } {
  const c = ctx.common!;
  const opt = ctx.opt!;
  const neq = ctx.neq;
  let t = tIn;

  let kflag: number;
  let jstart: number = 0;
  let ihit = 0;

  // Helper functions matching C macros
  function hardfailure(msg: string): { t: number; state: number } {
    ctx.error = msg;
    ctx.state = -3;
    return {t, state: ctx.state};
  }

  function softfailure(code: number, msg: string): { t: number; state: number } {
    ctx.error = msg;
    for (let i = 1; i <= neq; i++)
      y[i] = c.yh[1][i];
    t = c.tn;
    ctx.state = code;
    return {t, state: ctx.state};
  }

  function successreturn(): { t: number; state: number } {
    for (let i = 1; i <= neq; i++)
      y[i] = c.yh[1][i];
    t = c.tn;
    if (itask === 4 || itask === 5) {
      if (ihit)
        t = tcrit;
    }
    ctx.state = 2;
    return {t, state: ctx.state};
  }

  function intdyreturn(): { t: number; state: number } {
    const iflag = intdy(ctx, tout, 0, y);
    if (iflag !== 0) {
      ctx.error = `[lsoda] trouble from intdy, itask = ${itask}, tout = ${tout}`;
      for (let i = 1; i <= neq; i++)
        y[i] = c.yh[1][i];
      t = c.tn;
    }
    t = tout;
    ctx.state = 2;
    return {t, state: ctx.state};
  }

  if (c === null)
    return hardfailure('[lsoda] illegal common block did you call lsoda_prepare?');


  // Block b: initial/continuation call setup
  let h0: number = 0;
  let tcrit: number = 0;
  const rtol = opt.rtol;
  const atol = opt.atol;

  if (ctx.state === 1 || ctx.state === 3) {
    h0 = opt.h0;
    if (ctx.state === 1) {
      if ((tout - t) * h0 < 0.)
        return hardfailure(`[lsoda] tout = ${tout} behind t = ${t}. integration direction is given by ${h0}`);
    }
  }

  const itask = opt.itask;

  if (ctx.state === 3)
    jstart = -1;


  // Block c: initial call only (state == 1)
  if (ctx.state === 1) {
    c.meth = 1;
    c.tn = t;
    c.tsw = t;

    if (itask === 4 || itask === 5) {
      tcrit = opt.tcrit;
      if ((tcrit - tout) * (tout - t) < 0.)
        return hardfailure('[lsoda] itask = 4 or 5 and tcrit behind tout');

      if (h0 !== 0. && (t + h0 - tcrit) * h0 > 0.)
        h0 = tcrit - t;
    }

    jstart = 0;
    c.nq = 1;

    // Initial call to f
    ctx.func(t, y.subarray(1), c.yh[2].subarray(1), ctx.data);
    c.nfe = 1;

    // Load initial value vector
    for (let i = 1; i <= neq; i++)
      c.yh[1][i] = y[i];

    // Load and invert ewt
    ewset(ctx, y);
    for (let i = 1; i <= neq; i++) {
      if (c.ewt[i] <= 0.)
        return hardfailure(`[lsoda] ewt[${i}] = ${c.ewt[i]} <= 0.`);
    }

    // Compute initial step size h0 if not given
    if (h0 === 0.) {
      const tdist = Math.abs(tout - t);
      const w0 = Math.max(Math.abs(t), Math.abs(tout));
      if (tdist < 2. * ETA * w0)
        return hardfailure('[lsoda] tout too close to t to start integration');

      let tol = 0.;
      for (let i = 1; i <= neq; i++)
        tol = Math.max(tol, rtol[i]);
      if (tol <= 0.) {
        for (let i = 1; i <= neq; i++) {
          const atoli = atol[i];
          const ayi = Math.abs(y[i]);
          if (ayi !== 0.)
            tol = Math.max(tol, atoli / ayi);
        }
      }
      tol = Math.max(tol, 100. * ETA);
      tol = Math.min(tol, 0.001);
      const sum = vmnorm(neq, c.yh[2], c.ewt);
      const sumSq = 1. / (tol * w0 * w0) + tol * sum * sum;
      h0 = 1. / Math.sqrt(sumSq);
      h0 = Math.min(h0, tdist);
      h0 = h0 * ((tout - t >= 0.) ? 1. : -1.);
    }

    // Adjust h0 for hmax
    const rh = Math.abs(h0) * opt.hmxi;
    if (rh > 1.) h0 /= rh;

    // Load h and scale yh[2]
    c.h = h0;
    for (let i = 1; i <= neq; i++)
      c.yh[2][i] *= h0;
  }

  // Block d: continuation calls (state == 2 or 3)
  if (ctx.state === 2 || ctx.state === 3) {
    jstart = 1;
    c.nslast = c.nst;

    switch (itask) {
    case 1:
      if ((c.tn - tout) * c.h >= 0.)
        return intdyreturn();

      break;
    case 2:
      break;
    case 3: {
      const tp = c.tn - c.hu * (1. + 100. * ETA);
      if ((tp - tout) * c.h > 0.)
        return hardfailure(`[lsoda] itask = ${itask} and tout behind tcur - hu`);

      if ((c.tn - tout) * c.h < 0.) break;
      return successreturn();
    }
    case 4:
      tcrit = opt.tcrit;
      if ((c.tn - tcrit) * c.h > 0.)
        return hardfailure('[lsoda] itask = 4 or 5 and tcrit behind tcur');

      if ((tcrit - tout) * c.h < 0.)
        return hardfailure('[lsoda] itask = 4 or 5 and tcrit behind tout');

      if ((c.tn - tout) * c.h >= 0.)
        return intdyreturn();

      // fall through to case 5 logic
      // eslint-disable-next-line no-fallthrough
    case 5:
      if (itask === 5) {
        tcrit = opt.tcrit;
        if ((c.tn - tcrit) * c.h > 0.)
          return hardfailure('[lsoda] itask = 4 or 5 and tcrit behind tcur');
      }
      {
        const hmx = Math.abs(c.tn) + Math.abs(c.h);
        ihit = Math.abs(c.tn - tcrit) <= (100. * ETA * hmx) ? 1 : 0;
        if (ihit) {
          t = tcrit;
          return successreturn();
        }
        const tnext = c.tn + c.h * (1. + 4. * ETA);
        if ((tnext - tcrit) * c.h <= 0.) break;
        c.h = (tcrit - c.tn) * (1. - 4. * ETA);
        if (ctx.state === 2) jstart = -2;
      }
      break;
    }
  }

  // Block e: main integration loop
  while (true) {
    if (ctx.state !== 1 || c.nst !== 0) {
      if ((c.nst - c.nslast) >= opt.mxstep)
        return softfailure(-1, `[lsoda] ${opt.mxstep} steps taken before reaching tout`);

      ewset(ctx, c.yh[1]);
      for (let i = 1; i <= neq; i++) {
        if (c.ewt[i] <= 0.)
          return softfailure(-6, `[lsoda] ewt[${i}] = ${c.ewt[i]} <= 0.`);
      }
    }

    const tolsf = ETA * vmnorm(neq, c.yh[1], c.ewt);
    if (tolsf > 0.01) {
      const scaled = tolsf * 200.;
      if (c.nst === 0) {
        return hardfailure(
          'lsoda -- at start of problem, too much accuracy requested' +
          ` for precision of machine, suggested scaling factor = ${scaled}`,
        );
      }
      return softfailure(-2,
        `lsoda -- at t = ${t}, too much accuracy requested` +
        ` for precision of machine, suggested scaling factor = ${scaled}`,
      );
    }

    if ((c.tn + c.h) === c.tn) {
      c.nhnil++;
      if (c.nhnil <= opt.mxhnil) {
        console.error(`lsoda -- warning..internal t = ${c.tn} and h = ${c.h} are` +
          ' such that t + h = t on the next step');
        if (c.nhnil === opt.mxhnil) {
          console.error(`lsoda -- above warning has been issued ${c.nhnil}` +
            ' times, it will not be issued again for this problem');
        }
      }
    }

    // Call stoda
    kflag = stoda(ctx, y, jstart);

    if (kflag === 0) {
      // Block f: successful return from stoda
      // Dense output: capture snapshot
      if (ctx.snapshots)
        ctx.snapshots.push(captureSnapshot(c, neq));

      jstart = 1;
      if (c.meth !== c.mused) {
        c.tsw = c.tn;
        jstart = -1;
        if (opt.ixpr) {
          if (c.meth === 2) {
            console.error('[lsoda] a switch to the stiff method has occurred' +
              ` at t = ${c.tn}, tentative step size h = ${c.h}, step nst = ${c.nst}`);
          }
          if (c.meth === 1) {
            console.error('[lsoda] a switch to the nonstiff method has occurred' +
              ` at t = ${c.tn}, tentative step size h = ${c.h}, step nst = ${c.nst}`);
          }
        }
      }

      if (itask === 1) {
        if ((c.tn - tout) * c.h < 0.) continue;
        return intdyreturn();
      }
      if (itask === 2)
        return successreturn();

      if (itask === 3) {
        if ((c.tn - tout) * c.h >= 0.)
          return successreturn();

        continue;
      }
      if (itask === 4) {
        tcrit = opt.tcrit;
        if ((c.tn - tout) * c.h >= 0.)
          return intdyreturn();
        else {
          const hmx = Math.abs(c.tn) + Math.abs(c.h);
          ihit = Math.abs(c.tn - tcrit) <= (100. * ETA * hmx) ? 1 : 0;
          if (ihit)
            return successreturn();

          const tnext = c.tn + c.h * (1. + 4. * ETA);
          if ((tnext - tcrit) * c.h <= 0.) continue;
          c.h = (tcrit - c.tn) * (1. - 4. * ETA);
          jstart = -2;
          continue;
        }
      }
      if (itask === 5) {
        tcrit = opt.tcrit;
        const hmx = Math.abs(c.tn) + Math.abs(c.h);
        ihit = Math.abs(c.tn - tcrit) <= (100. * ETA * hmx) ? 1 : 0;
        return successreturn();
      }
    }

    // kflag == -1 or -2: error test / convergence failures
    if (kflag === -1 || kflag === -2) {
      let big = 0.;
      c.imxer = 1;
      for (let i = 1; i <= neq; i++) {
        const size = Math.abs(c.acor[i]) * c.ewt[i];
        if (big < size) {
          big = size;
          c.imxer = i;
        }
      }
      if (kflag === -1) {
        return softfailure(-4,
          `lsoda -- at t = ${c.tn} and step size h = ${c.h},` +
          ' the error test failed repeatedly or with abs(h) = hmin',
        );
      }
      if (kflag === -2) {
        return softfailure(-5,
          `lsoda -- at t = ${c.tn} and step size h = ${c.h},` +
          ' the corrector convergence failed repeatedly or with abs(h) = hmin',
        );
      }
    }
  }
}

// ---- High-level convenience API ----

/**
 * User-facing ODE function signature (0-based arrays).
 */
export type OdeFunction = (t: number, y: Float64Array, ydot: Float64Array, data?: any) => number;

export interface LsodaSolveResult {
  y: Float64Array;
  t: number;
}

/**
 * High-level LSODA solver class with 0-based public interface.
 */
export class Lsoda {
  private ctx: LsodaContext;
  private internalY: Float64Array;
  private denseEnabled: boolean = false;

  constructor(
    f: OdeFunction, neq: number,
    opt?: Partial<LsodaOpt & { dense: boolean }>, data?: any,
  ) {
    // Wrap user's 0-based function to internal 1-based convention
    const wrappedFunc: LsodaFunc = (t, y, ydot, data) => {
      return f(t, y, ydot, data);
    };

    this.ctx = new LsodaContext(wrappedFunc, neq, data);
    this.internalY = new Float64Array(neq + 1); // 1-indexed

    // Build full options with defaults
    const fullOpt: LsodaOpt = {
      ixpr: 0, mxstep: 0, mxhnil: 0, mxordn: 0, mxords: 0,
      tcrit: 0, h0: 0, hmax: 0, hmin: 0, hmxi: 0, itask: 1,
      rtol: new Float64Array(neq + 1),
      atol: new Float64Array(neq + 1),
      ...opt,
    };

    // If rtol/atol provided as partial, make sure they're 1-indexed Float64Arrays
    if (opt?.rtol) {
      if (opt.rtol.length === neq) {
        // User provided 0-based, convert to 1-based
        const r = new Float64Array(neq + 1);
        for (let i = 0; i < neq; i++) r[i + 1] = opt.rtol[i];
        fullOpt.rtol = r;
      } else
        fullOpt.rtol = opt.rtol;
    }
    if (opt?.atol) {
      if (opt.atol.length === neq) {
        const a = new Float64Array(neq + 1);
        for (let i = 0; i < neq; i++) a[i + 1] = opt.atol[i];
        fullOpt.atol = a;
      } else
        fullOpt.atol = opt.atol;
    }

    lsodaPrepare(this.ctx, fullOpt);

    if (opt?.dense) {
      this.denseEnabled = true;
      this.ctx.snapshots = [];
    }
  }

  /**
   * Integrate from current t to tout.
   * y: 0-based array of state values.
   * Returns updated { y, t }.
   */
  solve(y: ArrayLike<number>, t: number, tout: number): LsodaSolveResult {
    const neq = this.ctx.neq;
    // Copy 0-based user y into 1-based internal array
    for (let i = 0; i < neq; i++)
      this.internalY[i + 1] = y[i];

    const result = lsoda(this.ctx, this.internalY, t, tout);

    // Copy back to 0-based
    const out = new Float64Array(neq);
    for (let i = 0; i < neq; i++)
      out[i] = this.internalY[i + 1];

    return {y: out, t: result.t};
  }

  get state(): number {return this.ctx.state;}
  get error(): string | null {return this.ctx.error;}

  /** Returns a DenseOutput interpolator from collected snapshots. */
  getDenseOutput(): DenseOutput {
    if (!this.ctx.snapshots || this.ctx.snapshots.length === 0)
      throw new Error('[Lsoda] dense output not enabled or no steps taken — pass { dense: true } in options');
    return new DenseOutput(this.ctx.snapshots, this.ctx.neq);
  }

  reset(): void {
    lsodaReset(this.ctx);
    if (this.denseEnabled)
      this.ctx.snapshots = [];
  }
}
