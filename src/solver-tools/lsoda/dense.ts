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

import {NordsieckSnapshot} from './common';

/**
 * Dense output interpolator. Given a sequence of Nordsieck snapshots
 * captured during integration, evaluates the solution at arbitrary points
 * via Taylor/Horner interpolation.
 */
export class DenseOutput {
  private snaps: NordsieckSnapshot[];
  private neq: number;

  constructor(snapshots: NordsieckSnapshot[], neq: number) {
    if (snapshots.length === 0)
      throw new Error('[DenseOutput] no snapshots collected');
    this.snaps = snapshots;
    this.neq = neq;
  }

  /** Leftmost reachable time. */
  get tMin(): number {
    const s = this.snaps[0];
    return s.tn - s.hu;
  }

  /** Rightmost reachable time. */
  get tMax(): number {
    return this.snaps[this.snaps.length - 1].tn;
  }

  /**
   * Evaluate the solution at a single time point.
   * Uses Horner's method on the Nordsieck array of the enclosing step.
   */
  evaluateAt(t: number): Float64Array {
    const snap = this.findSnap(t);
    return this.horner(snap, t, new Float64Array(this.neq));
  }

  /**
   * Evaluate the solution at each point in a sorted ascending array.
   * Uses a linear scan over snapshots for efficiency.
   */
  solveAtTimes(tArray: Float64Array): Float64Array[] {
    if (tArray.length === 0) return [];
    const snaps = this.snaps;

    // Validate range
    const tFirst = tArray[0];
    const tLast = tArray[tArray.length - 1];
    if (tFirst < this.tMin - 1e-14 * Math.abs(this.tMin))
      throw new Error(`[DenseOutput] t=${tFirst} is before tMin=${this.tMin}`);
    if (tLast > this.tMax + 1e-14 * Math.abs(this.tMax))
      throw new Error(`[DenseOutput] t=${tLast} is after tMax=${this.tMax}`);

    const dky = new Float64Array(this.neq);
    const result: Float64Array[] = [];

    for (let i = 0; i < this.neq; ++i)
      result.push(new Float64Array(tArray.length));

    let si = 0; // snapshot pointer

    for (let qi = 0; qi < tArray.length; qi++) {
      const t = tArray[qi];
      // Advance snapshot pointer until we find an interval containing t
      while (si < snaps.length - 1 && t > snaps[si].tn + 1e-14 * Math.abs(snaps[si].tn))
        si++;

      this.horner(snaps[si], t, dky);

      for (let i = 0; i < this.neq; ++i)
        result[i][qi] = dky[i];
    }
    return result;
  }

  /**
   * Evaluate the solution on a uniform grid from tStart to tEnd with given step.
   */
  solveOnGrid(tStart: number, tEnd: number, step: number): { t: Float64Array, y: Float64Array[] } {
    if (step <= 0)
      throw new Error('[DenseOutput] step must be positive');

    const n = Math.floor((tEnd - tStart) / step) + 1;
    const tArr = new Float64Array(n);
    for (let i = 0; i < n - 1; i++)
      tArr[i] = tStart + i * step;
    tArr[n - 1] = tEnd;

    const t = tArr;
    const y = this.solveAtTimes(t);
    return {t, y};
  }

  /** Find the snapshot whose interval contains t via binary search. */
  private findSnap(t: number): NordsieckSnapshot {
    const snaps = this.snaps;

    // Quick boundary checks
    if (t <= snaps[0].tn) return snaps[0];
    if (t >= snaps[snaps.length - 1].tn) return snaps[snaps.length - 1];

    // Binary search: find first snap with tn >= t
    let lo = 0; let hi = snaps.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (snaps[mid].tn < t)
        lo = mid + 1;
      else
        hi = mid;
    }
    return snaps[lo];
  }

  /** Horner evaluation of the Nordsieck polynomial at time t. */
  private horner(snap: NordsieckSnapshot, t: number, dky: Float64Array): Float64Array {
    const {nq, h, yh} = snap;
    const neq = this.neq;
    const s = (t - snap.tn) / h;

    for (let i = 0; i < neq; i++) {
      let val = yh[nq][i];
      for (let j = nq - 1; j >= 0; j--)
        val = yh[j][i] + s * val;
      dky[i] = val;
    }
    return dky;
  }
}
