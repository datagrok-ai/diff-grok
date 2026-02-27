/**
 * TypeScript port of SUNDIALS CVODE v7.5.0 (https://github.com/LLNL/sundials)
 *
 * Copyright (c) 2025, Lawrence Livermore National Security, University of
 *   Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security and
 *   Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * Copyright (c) 2026 Datagrok.
 *
 * Licensed under the BSD-3-Clause License. See THIRD_PARTY_LICENSES
 * for the full license text and provenance chain of the original code.
 */

// CVODE TypeScript port - High-level Cvode class API
// Follows the Lsoda class pattern from src/lsoda.ts

import {
  CvodeMem,
  CvodeRhsFn,
  CvodeJacFn,
  CvodeRootFn,
  cvodeGetDky,
  CV_ADAMS, CV_BDF,
  CV_NORMAL, CV_ONE_STEP,
  CV_ROOT_RETURN,
} from './common';

import {cvodeCreate, cvodeInit, cvode} from './cvode';
import {cvodeSetLinearSolver, cvodeSetJacFn} from './cvode_ls';
import {cvodeRootInit} from './cvode_root';
import {
  CvodeStats,
  cvodeGetIntegratorStats,
  cvodeSStolerances,
  cvodeSVtolerances,
  cvodeSetMaxNumSteps,
  cvodeSetMaxOrd,
  cvodeSetMaxStep,
  cvodeSetMinStep,
  cvodeSetInitStep,
  cvodeSetStopTime,
  cvodeSetUserData,
} from './cvode_io';

// ============================================================================
// Public types
// ============================================================================

/**
 * User ODE function: y' = f(t, y).
 * Arrays are 0-based. Return value is void; wrapping to int is done internally.
 */
export type OdeFunction = (t: number, y: Float64Array, ydot: Float64Array) => void;

/**
 * Options for the Cvode solver.
 */
export interface CvodeOptions {
  /** Linear multistep method: 'adams' (non-stiff) or 'bdf' (stiff). Default: 'bdf'. */
  lmm?: 'adams' | 'bdf';
  /** Relative tolerance. Default: 1e-4. */
  rtol?: number;
  /** Absolute tolerance: scalar or per-component array. Default: 1e-6. */
  atol?: number | Float64Array;
  /** Maximum number of internal steps per cvode() call. */
  maxSteps?: number;
  /** Maximum method order. */
  maxOrder?: number;
  /** Maximum step size (0 = no limit). */
  maxStep?: number;
  /** Minimum step size. */
  minStep?: number;
  /** Initial step size (0 = auto-select). */
  initStep?: number;
  /** Stop time: integration will not proceed past this value. */
  stopTime?: number;
  /** Root function: g(t, y, gout). */
  rootFn?: (t: number, y: Float64Array, gout: Float64Array) => void;
  /** Number of root functions (required when rootFn is provided). */
  nRootFns?: number;
  /** User-supplied Jacobian function: J = df/dy. */
  jacFn?: (t: number, y: Float64Array, fy: Float64Array, J: Float64Array[]) => void;
  /** Arbitrary user data (passed as context, not forwarded to callbacks). */
  userData?: any;
}

/**
 * Result returned by solve() and step().
 */
export interface CvodeSolveResult {
  /** Time reached by the solver. */
  t: number;
  /** Solution vector at time t (0-based, copy). */
  y: Float64Array;
  /** Solver flag (CV_SUCCESS=0, CV_ROOT_RETURN=2, CV_TSTOP_RETURN=1, negative=error). */
  flag: number;
  /** Root indices found when flag === CV_ROOT_RETURN. */
  rootsFound?: Int32Array;
}

// ============================================================================
// Cvode class
// ============================================================================

/**
 * High-level CVODE solver class with a 0-based public interface.
 *
 * Usage:
 *   const solver = new Cvode(f, neq, t0, y0, { lmm: 'bdf', rtol: 1e-6, atol: 1e-9 });
 *   const result = solver.solve(tout);
 */
export class Cvode {
  private mem: CvodeMem;
  private neq: number;
  /** Output vector reused across calls (0-based, neq elements). */
  private yout: Float64Array;
  /** Wrapped RHS stored to allow re-initialization. */
  private wrappedF: CvodeRhsFn;

  constructor(
    f: OdeFunction,
    neq: number,
    t0: number,
    y0: Float64Array,
    options?: CvodeOptions,
  ) {
    this.neq = neq;
    this.yout = new Float64Array(neq);

    const opts = options ?? {};

    // --- 1. Create CvodeMem ---
    const lmmValue = (opts.lmm === 'adams') ? CV_ADAMS : CV_BDF;
    this.mem = cvodeCreate(lmmValue);

    // --- 2. Wrap user's OdeFunction to CvodeRhsFn (void -> number) ---
    this.wrappedF = (_t: number, _y: Float64Array, _ydot: Float64Array, _userData: any): number => {
      f(_t, _y, _ydot);
      return 0;
    };

    // --- 3. Initialize ---
    cvodeInit(this.mem, this.wrappedF, t0, y0);

    // --- 4. Set tolerances ---
    const rtol = opts.rtol ?? 1e-4;
    const atol = opts.atol ?? 1e-6;

    if (typeof atol === 'number')
      cvodeSStolerances(this.mem, rtol, atol);
    else
      cvodeSVtolerances(this.mem, rtol, atol);


    // --- 5. Set up dense linear solver (required for BDF; optional but useful for Adams) ---
    cvodeSetLinearSolver(this.mem);

    // --- 6. Set user-supplied Jacobian if provided ---
    if (opts.jacFn) {
      const userJac = opts.jacFn;
      const wrappedJacFn: CvodeJacFn = (
        t: number, y: Float64Array, fy: Float64Array,
        J: Float64Array[], userData: any,
      ): number => {
        userJac(t, y, fy, J);
        return 0;
      };
      cvodeSetJacFn(this.mem, wrappedJacFn);
    }

    // --- 7. Apply optional settings ---
    if (opts.maxSteps !== undefined) cvodeSetMaxNumSteps(this.mem, opts.maxSteps);
    if (opts.maxOrder !== undefined) cvodeSetMaxOrd(this.mem, opts.maxOrder);
    if (opts.maxStep !== undefined) cvodeSetMaxStep(this.mem, opts.maxStep);
    if (opts.minStep !== undefined) cvodeSetMinStep(this.mem, opts.minStep);
    if (opts.initStep !== undefined) cvodeSetInitStep(this.mem, opts.initStep);
    if (opts.stopTime !== undefined) cvodeSetStopTime(this.mem, opts.stopTime);
    if (opts.userData !== undefined) cvodeSetUserData(this.mem, opts.userData);

    // --- 8. Set up rootfinding if requested ---
    if (opts.rootFn && opts.nRootFns && opts.nRootFns > 0) {
      const userRoot = opts.rootFn;
      const wrappedRootFn: CvodeRootFn = (
        t: number, y: Float64Array, gout: Float64Array, _userData: any,
      ): number => {
        userRoot(t, y, gout);
        return 0;
      };
      cvodeRootInit(this.mem, opts.nRootFns, wrappedRootFn);
    }
  }

  // ==========================================================================
  // Public methods
  // ==========================================================================

  /**
   * Integrate from the current time to tout using CV_NORMAL mode
   * (integrator advances to tout and interpolates).
   *
   * @param tout - target output time
   * @returns CvodeSolveResult with { t, y, flag, rootsFound? }
   */
  solve(tout: number): CvodeSolveResult {
    const {flag, t} = cvode(this.mem, tout, this.yout, CV_NORMAL);
    return this._buildResult(t, flag);
  }

  /**
   * Take a single internal step using CV_ONE_STEP mode.
   * The solver advances by one internal step of whatever size it chooses.
   *
   * @returns CvodeSolveResult with { t, y, flag, rootsFound? }
   */
  step(): CvodeSolveResult {
    const {flag, t} = cvode(this.mem, 1e30, this.yout, CV_ONE_STEP);
    return this._buildResult(t, flag);
  }

  /**
   * Dense output: compute the k-th derivative of the interpolating polynomial
   * at time t using the current Nordsieck history array.
   *
   * @param t - interpolation time (must be within the last step interval)
   * @param k - derivative order (0 = solution, 1 = first derivative, ...)
   * @returns copy of the interpolated vector (length neq)
   */
  getDky(t: number, k: number): Float64Array {
    const dky = new Float64Array(this.neq);
    const flag = cvodeGetDky(this.mem, t, k, dky);
    if (flag !== 0)
      throw new Error(`[Cvode] getDky failed: t=${t} out of range or k=${k} invalid`);

    return dky;
  }

  /**
   * Return current integrator statistics.
   */
  getStats(): CvodeStats {
    return cvodeGetIntegratorStats(this.mem);
  }

  /**
   * Re-initialize the solver for a new initial-value problem with the same
   * ODE function and options. Resets internal state and counters.
   *
   * @param t0 - new initial time
   * @param y0 - new initial conditions (length neq)
   */
  reInit(t0: number, y0: Float64Array): void {
    cvodeInit(this.mem, this.wrappedF, t0, y0);
  }

  // ==========================================================================
  // Private helpers
  // ==========================================================================

  private _buildResult(t: number, flag: number): CvodeSolveResult {
    const y = new Float64Array(this.yout); // copy
    const result: CvodeSolveResult = {t, y, flag};

    if (flag === CV_ROOT_RETURN) {
      // Return a copy of the root direction info array
      result.rootsFound = new Int32Array(this.mem.cv_iroots);
    }

    return result;
  }
}
