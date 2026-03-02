# CVODE C-to-TypeScript Transformation Plan

This document describes the plan for transforming the CVODE solver (v7.5.0) from
C to TypeScript, following the patterns established in `LSODA-C-to-TS.md`.

**Source code location:** All TypeScript source files go into `src/cvode/`.

**Test files location:** All tests go into `test/`:
- `test/cvode.test.ts` — correctness tests (mirrors `test/test.test.ts`)
- `test/cvode.perf.test.ts` — benchmark/performance tests (mirrors `test/perf.test.ts`)

---

## 1. Scope

### 1.1 What to include

| Feature | Source | Notes |
|---------|--------|-------|
| BDF method (stiff) | `cvode.c` | Orders 1-5, primary use case |
| Adams-Moulton method (non-stiff) | `cvode.c` | Orders 1-12 |
| Nordsieck history array | `cvode.c` | Core representation, same concept as LSODA |
| Newton nonlinear solver | `cvode_nls.c` | Default for BDF/stiff |
| Fixed-point nonlinear solver | `cvode_nls.c` | For Adams/non-stiff |
| Dense direct linear solver | `cvode_ls.c` + `sunlinsol_dense` | LU with partial pivoting |
| Difference-quotient Jacobian | `cvode_ls.c` (`cvLsDQJac`) | Default when no user Jacobian |
| User-supplied Jacobian | `cvode_ls.c` | Optional callback |
| Error control & step/order selection | `cvode.c` | WRMS norm, adaptive h and q |
| Rootfinding | `cvode.c` (`cvRcheck*`, `cvRootfind`) | Zero-crossing detection |
| Dense output | `cvode.c` (`CVodeGetDky`) | Horner evaluation on Nordsieck array |
| BDF stability limit detection | `cvode.c` (`cvBDFStab`) | Optional safety feature |
| I/O: tolerances, step limits, stats | `cvode_io.c` (subset) | Configuration and statistics getters |

### 1.2 What to exclude

| Feature | Reason |
|---------|--------|
| Band matrix / band solver | Specialized; add later if needed |
| Sparse matrix / KLU / SuperLU | Specialized; add later if needed |
| Iterative solvers (SPGMR, SPBCGS, SPTFQMR, SPFGMR) | Krylov methods for large systems; add later if needed |
| Band preconditioner (`cvode_bandpre.c`) | Requires band solver |
| BBD preconditioner (`cvode_bbdpre.c`) | Requires band solver + MPI concepts |
| Diagonal solver (`cvode_diag.c`) | Rarely used; add later if needed |
| Projection (`cvode_proj.c`) | Constraint projection; add later if needed |
| Problem resizing (`cvode_resize.c`) | Unusual use case; add later if needed |
| Constraint enforcement (`cvode_constraints.c`) | Add later if needed |
| Fused kernel stubs (`cvode_fused_stubs.c`) | Performance optimization not relevant to TS |
| CLI utilities (`cvode_cli.c`) | Not applicable |
| Parallel N_Vector variants (OpenMP, CUDA, etc.) | Single-threaded TS runtime |
| SUNContext / SUNDIALS context | Heavyweight infrastructure; replace with simple config |

### 1.3 Size estimate

| C source (in scope) | Approx lines | TS estimate |
|----------------------|--------------|-------------|
| `cvode.c` (core step logic, rootfinding, dense output) | ~4,900 | ~2,500 |
| `cvode_nls.c` (nonlinear solver) | ~430 | ~250 |
| `cvode_ls.c` (linear solver, DQ Jacobian) | ~1,950 | ~800 |
| `cvode_io.c` (subset: tolerances, step config, stats) | ~600 | ~300 |
| `cvode_impl.h` (data structures, constants) | ~790 | ~200 |
| Dense linear algebra (LU factorize + solve) | ~300 | ~200 |
| **Total** | **~8,970** | **~4,250** |

Reduction factors: no `#include` guards, no memory management boilerplate, no
SUNDIALS abstraction layers, consolidated declarations, and TypeScript's more
concise syntax.

---

## 2. File & Module Mapping

### 2.1 C files → TypeScript modules

| C source | TS module | Purpose |
|----------|-----------|---------|
| `cvode_impl.h` | `src/cvode/common.ts` | `CvodeMem` class, constants, types, enums |
| `cvode.c` (core loop) | `src/cvode/cvode.ts` | `cvodeCreate`, `cvodeInit`, `cvode` (main entry), `cvStep`, predict/correct/error test |
| `cvode.c` (BDF coefficients) | `src/cvode/cvode_bdf.ts` | `cvSetBDF`, BDF coefficient tables |
| `cvode.c` (Adams coefficients) | `src/cvode/cvode_adams.ts` | `cvSetAdams`, Adams coefficient tables |
| `cvode.c` (initial step) | `src/cvode/cvode_hin.ts` | `cvHin` — initial step size estimation |
| `cvode.c` (rootfinding) | `src/cvode/cvode_root.ts` | `cvRcheck1/2/3`, `cvRootfind` — zero-crossing detection |
| `cvode_nls.c` | `src/cvode/cvode_nls.ts` | Nonlinear solver: residual, convergence test, Newton/FP wrappers |
| `cvode_ls.c` | `src/cvode/cvode_ls.ts` | Linear solver interface: Jacobian setup/solve, DQ Jacobian |
| `cvode_io.c` (subset) | `src/cvode/cvode_io.ts` | Tolerance setup, optional inputs, statistics getters |
| `sunlinsol_dense` + `sunmatrix_dense` | `src/cvode/dense_linalg.ts` | Dense LU factorization (`dgefa`) and solve (`dgesl`) |
| *(none)* | `src/cvode/cvode_class.ts` | High-level `Cvode` class (idiomatic TS API) |
| *(none)* | `src/cvode/index.ts` | Barrel re-exports for public API |

**Total: 12 TypeScript files** (vs ~10 C files in scope).

### 2.2 Rationale for splits

CVODE's `cvode.c` is ~4,900 lines — too large for a single TS module. Splitting
by functional area (BDF coefficients, Adams coefficients, initial step size,
rootfinding) produces manageable modules of 200-600 lines each, matching the
LSODA pattern where `stoda.ts`, `cfode.ts`, `methodswitch.ts`, etc. are separate
files.

The dense linear algebra is extracted into its own module (`dense_linalg.ts`)
because CVODE's SUNDIALS infrastructure wraps it in multiple abstraction layers
(SUNMatrix → SUNLinearSolver → CVLsMem). In TypeScript we flatten this: one
module with `dgefa` (LU factorize) and `dgesl` (LU solve) operating directly on
`Float64Array[]`.

### 2.3 Reuse from existing LSODA code

The existing `src/blas.ts` already contains `dgefa`, `dgesl`, `daxpy`, `ddot`,
`dscal`, `idamax` — the same BLAS routines CVODE needs. **Reuse `src/blas.ts`
directly** rather than creating `dense_linalg.ts` as a separate module if the
signatures are compatible.

**Check:** CVODE's `cvLsDQJac` and `cvLsSetup` use column-major dense matrix
indexing (`SM_ELEMENT_D(A, i, j)` = `A_data[j * M + i]`). The existing
`blas.ts` uses `Float64Array[]` (array of rows). **Decision:** Use the same
row-based `Float64Array[]` layout as LSODA's `blas.ts`. Translate CVODE's
column-major access patterns to row-major. This keeps the codebase consistent
and allows full `blas.ts` reuse.

If the existing `blas.ts` functions need minor signature adjustments (e.g.,
CVODE passes different offset patterns), extend rather than fork.

---

## 3. Data Structures

### 3.1 CvodeMem — Main solver state

Maps from `CVodeMemRec` in `cvode_impl.h`. Follow LSODA's `LsodaCommon`
pattern: mutable class with all fields initialized to defaults.

```typescript
// src/cvode/common.ts

export class CvodeMem {
  // --- Problem specification ---
  cv_f: CvodeRhsFn | null = null;     // RHS function: y' = f(t,y)
  cv_user_data: any = null;            // user data
  cv_lmm: number = 0;                 // CV_ADAMS or CV_BDF

  // --- Nordsieck history array ---
  cv_zn: Float64Array[] = [];          // zn[0..q]: Nordsieck array
  cv_ewt: Float64Array = new Float64Array(0);  // error weight vector

  // --- Step data ---
  cv_q: number = 0;                    // current method order
  cv_qprime: number = 0;              // order for next step
  cv_qwait: number = 0;               // steps before order change allowed
  cv_L: number = 0;                    // q + 1
  cv_h: number = 0;                    // current step size
  cv_hprime: number = 0;              // step size for next step
  cv_eta: number = 0;                  // h ratio: hprime / h
  cv_hscale: number = 0;              // step size at last Nordsieck rescale
  cv_tn: number = 0;                   // current internal time
  cv_tretlast: number = 0;            // last return time

  // --- Method coefficients ---
  cv_l: Float64Array = new Float64Array(0);       // polynomial coefficients [L_MAX]
  cv_tq: Float64Array = new Float64Array(0);      // error test quantities [NUM_TESTS+1]
  cv_tau: Float64Array = new Float64Array(0);      // step size history [L_MAX+1]

  // --- Implicit system ---
  cv_gamma: number = 0;               // gamma = h * rl1
  cv_gammap: number = 0;              // gamma at last setup
  cv_gamrat: number = 0;              // gamma / gammap
  cv_rl1: number = 0;                 // 1 / l[1]
  cv_crate: number = 0;               // NLS convergence rate

  // --- Solution vectors ---
  cv_y: Float64Array = new Float64Array(0);     // user's solution vector
  cv_acor: Float64Array = new Float64Array(0);  // accumulated correction
  cv_tempv: Float64Array = new Float64Array(0); // temporary vector
  cv_ftemp: Float64Array = new Float64Array(0); // temporary for f evaluations

  // --- Tolerances ---
  cv_reltol: number = 0;              // relative tolerance
  cv_Sabstol: number = 0;             // scalar absolute tolerance
  cv_Vabstol: Float64Array = new Float64Array(0); // vector absolute tolerance
  cv_itol: number = 0;                // tolerance type flag

  // --- Limits ---
  cv_qmax: number = 0;                // max method order
  cv_mxstep: number = 0;              // max internal steps
  cv_hmax_inv: number = 0;            // 1/hmax
  cv_hmin: number = 0;                // min step size
  cv_hin: number = 0;                 // initial step size (user)
  cv_nlscoef: number = 0;             // NLS convergence coefficient
  cv_mxhnil: number = 0;              // max t+h==t warnings

  // --- Linear solver ---
  cv_linit: ((mem: CvodeMem) => number) | null = null;
  cv_lsetup: ((mem: CvodeMem, convfail: number, ypred: Float64Array,
               fpred: Float64Array, jcurPtr: { value: boolean },
               tmp1: Float64Array, tmp2: Float64Array,
               tmp3: Float64Array) => number) | null = null;
  cv_lsolve: ((mem: CvodeMem, b: Float64Array, weight: Float64Array,
               ycur: Float64Array, fcur: Float64Array) => number) | null = null;
  cv_lfree: ((mem: CvodeMem) => void) | null = null;
  cv_lmem: CvLsMem | null = null;     // linear solver memory

  // --- Counters ---
  cv_nst: number = 0;                 // total steps
  cv_nfe: number = 0;                 // f evaluations
  cv_nni: number = 0;                 // nonlinear iterations
  cv_ncfn: number = 0;                // convergence failures
  cv_netf: number = 0;                // error test failures
  cv_nsetups: number = 0;             // linear solver setups

  // --- Rootfinding ---
  cv_gfun: CvodeRootFn | null = null; // root function
  cv_nrtfn: number = 0;               // number of root functions
  cv_glo: Float64Array = new Float64Array(0);   // g values at start of step
  cv_ghi: Float64Array = new Float64Array(0);   // g values at end of step
  cv_grout: Float64Array = new Float64Array(0); // g values at root
  cv_iroots: Int32Array = new Int32Array(0);     // root direction info
  cv_rootdir: Int32Array = new Int32Array(0);    // root direction constraints
  cv_gactive: Uint8Array = new Uint8Array(0);    // active root functions
  cv_tlo: number = 0;                 // t at low end of bracket
  cv_thi: number = 0;                 // t at high end
  cv_trout: number = 0;               // t at root
  cv_taskc: number = 0;               // root task counter
  cv_irfnd: number = 0;               // root found flag
  cv_nge: number = 0;                 // g evaluations

  // --- Flags ---
  cv_MallocDone: boolean = false;
  cv_SetupDone: boolean = false;
  cv_VabstolMallocDone: boolean = false;
  cv_jcur: boolean = false;           // is Jacobian current?
  cv_tstopset: boolean = false;
  cv_tstop: number = 0;

  // --- Problem size ---
  cv_N: number = 0;                   // number of equations
}
```

### 3.2 CvLsMem — Linear solver state

Maps from `CVLsMemRec` in `cvode_ls_impl.h`:

```typescript
export class CvLsMem {
  // Jacobian matrix (dense, row-major Float64Array[])
  A: Float64Array[] = [];              // N x N matrix
  pivots: Int32Array = new Int32Array(0);  // LU pivot indices
  savedJ: Float64Array[] = [];         // saved copy of Jacobian

  // User-supplied Jacobian (optional)
  jacFn: CvodeJacFn | null = null;

  // Counters
  nje: number = 0;                     // Jacobian evaluations
  nfeDQ: number = 0;                   // f evals for DQ Jacobian

  // State
  nstlj: number = 0;                  // nst at last Jacobian eval
  jacCurrent: boolean = false;         // is Jacobian current?
}
```

### 3.3 Function types

```typescript
// RHS function: y' = f(t, y)
export type CvodeRhsFn = (
  t: number, y: Float64Array, ydot: Float64Array, userData: any
) => number;

// Jacobian function: J = df/dy
export type CvodeJacFn = (
  t: number, y: Float64Array, fy: Float64Array,
  J: Float64Array[], userData: any
) => number;

// Root function: g(t, y)
export type CvodeRootFn = (
  t: number, y: Float64Array, gout: Float64Array, userData: any
) => number;
```

### 3.4 Constants and enums

```typescript
// Method types
export const CV_ADAMS = 1;
export const CV_BDF = 2;

// Task types
export const CV_NORMAL = 1;        // integrate to tout, interpolate
export const CV_ONE_STEP = 2;      // take one internal step

// Return values
export const CV_SUCCESS = 0;
export const CV_TSTOP_RETURN = 1;
export const CV_ROOT_RETURN = 2;
export const CV_TOO_MUCH_WORK = -1;
export const CV_TOO_MUCH_ACC = -2;
export const CV_ERR_FAILURE = -3;
export const CV_CONV_FAILURE = -4;
export const CV_RHSFUNC_FAIL = -8;

// Internal constants
export const L_MAX = 13;           // max order + 1 (Adams max = 12)
export const NUM_TESTS = 5;
export const HMIN_DEFAULT = 0.0;
export const HMAX_INV_DEFAULT = 0.0;
export const MXHNIL_DEFAULT = 10;
export const MXSTEP_DEFAULT = 500;
export const NLS_MAXCOR = 3;       // max corrector iterations
export const CRDOWN = 0.3;         // convergence rate damping
export const RDIV = 2.0;           // divergence threshold
export const DGMAX_LSETUP = 0.3;   // gamma change threshold for Jacobian redo
export const MSBP = 20;            // max steps between Jacobian setups
export const ADDON = 1.0e-6;       // eta floor addend
export const BIAS1 = 6.0;          // step size bias factors
export const BIAS2 = 6.0;
export const BIAS3 = 10.0;
export const ETAMXF = 0.2;         // max eta on error test fail
export const ETAMX1 = 10000.0;     // max eta on first step
export const ETAMX2 = 10.0;        // max eta on subsequent steps
export const ETAMX3 = 10.0;
export const ETAMIN = 0.1;         // min eta
export const THRESH = 1.5;         // eta threshold to skip rescale
export const SMALL_NST = 10;       // "small" step count
```

### 3.5 Indexing convention

**Use 0-based indexing** for CVODE, unlike LSODA which preserved 1-based.

LSODA inherited Fortran 1-based indexing through the C port and it was too
entangled to change. CVODE's C code uses **native 0-based indexing** (`NV_Ith_S`
macros use 0-based access, `N_VGetArrayPointer` returns 0-based arrays). This
means:

- All internal arrays are 0-based: `cv_zn[j][i]` where `i = 0..N-1`
- No `+1` allocation padding needed
- No 0↔1 conversion at the public API boundary
- The user-facing API is directly 0-based

This is simpler than the LSODA approach and avoids the conversion overhead.

---

## 4. Implementation Phases

### Phase 1: Foundation — Types, constants, and linear algebra

**Files:** `src/cvode/common.ts`, `src/cvode/dense_linalg.ts`

**Tasks:**
1. Define `CvodeMem`, `CvLsMem` classes with all fields
2. Define function types (`CvodeRhsFn`, `CvodeJacFn`, `CvodeRootFn`)
3. Define all constants and return codes
4. Implement or verify reuse of dense LU factorization (`dgefa`) and solve
   (`dgesl`) from existing `src/blas.ts`
   - If `blas.ts` uses 1-based indexing internally, create 0-based wrappers in
     `dense_linalg.ts`
   - Alternatively, write fresh 0-based `dgefa`/`dgesl` in `dense_linalg.ts`
5. Implement WRMS norm: `wrmsNorm(n, v, w)` — a core utility used everywhere

**Validation:** Unit tests for LU decomposition and WRMS norm.

### Phase 2: Core step algorithm

**Files:** `src/cvode/cvode.ts`, `src/cvode/cvode_bdf.ts`,
`src/cvode/cvode_adams.ts`, `src/cvode/cvode_hin.ts`

**Tasks:**
1. `cvodeCreate(lmm)` — allocate and initialize `CvodeMem`
2. `cvodeInit(mem, f, t0, y0)` — set RHS, allocate Nordsieck array and vectors
3. `cvHin(mem)` — compute initial step size `h0` (from `cvode_hin.ts`)
4. `cvInitialSetup(mem)` — validate inputs, compute `ewt`, call `cvHin`
5. `cvPredict(mem)` — advance Nordsieck array (Pascal triangle summation)
6. `cvSetBDF(mem)` — compute BDF coefficients `l[]`, `tq[]` (from `cvode_bdf.ts`)
7. `cvSetAdams(mem)` — compute Adams coefficients (from `cvode_adams.ts`)
8. `cvSet(mem)` — dispatch to `cvSetBDF` or `cvSetAdams`
9. `cvDoErrorTest(mem, ...)` — WRMS error test, step size reduction on failure
10. `cvCompleteStep(mem)` — update Nordsieck array with correction, update `tau`
11. `cvPrepareNextStep(mem)` — select new order `qprime` and step ratio `eta`
12. `cvAdjustParams(mem)` / `cvRescale(mem)` — rescale Nordsieck when `h` changes
13. `cvStep(mem)` — single internal step: predict → set → NLS → error test → complete
14. `cvHandleNFlag(mem, ...)` — handle nonlinear solver failure/retry

**Validation:** Integrate a simple non-stiff ODE (exponential decay) with Adams
method, then a stiff ODE (Robertson) with BDF. Compare against reference output
from `CVODE.md`.

### Phase 3: Nonlinear solver

**File:** `src/cvode/cvode_nls.ts`

**Tasks:**
1. `cvNlsResidual(mem, ycor, res)` — compute `G = rl1*zn[1] + ycor - gamma*f(tn, zn[0]+ycor)`
2. `cvNlsFPFunction(mem, ycor, res)` — fixed-point function for non-stiff
3. `cvNlsLSetup(mem, convfail)` — wrapper to call linear solver setup
4. `cvNlsLSolve(mem, b)` — wrapper to call linear solver solve
5. `cvNlsConvTest(mem, del, delp, tol)` — convergence test:
   `del * min(1, crate) / tol < 1`
6. `cvNls(mem)` — top-level: run Newton or fixed-point iteration (max
   `NLS_MAXCOR` = 3 iterations)

**Design note:** CVODE uses a pluggable `SUNNonlinearSolver` interface. For the
TS port, implement Newton and fixed-point iterations directly inside
`cvode_nls.ts` without the abstraction layer. Use a simple flag
(`mem.cv_lmm === CV_BDF` → Newton, `CV_ADAMS` → fixed-point) to dispatch.

**Validation:** Verify Newton convergence on stiff Robertson problem; verify
fixed-point convergence on non-stiff test.

### Phase 4: Linear solver

**File:** `src/cvode/cvode_ls.ts`

**Tasks:**
1. `cvodeSetLinearSolver(mem)` — initialize `CvLsMem`, set function pointers
2. `cvodeSetJacFn(mem, jac)` — set user-supplied Jacobian
3. `cvLsSetup(mem, convfail, ypred, fpred, ...)` — build `A = I - gamma*J`,
   factorize (LU)
   - If Jacobian is stale: evaluate Jacobian (user-supplied or DQ)
   - Scale: `A[i][j] = -gamma * J[i][j]`, then `A[i][i] += 1`
   - Factorize: call `dgefa`
4. `cvLsSolve(mem, b, weight, ycur, fcur)` — solve `A*x = b` via `dgesl`
5. `cvLsDQJac(mem, t, y, fy, J)` — difference-quotient Jacobian approximation
   - For each column j: perturb `y[j]`, evaluate `f`, compute `J[:,j]` by finite
     difference
6. `cvLsLinSys(mem, ...)` — default linear system function `A = I - gamma*J`

**Design note:** No SUNMatrix/SUNLinearSolver abstraction. The matrix is a plain
`Float64Array[]` (row-major). The LU factorize/solve are direct function calls.
This eliminates ~1,000 lines of SUNDIALS abstraction code.

**Validation:** Solve Robertson problem with both DQ Jacobian and user-supplied
Jacobian; verify matching statistics.

### Phase 5: Main integration loop

**File:** `src/cvode/cvode.ts` (extend)

**Tasks:**
1. `cvode(mem, tout, yout, itask)` — main entry point
   - `CV_NORMAL`: integrate to `tout`, interpolate final answer
   - `CV_ONE_STEP`: take one internal step, return
2. Handle `tstop` (stop time constraint)
3. `cvodeGetDky(mem, t, k, dky)` — dense output via Horner evaluation on
   Nordsieck array
4. `cvEwtSet(mem, ycur)` — compute error weight vector
5. Step loop: call `cvStep` repeatedly until `tout` (or `tstop`) is reached
6. Interpolation: use `cvodeGetDky` for output at exact `tout`

**Validation:** Full Robertson problem integration; match reference output from
`CVODE.md` (542 steps, 754 f evals, 11 Jacobian evals).

### Phase 6: Rootfinding

**File:** `src/cvode/cvode_root.ts`

**Tasks:**
1. `cvodeRootInit(mem, nrtfn, g)` — initialize rootfinding state
2. `cvRcheck1(mem)` — check for roots at initial time
3. `cvRcheck2(mem)` — check for root in current step (sign changes in `g`)
4. `cvRcheck3(mem)` — check for root after successful step
5. `cvRootfind(mem)` — Illinois/bisection root isolation within a step
6. Integrate with main loop: `cvode()` returns `CV_ROOT_RETURN` when root found

**Validation:** Robertson problem with root functions from `CVODE.md` example;
verify root times match reference (`t = 2.6391e-01` for y3 reaching 0.01,
`t = 2.0790e+07` for y1 dropping to 1e-4).

### Phase 7: Configuration and statistics

**File:** `src/cvode/cvode_io.ts`

**Tasks:**
1. Tolerance functions:
   - `cvodeSStolerances(mem, rtol, atol)` — scalar tolerances
   - `cvodeSVtolerances(mem, rtol, atol_vec)` — vector tolerances
2. Optional inputs:
   - `cvodeSetMaxNumSteps(mem, mxstep)`
   - `cvodeSetMaxOrd(mem, maxord)`
   - `cvodeSetMaxStep(mem, hmax)` / `cvodeSetMinStep(mem, hmin)`
   - `cvodeSetInitStep(mem, hin)`
   - `cvodeSetStopTime(mem, tstop)`
   - `cvodeSetUserData(mem, data)`
3. Statistics getters:
   - `cvodeGetNumSteps(mem)` → `cv_nst`
   - `cvodeGetNumRhsEvals(mem)` → `cv_nfe`
   - `cvodeGetNumLinSolvSetups(mem)` → `cv_nsetups`
   - `cvodeGetNumErrTestFails(mem)` → `cv_netf`
   - `cvodeGetNumNonlinSolvIters(mem)` → `cv_nni`
   - `cvodeGetNumNonlinSolvConvFails(mem)` → `cv_ncfn`
   - `cvodeGetNumJacEvals(mem)` → `cv_lmem.nje`
   - `cvodeGetNumGEvals(mem)` → `cv_nge`
   - `cvodeGetIntegratorStats(mem)` — return all stats as an object

**Validation:** Verify stats match reference output for Robertson problem.

### Phase 8: High-level class API

**File:** `src/cvode/cvode_class.ts`

**Tasks:**
1. Design `Cvode` class following `Lsoda` class pattern:

```typescript
export interface CvodeOptions {
  lmm?: 'adams' | 'bdf';           // method (default: 'bdf')
  rtol?: number;                     // relative tolerance
  atol?: number | Float64Array;     // absolute tolerance (scalar or vector)
  maxSteps?: number;                 // max internal steps
  maxOrder?: number;                 // max method order
  maxStep?: number;                  // max step size
  minStep?: number;                  // min step size
  initStep?: number;                // initial step size
  stopTime?: number;                 // stop time
  rootFn?: (t: number, y: Float64Array, gout: Float64Array) => void;
  nRootFns?: number;                 // number of root functions
  jacFn?: (t: number, y: Float64Array, fy: Float64Array, J: Float64Array[]) => void;
}

export interface CvodeSolveResult {
  t: number;                         // time reached
  y: Float64Array;                   // solution at t (0-based)
  flag: number;                      // CV_SUCCESS, CV_ROOT_RETURN, CV_TSTOP_RETURN
  rootsFound?: Int32Array;           // which roots were found (if flag == CV_ROOT_RETURN)
}

export interface CvodeStats {
  nSteps: number;
  nRhsEvals: number;
  nLinSolvSetups: number;
  nErrTestFails: number;
  nNonlinSolvIters: number;
  nNonlinSolvConvFails: number;
  nJacEvals: number;
  nGEvals: number;
  lastOrder: number;
  lastStep: number;
  currentTime: number;
}

export class Cvode {
  constructor(f: OdeFunction, neq: number, options?: CvodeOptions)
  solve(y: Float64Array, t: number, tout: number): CvodeSolveResult
  getDky(t: number, k: number): Float64Array   // dense output
  getStats(): CvodeStats
  reInit(t0: number, y0: Float64Array): void   // re-initialize for new IC
  free(): void
}
```

2. Constructor handles all setup: create solver, set tolerances, attach linear
   solver, set Jacobian, initialize rootfinding
3. `solve()` calls `cvode()` internally, returns clean result object
4. No 0↔1 conversion needed (CVODE is natively 0-based)

**Validation:** Run Robertson problem through high-level API; verify identical
results to low-level API.

### Phase 9: Barrel exports and integration

**File:** `src/cvode/index.ts`

**Tasks:**
1. Export public API:
   ```typescript
   // High-level class API
   export { Cvode } from './cvode_class';
   export type { CvodeOptions, CvodeSolveResult, CvodeStats } from './cvode_class';

   // Low-level procedural API
   export { cvodeCreate, cvodeInit, cvode, cvodeGetDky } from './cvode';
   export { cvodeSStolerances, cvodeSVtolerances,
            cvodeSetMaxNumSteps, cvodeSetMaxOrd, /* ... */
            cvodeGetIntegratorStats } from './cvode_io';
   export { cvodeSetLinearSolver, cvodeSetJacFn } from './cvode_ls';
   export { cvodeRootInit } from './cvode_root';

   // Types and constants
   export { CvodeMem, CvodeRhsFn, CvodeJacFn, CvodeRootFn,
            CV_ADAMS, CV_BDF, CV_NORMAL, CV_ONE_STEP,
            CV_SUCCESS, CV_ROOT_RETURN, CV_TSTOP_RETURN } from './common';
   ```

2. Update root `src/index.ts` to re-export CVODE module:
   ```typescript
   export * from './cvode';
   ```

### Phase 10: Testing

Tests are split across two files, mirroring the existing LSODA test structure.

#### 10.1 `test/cvode.test.ts` — Correctness tests (mirrors `test/test.test.ts`)

**Structure:** `describe('CVODE solver', () => { ... })` top-level block.

**Robertson chemical kinetics** (stiff, BDF):
- Integrate from `t=0` to `t=4e10` with output at `t = 0.4, 4, 40, ..., 4e10`
  (12 steps, each `tout *= 10`)
- `y(0) = [1, 0, 0]`, `rtol = [1e-4, 1e-4, 1e-4]`, `atol = [1e-6, 1e-10, 1e-6]`
- Verify each output point against reference values (same as `test.test.ts`):
  ```
  t = 4.0000e-01   y ≈ [9.851712e-01, 3.386380e-05, 1.479493e-02]
  t = 4.0000e+00   y ≈ [9.055333e-01, 2.240655e-05, 9.444430e-02]
  ...
  t = 4.0000e+10   y ≈ [1.431100e-08, 5.724404e-14, 1.000000e+00]
  ```
- Relative error check: `relError <= 5e-3` per component at each output point

**Non-stiff problems** (`describe('non-stiff problems', ...)`):

1. **Non-stiff 1D** — `y' = 4*exp(0.8*t) - 0.5*y`, `y(0) = 2`,
   `t ∈ [0, 4]`, step `0.01`, `tol = 1e-5`
   Exact: `y(t) = (exp(0.8t) - exp(-0.5t)) * 4/1.3 + 2*exp(-0.5t)`

2. **Non-stiff 2D** — `y1' = y1+y2`, `y2' = y2-y1`, `y(0) = [1, 1]`,
   `t ∈ [0, 4]`, step `0.01`, `tol = 1e-9`
   Exact: `y1 = exp(t)*(cos(t)+sin(t))`, `y2 = exp(t)*(cos(t)-sin(t))`

3. **Non-stiff 3D** — `y1' = 5*y1+2*y2+sin(t)`, `y2' = -4*y1-y2+exp(2t)`,
   `y3' = 5*t^4-3*t^2+2*t`, `y(0) = [0.3, -0.8, 0]`,
   `t ∈ [0, 2]`, step `0.001`, `tol = 1e-8`
   Exact solution via matrix exponential + particular integral

**Stiff problems** (`describe('stiff problems', ...)`):

1. **Stiff 1D** — `y' = -1000*y + 3000 - 2000*exp(-t)`, `y(0) = 0`,
   `t ∈ [0, 4]`, step `0.01`, `tol = 5e-7`
   Exact: `y(t) = 3 - 0.998*exp(-1000t) - 2.002*exp(-t)`

2. **Stiff 2D** — `y1' = -5*y1+3*y2`, `y2' = 100*y1-301*y2`,
   `y(0) = [52.29, 83.82]`, `t ∈ [0, 4]`, step `0.01`, `tol = 5e-7`
   Exact: two-exponential solution

3. **Stiff 3D** — `y1' = -100*y1+100*y2+y3`, `y2' = y3`, `y3' = -y2`,
   `y(0) = [2, 1, 0]`, `t ∈ [0, 4]`, step `0.01`, `tol = 1e-8`
   Exact: `y1 = exp(-100t)+cos(t)`, `y2 = cos(t)`, `y3 = -sin(t)`

**All non-stiff and stiff tests:** integrate to final time, compare against
exact solution with combined tolerance: `|computed - exact| < checkTol * |exact| + checkAbsTol`

**Rootfinding** (`describe('rootfinding', ...)`):
- Robertson problem with root functions `g1 = y1 - 1e-4`, `g2 = y3 - 0.01`
- Verify root at `t ≈ 2.6391e-01` (y3 reaches 0.01)
- Verify root at `t ≈ 2.0790e+07` (y1 drops to 1e-4)
- Verify `CV_ROOT_RETURN` flag and `rootsFound` array

**Dense output** (`describe('dense output', ...)`):

1. **Consistency** — solve Robertson with dense output, query at same output
   points, verify match against direct solver results (`relError <= 5e-3`)

2. **Monotonicity** — query 1000 uniform points across `[tMin, tMax]`, verify
   all values are finite

3. **Boundary** — query at `tMin` and `tMax`, verify finite results

4. **Dense output vs exact** — for each non-stiff and stiff problem:
   solve to final time with dense output, query at 100 uniform points,
   compare against exact solution with combined tolerance

**Error handling** (`describe('error handling', ...)`):
- Invalid tolerances (negative rtol, zero atol)
- `CV_TOO_MUCH_WORK` (mxstep exceeded)
- `CV_RHSFUNC_FAIL` (user function returns error)

**High-level API** (`describe('Cvode class', ...)`):
- Same Robertson, non-stiff, and stiff problems through `Cvode` class
- Verify `getStats()` returns valid statistics
- Verify `reInit()` for solving with new initial conditions

#### 10.2 `test/cvode.perf.test.ts` — Benchmark tests (mirrors `test/perf.test.ts`)

**Structure:** Each problem is a separate `describe` block. Same `ODEs` interface,
same `solveBenchmarkAndCheck` helper pattern as `perf.test.ts`. Each test has
`60000` ms timeout.

**Problems (all stiff, using BDF):**

1. **Robertson** — 3 equations, `t ∈ [0, 10e11]`, step `2.5e6`, `tol = 1e-7`
   Reference: `[0.2083340149701255e-7, 0.8333360770334713e-13, 0.9999999791665050]`
   Options: `mxstep = 50000`, `atol = [1e-8, 1e-14, 1e-8]`, `rtol = [1e-7, 1e-7, 1e-7]`

2. **HIRES** (High Irradiance Responses) — 8 equations, `t ∈ [0, 321.8122]`,
   step `0.01`, `tol = 1e-10`
   Reference: `[0.7371e-3, 0.1442e-3, 0.5889e-4, 0.1176e-2, 0.2386e-2, 0.6239e-2, 0.2850e-2, 0.2850e-2]`

3. **Van der Pol** (mu=1000) — 2 equations, `t ∈ [0, 2000]`, step `0.1`,
   `tol = 1e-12`
   Reference: `[1.706167732170469, -0.8928097010248125e-3]`
   Options: `mxstep = 50000`

4. **OREGO** (Oregonator) — 3 equations, `t ∈ [0, 360]`, step `0.01`,
   `tol = 1e-8`
   Reference: `[1.000814870318523, 1228.178521549917, 132.0554942846706]`

5. **E5** (Chemical pyrolysis) — 4 equations, `t ∈ [0, 1e13]`, step `2.5e8`,
   `tol = 1e-6`, extremely stiff (rate constants span 20 orders)
   Reference: `[0.1153e-290, 0.8868e-22, 0.8855e-22, 0]`
   Options: `mxstep = 50000`, `warmupTout = 1e-5`

6. **Pollution** (air pollution model) — 20 equations, `t ∈ [0, 60]`,
   step `0.002`, `tol = 1e-6`, 25 reaction rate constants
   Reference: 20-component vector (see `perf.test.ts`)

**Each benchmark test:**
- Integrates from `start` to `finish` in steps of `step`
- Checks solver state > 0 at each step (no failure)
- Verifies final solution against reference point with problem-specific
  `relTol` and `absTol`

---

## 5. Transformation Patterns (CVODE-specific)

### 5.1 SUNDIALS abstractions → Direct calls

CVODE's C code goes through multiple abstraction layers:

```c
// C: SUNMatrix + SUNLinearSolver + N_Vector abstractions
N_Vector y = N_VNew_Serial(N, ctx);
NV_Ith_S(y, i) = value;                       // macro for y->data[i]
SUNMatrix A = SUNDenseMatrix(N, N, ctx);
SM_ELEMENT_D(A, i, j) = value;                // macro for A->data[j*N+i]
SUNLinearSolver LS = SUNLinSol_Dense(y, A, ctx);
retval = SUNLinSolSetup(LS, A);
retval = SUNLinSolSolve(LS, A, x, b, tol);
```

```typescript
// TypeScript: plain arrays, direct function calls
const y = new Float64Array(N);
y[i] = value;
const A: Float64Array[] = Array.from({length: N}, () => new Float64Array(N));
A[i][j] = value;                               // row-major
const pivots = new Int32Array(N);
dgefa(A, N, pivots);                           // direct LU call
dgesl(A, N, pivots, b, 0);                     // direct solve call
```

**Rule:** Eliminate all SUNDIALS wrapper types. `N_Vector` → `Float64Array`.
`SUNMatrix` → `Float64Array[]`. `SUNLinearSolver` → direct `dgefa`/`dgesl`
calls. `SUNNonlinearSolver` → inline Newton/fixed-point code.

### 5.2 N_Vector operations → Typed array operations

| C (N_Vector) | TypeScript |
|--------------|------------|
| `N_VNew_Serial(N, ctx)` | `new Float64Array(N)` |
| `NV_Ith_S(v, i)` | `v[i]` |
| `N_VGetArrayPointer(v)` | `v` (already a Float64Array) |
| `N_VLinearSum(a, x, b, y, z)` | `for (i) z[i] = a*x[i] + b*y[i]` or helper |
| `N_VScale(c, x, z)` | `for (i) z[i] = c * x[i]` or helper |
| `N_VConst(c, z)` | `z.fill(c)` |
| `N_VClone(v)` | `new Float64Array(v)` (copy) or `new Float64Array(v.length)` |
| `N_VDestroy(v)` | *(no-op, GC)* |
| `N_VAbs(x, z)` | `for (i) z[i] = Math.abs(x[i])` |
| `N_VWrmsNorm(x, w)` | `wrmsNorm(N, x, w)` (helper function) |

Create small inline helpers for `N_VLinearSum`, `N_VScale` etc. in `common.ts`,
or inline them at call sites when used only once.

### 5.3 CVodeMem access pattern

```c
// C: void* cast + field access
CVodeMem cv_mem = (CVodeMem)cvode_mem;
cv_mem->cv_tn += cv_mem->cv_h;
```

```typescript
// TypeScript: typed class, direct access
const m = mem;  // CvodeMem, already typed
m.cv_tn += m.cv_h;
```

Unlike LSODA's `_C(x)` macro pattern, CVODE accesses fields directly on the
`CVodeMem` struct. The TS translation is straightforward: `cv_mem->cv_field` →
`mem.cv_field`.

Consider dropping the `cv_` prefix in TypeScript since the class name already
provides context:

```typescript
// Option A: keep prefix (closer to C, easier to cross-reference)
mem.cv_tn += mem.cv_h;

// Option B: drop prefix (more idiomatic TS)
mem.tn += mem.h;
```

**Decision:** Keep `cv_` prefix for the low-level `CvodeMem` class fields. This
makes cross-referencing with the C source trivial during translation and
debugging. The high-level `Cvode` class uses clean names (`getStats()`,
`solve()`, etc.).

### 5.4 Error handling

CVODE uses integer return codes extensively:

```c
// C
retval = cv_mem->cv_f(tn, y, ftemp, cv_mem->cv_user_data);
if (retval < 0) { cvProcessError(..., CV_RHSFUNC_FAIL, ...); return CV_RHSFUNC_FAIL; }
if (retval > 0) { /* recoverable error */ }
```

```typescript
// TypeScript
const retval = mem.cv_f!(mem.cv_tn, y, ftemp, mem.cv_user_data);
if (retval < 0) return CV_RHSFUNC_FAIL;
if (retval > 0) { /* recoverable error */ }
```

**Rule:** Drop `cvProcessError` (which formats error messages and calls user
error handler). Instead, set an error string on the mem object for debugging:

```typescript
if (retval < 0) {
  mem.cv_error = `RHS function failed at t = ${mem.cv_tn}`;
  return CV_RHSFUNC_FAIL;
}
```

### 5.5 `goto` elimination

CVODE's `cvode.c` uses `goto` for error cleanup paths. TypeScript has no `goto`.
Replace with early returns, labeled `break`, or restructured control flow:

```c
// C: goto for cleanup
retval = cvStep(cv_mem);
if (retval != CV_SUCCESS) goto error_return;
// ...
error_return:
  // cleanup
  return retval;
```

```typescript
// TypeScript: early return
const retval = cvStep(mem);
if (retval !== CV_SUCCESS) {
  // cleanup
  return retval;
}
```

For complex multi-level cleanup, use a try/finally or break from a labeled block.

---

## 6. Key Algorithmic Details to Preserve

### 6.1 Nordsieck array operations

The Nordsieck array `zn[j]` for `j = 0..q` stores scaled derivatives. All
operations must be translated exactly:

- **Predict** (Pascal triangle):
  ```
  for k = 1..q:
    for j = q down to k:
      zn[j-1][i] += zn[j][i]   (for all i = 0..N-1)
  ```

- **Correct** (after NLS converges):
  ```
  for j = 0..q:
    zn[j][i] += l[j] * acor[i]   (for all i = 0..N-1)
  ```

- **Rescale** (when h changes by factor eta):
  ```
  for j = 1..q:
    scale = eta^j
    zn[j][i] *= scale   (for all i = 0..N-1)
  ```

### 6.2 BDF coefficient computation

The BDF method coefficients `l[]` and `tq[]` depend on the step size history
`tau[]`. The computation in `cvSetBDF` involves building the characteristic
polynomials from past step ratios. This is numerically sensitive and must be
translated exactly.

### 6.3 Adams coefficient computation

Adams-Moulton coefficients require computing `l[]` and `tq[]` from the step
size history. The `cvSetAdams` function builds these using a more complex
algorithm involving polynomial products. The static data arrays
(`cv_AdamsStart`, etc.) should be translated as module-level `const` arrays.

### 6.4 Initial step size selection

`cvHin` estimates the initial step size from the RHS and its rate of change.
This involves:
1. Computing `||y'||_wrms`
2. Estimating `||y''||_wrms` by finite difference
3. Combining: `h = (1/max(ydd, yd^2))^(1/(q+1))`

### 6.5 Rootfinding algorithm

The rootfinding in CVODE uses a modified Illinois method (regula falsi variant)
with bisection fallback. The algorithm:
1. Detect sign changes in `g(t,y)` between step boundaries
2. Isolate root to within step tolerance using Illinois/bisection
3. Handle multiple simultaneous roots
4. Handle root direction constraints

---

## 7. Dependency Graph (Build Order)

```
common.ts                  ← Phase 1 (no dependencies)
    │
    ├── dense_linalg.ts    ← Phase 1 (depends on common for types)
    │       │
    ├── cvode_bdf.ts       ← Phase 2 (depends on common)
    ├── cvode_adams.ts     ← Phase 2 (depends on common)
    ├── cvode_hin.ts       ← Phase 2 (depends on common)
    │       │
    ├── cvode_nls.ts       ← Phase 3 (depends on common, cvode_ls)
    │       │
    ├── cvode_ls.ts        ← Phase 4 (depends on common, dense_linalg)
    │       │
    ├── cvode_io.ts        ← Phase 7 (depends on common)
    │       │
    └── cvode.ts           ← Phase 2+5 (depends on all above)
            │
    ├── cvode_root.ts      ← Phase 6 (depends on common, cvode)
    │       │
    ├── cvode_class.ts     ← Phase 8 (depends on all above)
    │       │
    └── index.ts           ← Phase 9 (depends on all above)
```

---

## 8. Testing Strategy

**Test files:**
- `test/cvode.test.ts` — correctness tests (see Phase 10.1 for full spec)
- `test/cvode.perf.test.ts` — benchmark tests (see Phase 10.2 for full spec)

### 8.1 Bottom-up unit tests (in `test/cvode.test.ts`)

Each phase produces testable units:

| Phase | Unit tests |
|-------|------------|
| 1 | LU factorize/solve on known matrices; WRMS norm |
| 2 | `cvPredict` on known Nordsieck array; BDF/Adams coefficients vs tabulated values |
| 3 | Newton convergence on simple system; fixed-point convergence |
| 4 | DQ Jacobian vs analytical Jacobian; `I - gamma*J` assembly |
| 5 | Full integration: non-stiff 1D/2D/3D, stiff 1D/2D/3D, Robertson |
| 6 | Rootfinding on Robertson with known root times |
| 7 | Statistics match reference |
| 8 | High-level `Cvode` class produces same results as low-level API |

### 8.2 Benchmark tests (in `test/cvode.perf.test.ts`)

Six stiff benchmark problems with published reference solutions:
Robertson, HIRES, van der Pol (mu=1000), OREGO, E5, Pollution.
Same problems, reference points, tolerances, and helper infrastructure as the
existing `test/perf.test.ts` — but importing from `src/cvode/` instead of
`src/`.

### 8.3 Reference test: Robertson problem

The primary validation target. From `CVODE.md`:

```
Input:
  y' = [-0.04*y1 + 1e4*y2*y3,
         0.04*y1 - 1e4*y2*y3 - 3e7*y2^2,
         3e7*y2^2]
  y(0) = [1, 0, 0]
  rtol = 1e-4, atol = [1e-8, 1e-14, 1e-6]

Expected output (selected points):
  t = 0.4:     y ≈ [9.8516e-01, 3.3862e-05, 1.4802e-02]
  t = 4e10:    y ≈ [6.9345e-08, 2.7738e-13, 1.0000e+00]

Expected statistics:
  Steps: 542, f evals: 754, Jacobian evals: 11
  Roots found: t ≈ 2.6391e-01 (y3=0.01), t ≈ 2.0790e+07 (y1=1e-4)
```

### 8.4 Comparison with LSODA

Since both LSODA and CVODE solve the same class of problems, cross-validate:
- Solve Robertson with both solvers (BDF mode)
- Solutions should match to within tolerances
- CVODE may take fewer/more steps due to different coefficient strategies

---

## 9. Summary Checklist

When transforming each CVODE C function to TypeScript:

1. **File:** Map function to appropriate `.ts` module per Section 2.
2. **Imports:** Replace `#include` with `import { ... } from '...'`.
3. **Signature:** `static int cvFunc(CVodeMem mem, ...)` →
   `export function cvFunc(mem: CvodeMem, ...): number`.
4. **SUNDIALS types:** `N_Vector` → `Float64Array`, `SUNMatrix` →
   `Float64Array[]`, `sunrealtype` → `number`.
5. **N_Vector ops:** `NV_Ith_S(v,i)` → `v[i]`, `N_VLinearSum` → loop,
   `N_VWrmsNorm` → `wrmsNorm()`.
6. **Memory:** `N_VNew_Serial(N)` → `new Float64Array(N)`, `N_VDestroy` →
   *(nothing)*, `SUNMatDestroy` → *(nothing)*.
7. **Access:** `cv_mem->cv_field` → `mem.cv_field`.
8. **Vars:** Declare at point of use; `const` for immutables, `let` for mutated.
9. **Math:** `SUNRabs` → `Math.abs`, `SUNRsqrt` → `Math.sqrt`,
   `SUNRpowerI(x,n)` → `Math.pow(x,n)`.
10. **Macros/constants:** `#define` → `export const`.
11. **Errors:** `cvProcessError(...)` → set `mem.cv_error` string + return code.
12. **`goto`:** Restructure to early returns or labeled breaks.
13. **Context:** Drop `SUNContext` parameter from all functions.
14. **Indexing:** Preserve 0-based (native CVODE convention).
