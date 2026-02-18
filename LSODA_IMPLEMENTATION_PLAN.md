# LSODA ODE Solver — Implementation Plan

## Goal

Implement the LSODA solver (automatic switching between Adams and BDF) fitting the existing `diff-grok` library conventions. The solver switches between Adams-Bashforth-Moulton methods (orders 1–12) for non-stiff regions and BDF methods (orders 1–5) for stiff regions, using the Nordsieck array representation internally.

---

## Library Conventions (must follow)

These patterns are extracted from `mrt-method.ts` and must be matched exactly:

```typescript
// Signature
export function lsoda(odes: ODEs, callback?: Callback): Float64Array[];

// Imports
import {ODEs, max, abs, SAFETY, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, EPS, jacobian, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {luDecomp, luSolve, solve1d2d} from './lin-alg-tools';
```

**Hard constraints:**
- All buffers are preallocated `Float64Array` — zero allocations inside the main loop.
- `Func` writes into a preallocated output array: `f(t, y, output)`.
- Use `jacobian()` from `solver-defs` for finite-difference Jacobian.
- Use `luDecomp` / `luSolve` from `lin-alg-tools` for dim > 2, `solve1d2d` for dim ≤ 2.
- Error norm: `errmax = max_i |yErr[i] / yScale[i]| / tolerance`.
- Step adaptation: `SAFETY`, `REDUCE_COEF`, `GROW_COEF`, `ERR_CONTR` constants.
- Output on a uniform grid via linear interpolation.
- Callback: `callback.onIterationStart()` every outer iteration, `callback.onComputationsCompleted()` at end.
- On failure: `throw new Error(ERROR_MSG.LSODA_FAILS)`.
- Return: `Float64Array[]` where `[0]` = t, `[1..dim]` = solution components.

---

## Architecture Overview

```
lsoda-method.ts          — main solver function (single file, follows MRT pattern)
solver-defs.ts           — existing, no changes (already has ERROR_MSG.LSODA_FAILS)
lin-alg-tools.ts         — existing, no changes
callbacks/callback-base  — existing, no changes
```

Everything goes into a single file `lsoda-method.ts`, just like `mrt-method.ts`.

---

## Phase 1: Constants and Coefficients

Define all constants at the top of the file, outside the function (same pattern as MRT's `D`, `E32`).

### 1.1 Adams Coefficients

Precompute correction vectors `l[0..q]` in Nordsieck form for Adams-Moulton corrector, orders 1–12.

```typescript
// ADAMS_L[q] = correction vector of length q+1 for Adams order q
// ADAMS_ERR_COEFF[q] = error constant C(q) for Adams order q
const ADAMS_MAX_ORDER = 12;
const ADAMS_L: number[][] = [/* 12 arrays */];
const ADAMS_ERR_COEFF: number[] = [/* 12 values */];
```

These are the `l` vectors from the Nordsieck formulation. They determine how the corrector updates each component of the Nordsieck array. Derive them from the standard Adams-Moulton generating polynomials or take from the LSODE Fortran source (subroutine `cfode`).

### 1.2 BDF Coefficients

Same structure for BDF, orders 1–5:

```typescript
const BDF_MAX_ORDER = 5;
const BDF_L: number[][] = [/* 5 arrays */];
const BDF_ERR_COEFF: number[] = [/* 5 values */];
const BDF_GAMMA: number[] = [/* β_s for each order: 1, 2/3, 6/11, 12/25, 60/137 */];
```

### 1.3 Solver Tuning Constants

```typescript
// Newton iteration
const MAX_NEWTON_ITER = 4;
const NEWTON_CONV_TOL = 0.33;

// Step size adaptation exponents (computed per order, but define helpers)
// For order q: PSHRNK = -1/(q+1), PSGROW = -1/(q+1) with safety

// Stiffness detection
const STIFF_DETECT_THRESHOLD = 5;
const NONSTIFF_TEST_INTERVAL = 30;

// Order change constraints
const MIN_STEPS_BEFORE_ORDER_UP = (q: number) => q + 1;

// Max step growth per step
const MAX_STEP_GROW = 10.0;
const MIN_STEP_SHRINK = 0.2;

// Jacobian refresh
const MAX_JACOBIAN_AGE = 20;
```

### 1.4 Testing Checkpoint

- Verify `BDF_L` and `BDF_GAMMA` reproduce the known BDF formulas for orders 1–5.
  For example, BDF2: `y_n = (4/3)y_{n-1} - (1/3)y_{n-2} + (2/3)h·f(t_n, y_n)`.
- Verify `ADAMS_L` vectors for orders 1–4 against textbook Adams-Moulton coefficients.

---

## Phase 2: Nordsieck Array Operations

Implement as **pure functions** at the bottom of the file (or as closures if they need solver buffers). Since all arrays are preallocated, these functions operate in-place.

### 2.1 Nordsieck Storage

The Nordsieck array is stored as a flat `Float64Array` of size `(maxOrder + 1) * dim`, with row-major layout:

```typescript
// nordsieck[i * dim + j] = i-th scaled derivative of j-th component
// i = 0: y_n
// i = 1: h * y'_n
// i = 2: h^2 * y''_n / 2!
// ...
// i = q: h^q * y^(q)_n / q!
```

Preallocate for maximum possible order:

```typescript
const maxNordsieckRows = ADAMS_MAX_ORDER + 1; // 13
const nordsieck = new Float64Array(maxNordsieckRows * dim);
const nordsieckPred = new Float64Array(maxNordsieckRows * dim);
```

### 2.2 Predict — Apply Pascal Shift

```typescript
function nordsieckPredict(z: Float64Array, zPred: Float64Array, q: number, dim: number): void
```

Compute `zPred = P · z` where `P[i][j] = C(j, i)` (binomial coefficient), for `0 ≤ i ≤ j ≤ q`.

This is a Taylor shift: predicts all Nordsieck components at `t + h` from values at `t`. Operate in-place into `zPred`. Use the recurrence for binomial coefficients (no factorials).

Implementation detail: iterate `i` from 0 to `q`, for each `i` sum `C(j, i) * z[j]` for `j = i..q`. The binomial coefficients `C(j, i)` can be precomputed in a small array at the start of the function.

### 2.3 Correct — Apply Correction Vector

```typescript
function nordsieckCorrect(
  zPred: Float64Array, z: Float64Array, l: number[], delta: Float64Array,
  q: number, dim: number
): void
```

After the corrector converges with correction `delta = y_corrected - y_predicted`:

```
z[i][j] = zPred[i][j] + l[i] * delta[j],   for i = 0..q, j = 0..dim-1
```

### 2.4 Rescale — After Step Size Change

```typescript
function nordsieckRescale(z: Float64Array, eta: number, q: number, dim: number): void
```

When step size changes by factor `eta = h_new / h_old`:

```
z[i][j] *= eta^i,   for i = 0..q
```

Precompute `eta^i` in a small local array of length `q + 1`.

### 2.5 Initialize Nordsieck Array

```typescript
function nordsieckInit(
  y0: Float64Array, f0: Float64Array, h: number, z: Float64Array, dim: number
): void
```

Set order 1: `z[0] = y0`, `z[1] = h * f0`. Zero out higher rows.

### 2.6 Testing Checkpoint

- Initialize Nordsieck for `y' = -y` at `y(0) = 1`, `h = 0.1`.
- Predict → should give Taylor approximation: `y_pred ≈ 1 + 0.1·(-1) = 0.9`.
- Verify rescale: after `eta = 2`, `z[1]` should double.

---

## Phase 3: Single Step — Predict-Correct Cycle

### 3.1 Predictor

1. Call `nordsieckPredict(nordsieck, nordsieckPred, q, dim)`.
2. Extract predicted y: `yPred[j] = nordsieckPred[0 * dim + j]` for `j = 0..dim-1`.

### 3.2 Corrector

**Adams mode (non-stiff):** functional iteration (no Jacobian needed).

Use **PECE mode** (Predict-Evaluate-Correct-Evaluate) — a single correction pass, which is standard in LSODE/LSODA:

1. Evaluate `f(t + h, yPred)` → `fNewton`.
2. Compute delta: `delta[j] = h * fNewton[j] - nordsieckPred[1 * dim + j]`.
3. Corrected y: `yTemp[j] = nordsieckPred[j] + l[0] * delta[j]`.
4. (Optionally: second evaluation at corrected y for PECE — evaluate f at corrected y, recompute delta.)

**BDF mode (stiff):** modified Newton iteration (needs Jacobian).

1. Form iteration matrix `W = I - l[0] * h * J` where `l[0]` absorbs the BDF gamma.
2. LU-decompose `W` (reuse if Jacobian is still valid).
3. Newton loop (max `MAX_NEWTON_ITER` iterations):
   a. Evaluate `f(t + h, y^(k))` → `fNewton`.
   b. Compute correction delta: `delta[j] = h * fNewton[j] - nordsieckPred[1 * dim + j]`.
   c. Compute residual: `G[j] = -(y^(k)[j] - nordsieckPred[j] - l[0] * delta[j])`.
   d. Solve `W · Δe = G` using LU (or `solve1d2d` for dim ≤ 2).
   e. Update: `y^(k+1)[j] += Δe[j]`.
   f. Convergence check: `max_j |Δe[j] / yScale[j]| < NEWTON_CONV_TOL`.

This maps directly to the existing library pattern:
- `W = I - l[0]*h*J` → same as MRT's `W = I - h*d*J`.
- `luDecomp(W, L, U, dim)` for dim > 2, `solve1d2d` for dim ≤ 2.
- `jacobian(t, y, f, EPS, f0Buf, f1Buf, Jac)` from `solver-defs`.

### 3.3 Error Estimation

After the corrector converges with correction `delta`:

```typescript
// errCoeff = method error constant for current order q
const errCoeff = isStiff ? BDF_ERR_COEFF[q - 1] : ADAMS_ERR_COEFF[q - 1];

for (let i = 0; i < dim; ++i)
  yErr[i] = errCoeff * delta[i];

errmax = 0;
for (let i = 0; i < dim; ++i)
  errmax = max(errmax, abs(yErr[i] / yScale[i]));
errmax /= tolerance;
```

Accept if `errmax ≤ 1`, reject otherwise. Same pattern as MRT.

### 3.4 Nordsieck Update on Accept

Call `nordsieckCorrect(nordsieckPred, nordsieck, l, delta, q, dim)` to update all Nordsieck components.

### 3.5 Testing Checkpoint

- Single Adams step (order 1) on `y' = -y`, verify against implicit midpoint.
- Single BDF step (order 1 = implicit Euler) on `y' = -y`, verify Newton converges in 1–2 iterations.
- Check that error estimate is O(h²) for order 1.

---

## Phase 4: Step Size and Order Selection

### 4.1 Step Size from Error Norm

For the accepted step with order `q` and error norm `errmax`:

```typescript
const psgrow = -1.0 / (q + 1);
let eta: number;
if (errmax > ERR_CONTR)
  eta = SAFETY * errmax ** psgrow;
else
  eta = GROW_COEF;
eta = Math.min(eta, MAX_STEP_GROW);
```

For rejected steps:

```typescript
const pshrnk = -1.0 / (q + 1);
let eta = SAFETY * errmax ** pshrnk;
eta = Math.max(eta, MIN_STEP_SHRINK);
```

### 4.2 Order Selection

After a successful step, evaluate three candidates. The error estimates come from existing Nordsieck components (cheap):

**Error for order q (current):** already computed as `errmax`.

**Error for order q−1** (if q > 1):

```typescript
// The q-th Nordsieck component (after correction) estimates the order-q-1 error
let errDown = 0;
for (let j = 0; j < dim; ++j) {
  const val = nordsieck[q * dim + j];
  errDown = max(errDown, abs(val / yScale[j]));
}
errDown /= tolerance;
const etaDown = (errDown > TINY) ? SAFETY * errDown ** (-1.0 / q) : GROW_COEF;
```

**Error for order q+1** (if `q < maxOrder` and enough steps at current order):

```typescript
// Requires storing the previous step's delta to estimate the (q+1)-th derivative
// errUp ≈ ‖delta_current - delta_previous‖ * C(q+1) / (h ratio factor)
let errUp = 0;
for (let j = 0; j < dim; ++j) {
  const val = (delta[j] - deltaPrev[j]) * ADAMS_ERR_COEFF[q]; // or BDF
  errUp = max(errUp, abs(val / yScale[j]));
}
errUp /= tolerance;
const etaUp = (errUp > TINY) ? SAFETY * errUp ** (-1.0 / (q + 2)) : GROW_COEF;
```

**Decision:**

```typescript
let qNew = q;
let etaBest = etaSame;

if (q > 1 && etaDown > etaBest) { qNew = q - 1; etaBest = etaDown; }
if (q < maxOrder && stepsAtCurrentOrder >= q + 1 && !wasRejected && etaUp > etaBest) {
  qNew = q + 1; etaBest = etaUp;
}
```

Apply the winning `eta` to compute `h_new = h * etaBest`. Then call `nordsieckRescale` if step size changes, and adjust the Nordsieck array if order changes (add or remove the last row).

### 4.3 State Variables for Order Tracking

```typescript
let q = 1;                    // current order
let stepsAtCurrentOrder = 0;  // steps since last order change
const deltaPrev = new Float64Array(dim);  // previous correction for q+1 estimate
let prevDeltaValid = false;
```

### 4.4 Testing Checkpoint

- On `y' = -y`: verify order ramps from 1 up to ~5–6 within 10–15 steps.
- On `y' = cos(t)` (smooth, non-stiff): verify Adams reaches high order (8+).
- On `y' = -1000y` (stiff): verify BDF stays at order 1–2.

---

## Phase 5: Stiffness Detection and Method Switching

### 5.1 Adams → BDF (Stiffness Detected)

Monitor during Adams mode. Two complementary signals:

**Signal A — Repeated step rejections:**
If `consecutiveAdamsRejects >= STIFF_DETECT_THRESHOLD`, switch to BDF.

**Signal B — Newton-like convergence rate check (optional, stronger):**
After each Adams step, estimate the dominant eigenvalue magnitude via the ratio of successive corrections or via `‖f(t+h, y) - f(t, y)‖ / ‖y_new - y‖`. If `h * |λ_max|` exceeds the Adams stability boundary for the current order, increment a stiffness counter.

Start with Signal A (simpler). Add Signal B later if needed.

### 5.2 BDF → Adams (Non-Stiffness Detected)

Every `NONSTIFF_TEST_INTERVAL` accepted BDF steps, perform a trial explicit step (RKF45). If the trial step's error is within tolerance at the current step size, the problem is no longer stiff — switch back to Adams.

Reuse the RKF45 Butcher tableau constants (same as in the bootstrap/trial step pattern).

### 5.3 Switch Procedure

```typescript
function switchToStiff(): void {
  isStiff = true;
  q = 1;
  stepsAtCurrentOrder = 0;
  prevDeltaValid = false;
  jacobianValid = false;
  acceptedBdfSteps = 0;
  // Reinitialize Nordsieck at order 1:
  nordsieckInit(y, dydt, h, nordsieck, dim);
  stepsSinceSwitch = 0;
}

function switchToAdams(): void {
  isStiff = false;
  q = 1;
  stepsAtCurrentOrder = 0;
  prevDeltaValid = false;
  consecutiveAdamsRejects = 0;
  nordsieckInit(y, dydt, h, nordsieck, dim);
  stepsSinceSwitch = 0;
}
```

### 5.4 Cooldown After Switch

For the first `MAX_SWITCH_COOLDOWN` (3) steps after a switch, clamp step growth to `MAX_SWITCH_GROW` (2.0×). This prevents oscillation between modes.

### 5.5 Testing Checkpoint

- Van der Pol oscillator (μ = 1000): confirm Adams → BDF switching during fast transient.
- `y' = -1000y + sin(t)`: confirm quick switch to BDF, then gradual return to Adams as the transient dies.

---

## Phase 6: Initial Step Size

Use `h = odes.arg.step` as initial step, matching the MRT convention. The adaptive loop will shrink/grow as needed.

Optional enhancement (Hairer-Wanner algorithm) for stiff problems where the user-provided step may be too large:

```typescript
// 1. d0 = max_i(|y0[i]| / yScale[i]), d1 = max_i(|f0[i]| / yScale[i])
// 2. h0 = 0.01 * d0 / d1 (or small default if d0 or d1 ≈ 0)
// 3. Trial Euler: y1 = y0 + h0 * f0, f1 = f(t0 + h0, y1)
// 4. d2 = max_i(|(f1[i] - f0[i])| / yScale[i]) / h0
// 5. h1 = (0.01 / max(d1, d2))^(1/2)
// 6. h = min(100*h0, h1, hDataframe, t1-t0)
```

Implement only if tests show the user-provided step causes issues.

---

## Phase 7: Main Integration Loop

### 7.1 Buffer Preallocation

Follow MRT pattern — all buffers declared before the main loop:

```typescript
// --- Common buffers (same as MRT) ---
const y = new Float64Array(odes.initial);
const yPrev = new Float64Array(odes.initial);
const dydt = new Float64Array(dim);
const yScale = new Float64Array(dim);
const yTemp = new Float64Array(dim);
const yErr = new Float64Array(dim);

// --- Nordsieck buffers ---
const maxRows = ADAMS_MAX_ORDER + 1; // 13
const nordsieck = new Float64Array(maxRows * dim);
const nordsieckPred = new Float64Array(maxRows * dim);
const delta = new Float64Array(dim);
const deltaPrev = new Float64Array(dim);

// --- BDF / Newton buffers ---
const Ident = new Float64Array(dimSquared);
for (let i = 0; i < dim; ++i) Ident[i + i * dim] = 1.0;
const W = new Float64Array(dimSquared);
const Jac = new Float64Array(dimSquared);
const L = new Float64Array(dimSquared);
const U = new Float64Array(dimSquared);
const luBuf = new Float64Array(dim);
const f0Buf = new Float64Array(dim);
const f1Buf = new Float64Array(dim);
const G = new Float64Array(dim);
const yNewton = new Float64Array(dim);
const fNewton = new Float64Array(dim);
const toUseLU = dim > 2;

// --- RKF45 buffers (for trial non-stiffness test) ---
const k1 = new Float64Array(dim);
const k2 = new Float64Array(dim);
const k3 = new Float64Array(dim);
const k4 = new Float64Array(dim);
const k5 = new Float64Array(dim);
const k6 = new Float64Array(dim);
```

### 7.2 Main Loop Pseudocode

```
// Initialize
f(t, y, dydt);
nordsieckInit(y, dydt, h, nordsieck, dim);
tArr[0] = t0; for (i) yArrs[i][0] = y[i];

while (flag) {
  f(t, y, dydt);
  if (callback) callback.onIterationStart();
  for (i) yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;
  if (t + h > t1) { h = t1 - t; flag = false; }

  while (true) {  // adaptive step
    // PREDICT
    nordsieckPredict(nordsieck, nordsieckPred, q, dim);

    if (!isStiff) {
      // ADAMS PECE
      // ... evaluate, correct, estimate error ...
    } else {
      // BDF NEWTON
      // ... jacobian check, W = I - l[0]*h*J, LU, iterate ...
    }

    if (errmax > 1) {
      // REJECT — shrink h, maybe switch to stiff
      continue;
    }

    // ACCEPT
    nordsieckCorrect(nordsieckPred, nordsieck, l, delta, q, dim);
    for (j) y[j] = nordsieck[j];
    t += h;

    // ORDER SELECTION → qNew, etaBest
    // RESCALE nordsieck if h changes
    // STIFFNESS CHECK (periodic)
    // POST-ACCEPT GUARDS
    break;
  }

  // INTERPOLATION to output grid (identical to MRT)
  while (timeDataframe < t) { ... }

  h = hNext; tPrev = t;
  for (i) yPrev[i] = y[i];
}

// FINAL POINT + CALLBACK + RETURN (identical to MRT)
```

### 7.3 Testing Checkpoint

- `y' = -y`, y(0)=1 on [0, 10]: error within tolerance.
- Robertson problem on [0, 1e7]: method switching, reasonable step count.
- Lotka-Volterra (non-stiff oscillatory): Adams stays active, high order.

---

## Phase 8: Edge Cases and Robustness

### 8.1 Step Size Floor

```typescript
tNew = t + h;
if (tNew == t)
  throw new Error(ERROR_MSG.LSODA_FAILS);
```

### 8.2 Newton Failure

- Halve `h`, invalidate Jacobian, drop order to 1 if order > 1, retry.

### 8.3 Repeated Failures

- 3 consecutive rejects → force Jacobian recomputation.
- 7 consecutive rejects → drop order to 1.
- 10 consecutive rejects → `throw new Error(ERROR_MSG.LSODA_FAILS)`.

### 8.4 Step Size Cap

```typescript
const hMax = hDataframe * 10;
hNext = Math.min(hNext, hMax);
```

---

## Implementation Order

Execute sequentially. Each phase ends with a testing checkpoint — **do not proceed until tests pass.**

| Step | Phase | Deliverable |
|------|-------|-------------|
| 1 | Phase 1 | Constants and coefficient arrays at top of file |
| 2 | Phase 2 | Nordsieck helper functions |
| 3 | Phase 3 | Single step logic (Adams PECE + BDF Newton) |
| 4 | Phase 7 | Main loop with **fixed order 1 only** — test on `y' = -y` |
| 5 | Phase 4 | Order selection — verify order ramps up |
| 6 | Phase 5 | Stiffness detection + switching — test on Van der Pol |
| 7 | Phase 6 | Initial step size refinement (if needed) |
| 8 | Phase 8 | Edge cases and robustness hardening |

---

## Key Numerical Pitfalls

- **Nordsieck rescaling exponents:** `z[i] *= eta^i` — get the exponent right for each row.
- **Pascal matrix direction:** use separate `zPred` buffer to avoid overwriting values still needed.
- **BDF order > 2 stability:** BDF3–5 are not A-stable. If eigenvalues are near the imaginary axis, limit BDF order.
- **Adams functional iteration divergence:** if the problem is stiff, functional iteration diverges. Count divergent iterations as a stiffness signal.
- **Jacobian perturbation:** the library's `jacobian()` uses a fixed `EPS = 1e-10`. Accept this for now.
- **`solve1d2d` for dim ≤ 2:** always use this fast path instead of LU, matching the library convention.
- **Zero allocations in the loop:** no `new Float64Array(...)` inside `while(flag)`.
