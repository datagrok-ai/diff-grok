# LSODA ODE Solver — TypeScript Implementation Plan

## Goal

Implement the LSODA solver from scratch in TypeScript with zero runtime dependencies (no external math libraries). The solver automatically switches between Adams-Bashforth-Moulton methods (orders 1–12) for non-stiff regions and BDF methods (orders 1–5) for stiff regions.

## Target API

```typescript
interface ODEResult {
  t: number[];          // time points
  y: number[][];        // solution at each time point, y[i] = state at t[i]
  nSteps: number;
  nFEval: number;
  nJEval: number;
  success: boolean;
  message: string;
}

interface ODEOptions {
  rtol?: number;            // relative tolerance, default 1e-6
  atol?: number | number[]; // absolute tolerance, default 1e-9
  maxSteps?: number;        // default 5000
  h0?: number;              // initial step size (auto if omitted)
  jac?: (t: number, y: number[]) => number[][]; // optional analytic Jacobian
  tEval?: number[];         // optional output time points
  maxOrder?: { adams: number; bdf: number }; // cap orders (defaults: 12, 5)
}

function lsoda(
  f: (t: number, y: number[]) => number[],
  tSpan: [number, number],
  y0: number[],
  options?: ODEOptions
): ODEResult;
```

---

## Phase 1: Linear Algebra Utilities

> No external dependencies. Implement the minimal set of dense linear algebra routines needed by the solver.

### 1.1 File: `src/linalg.ts`

| Function | Purpose |
|---|---|
| `luDecompose(A: number[][]): { L, U, P }` | PA = LU with partial pivoting |
| `luSolve(L, U, P, b: number[]): number[]` | Solve Ax = b given LU factors |
| `matVec(A: number[][], x: number[]): number[]` | Matrix-vector product |
| `matScale(A: number[][], s: number): number[][]` | Scalar × matrix |
| `matAdd(A, B): number[][]` | A + B |
| `eye(n: number): number[][]` | Identity matrix |
| `vecAxpy(a, x, y): number[]` | a·x + y (BLAS-style) |
| `vecNorm(x: number[]): number` | Euclidean norm |
| `vecScale(a, x): number[]` | a·x |

### 1.2 Testing Checkpoint
- Unit-test LU decomposition against known 3×3 and 5×5 systems.
- Verify `luSolve` round-trips correctly.
- Edge cases: singular matrix detection, 1×1 system.

---

## Phase 2: Core Data Structures

### 2.1 File: `src/types.ts`

```typescript
enum MethodType { ADAMS, BDF }

interface SolverState {
  t: number;
  h: number;
  q: number;               // current order
  method: MethodType;
  nordsieck: number[][];    // (q+1) × n Nordsieck array
  nSteps: number;
  nFEval: number;
  nJEval: number;
  nLU: number;
  // Jacobian cache
  jac: number[][] | null;
  jacTime: number;          // t at which Jacobian was last computed
  jacStale: boolean;
  // LU cache for iteration matrix M = I - h·β·J
  iterLU: { L: number[][]; U: number[][]; P: number[] } | null;
  iterH: number;            // h used when LU was computed
  iterGamma: number;        // h·β used when LU was computed
  // Error history for order selection
  errorHistory: number[];   // recent WRMS error norms
  // Stiffness detection
  stiffDetectCount: number;
  nSinceMethodSwitch: number;
}
```

### 2.2 File: `src/nordsieck.ts`

Nordsieck array helpers:

| Function | Purpose |
|---|---|
| `nordsieckInit(y0, f0, h, n)` | Build initial order-1 Nordsieck array `[y0, h·f0]` |
| `nordsieckPredict(z, q)` | Multiply by Pascal matrix: `z_pred = P·z` |
| `nordsieckCorrect(zPred, l, delta)` | Apply correction: `z[i] += l[i] * delta` |
| `nordsieckRescale(z, eta, q)` | Rescale after step size change: `z[i] *= eta^i` |
| `nordsieckChangeOrder(z, qOld, qNew)` | Adjust array dimensions for order change |

Pascal matrix multiplication should be done in-place without forming the full matrix:
```
z_pred[i] = Σ_{j=i}^{q} C(j,i) · z[j]
```
where `C(j,i)` is the binomial coefficient.

---

## Phase 3: Method Coefficients

### 3.1 File: `src/coefficients.ts`

Precompute all coefficients at module load time.

**Adams corrector (Moulton) coefficients and error constants:**
- For each order q = 1..12, store the correction vector `l[0..q]` in Nordsieck form.
- Store error constant `C_adams(q)`.
- Reference: Hairer & Wanner, or the original LSODE source (opkda1.f).

**BDF coefficients and error constants:**
- For each order q = 1..5, store:
  - `alpha[0..q]` — the BDF formula coefficients
  - `beta_s` — the coefficient on `h·f(t_n, y_n)`
  - The correction vector `l[0..q]` in Nordsieck form
  - Error constant `C_bdf(q)`
- Known BDF values for validation:

| Order | β_s | Error const |
|---|---|---|
| 1 | 1 | 1/2 |
| 2 | 2/3 | 2/9 |
| 3 | 6/11 | 3/22 |
| 4 | 12/25 | 12/125 |
| 5 | 60/137 | 10/137 |

**Provide a function:**
```typescript
function getCoefficients(method: MethodType, q: number): {
  l: number[];         // correction vector
  errCoeff: number;    // error constant C(q)
  gamma: number;       // β_s (the implicit coefficient h multiplier)
}
```

### 3.2 Testing Checkpoint
- Verify BDF coefficients reproduce known formulas for orders 1–5.
- Check that Adams correction vectors match reference values.

---

## Phase 4: Single Step — Predict-Correct Cycle

### 4.1 File: `src/step.ts`

This is the core computational kernel.

#### 4.1.1 `predict(state: SolverState): number[][]`
- Apply Pascal matrix to Nordsieck array → `zPred`.
- Return predicted Nordsieck array. `zPred[0]` is the predicted `y_n`.

#### 4.1.2 `correct(state, zPred, f, jac, options): StepResult`

```typescript
interface StepResult {
  success: boolean;
  yNew: number[];
  errNorm: number;          // WRMS norm of local error
  delta: number[];          // correction vector (y_corrected - y_predicted)
  converged: boolean;       // did Newton converge?
  nNewtonIter: number;
}
```

Newton iteration loop (max 4 iterations):

1. Compute residual: `g = y^(k) - zPred[0] - h·γ·f(t_new, y^(k))` where `γ` = β_s for BDF or the Adams corrector coefficient.
2. If Jacobian is stale or iteration matrix doesn't match current `h·γ`:
   - Recompute Jacobian (analytic or finite differences).
   - Form iteration matrix `M = I - h·γ·J`.
   - LU-decompose `M`. Cache it.
3. Solve `M · Δy = -g` using cached LU.
4. Update `y^(k+1) = y^(k) + Δy`.
5. Check convergence: `‖Δy‖_WRMS < 0.3`.
6. If converged: compute local error estimate `e = errCoeff · (y_corrected - zPred[0])`.

#### 4.1.3 WRMS Norm

```typescript
function wrmsNorm(v: number[], y: number[], atol: number[], rtol: number): number {
  let sum = 0;
  for (let i = 0; i < v.length; i++) {
    const wt = atol[i] + rtol * Math.abs(y[i]);
    sum += (v[i] / wt) ** 2;
  }
  return Math.sqrt(sum / v.length);
}
```

### 4.2 Testing Checkpoint
- Test a single step of BDF1 (implicit Euler) on `y' = -y`, verify against analytic solution.
- Test a single step of Adams order 1 on `y' = -y`.
- Confirm Newton convergence in ≤ 3 iterations for a smooth problem.

---

## Phase 5: Jacobian Computation

### 5.1 File: `src/jacobian.ts`

#### 5.1.1 Finite Difference Jacobian

```typescript
function finiteDiffJac(
  f: (t: number, y: number[]) => number[],
  t: number,
  y: number[],
  f0: number[],       // f(t, y) already evaluated
  atol: number[]
): number[][]
```

Column-by-column forward differences:
```
δ_j = sqrt(ε) · max(|y_j|, 1 / wt_j)
J[:, j] = (f(t, y + δ_j·e_j) - f0) / δ_j
```
where `ε = Number.EPSILON ≈ 2.2e-16`.

Cost: `n` extra function evaluations per Jacobian computation.

#### 5.1.2 Jacobian Refresh Strategy
- Recompute when:
  - Newton iteration fails to converge.
  - Step is rejected more than once in a row.
  - More than 20 steps since last Jacobian evaluation.
  - `h·γ` has changed by more than 30% since last LU factorization.
- Never recompute more than once per step attempt.

---

## Phase 6: Step Size and Order Selection

### 6.1 File: `src/adapt.ts`

#### 6.1.1 Step Size from Error Norm

For method of order `q` with error norm `errNorm`:

```typescript
function computeEta(errNorm: number, q: number, safety: number = 0.9): number {
  // eta = safety * (1 / errNorm)^(1/(q+1))
  return safety * Math.pow(1.0 / errNorm, 1.0 / (q + 1));
}
```

Clamp `eta` to `[0.2, 10.0]` for accepted steps, `[0.2, 0.9]` for rejected steps.

#### 6.1.2 Order Selection (after successful step)

Evaluate three candidates (if available):

| Candidate | Error estimate source | Formula |
|---|---|---|
| Order q−1 | From Nordsieck component `z[q]` | `eta_down = (1/errDown)^(1/q)` |
| Order q   | From step error (already computed) | `eta_same = (1/errNorm)^(1/(q+1))` |
| Order q+1 | From difference of consecutive corrections | `eta_up = (1/errUp)^(1/(q+2))` |

Rules:
- Order q+1 can only be considered if at least `q+1` steps have been taken at current order.
- Order q+1 cannot be considered immediately after a failed step.
- Choose the candidate with the largest `eta` (i.e., largest possible new step size).
- Apply safety factor to all candidates.

#### 6.1.3 Error Estimate for q−1

```typescript
// errDown ≈ ‖z[q]‖_WRMS / (C(q-1) error constant ratio)
```

This is cheap — just read an existing Nordsieck component.

#### 6.1.4 Error Estimate for q+1

```typescript
// errUp ≈ ‖currentCorrection - previousCorrection‖_WRMS scaled by C(q+1)
```

Requires storing the previous step's correction delta.

### 6.2 Testing Checkpoint
- On `y' = -y`, verify that the solver ramps up the order from 1 to ~4-5 over the first few steps.
- Verify step size grows smoothly on a non-stiff problem.

---

## Phase 7: Stiffness Detection and Method Switching

### 7.1 File: `src/stiffness.ts`

#### 7.1.1 Detecting Stiffness (Adams → BDF)

While running in Adams mode, monitor a stiffness indicator every ~10–15 steps:

- After each step, estimate the dominant eigenvalue of the Jacobian using the Newton convergence rate: `ρ ≈ ‖Δy^(k+1)‖ / ‖Δy^(k)‖`.
- If `h · ρ` exceeds the Adams stability boundary for the current order, increment a stiffness counter.
- Alternative / complementary: if the step size is being restricted far below what the error alone would allow, suspect stiffness.
- If the stiffness counter reaches a threshold (e.g., 5), switch to BDF.

#### 7.1.2 Detecting Non-Stiffness (BDF → Adams)

While running in BDF mode, periodically (every ~20 steps) check:

- Estimate `h · ρ` as above.
- If `h · ρ` is well within the Adams stability region for a reasonable order, increment a non-stiffness counter.
- If counter reaches threshold, switch to Adams.

#### 7.1.3 Method Switch Procedure

1. Save current `t`, `y_n`, `h`.
2. Reset order to `q = 1`.
3. Reinitialize the Nordsieck array for the new method at order 1: `z = [y_n, h·f(t, y_n)]`.
4. Set the method flag.
5. Allow the order selection mechanism to ramp up naturally.
6. Reset stiffness/non-stiffness counters.

### 7.2 Testing Checkpoint
- Van der Pol oscillator with μ = 1000: confirm solver starts in Adams mode for the slow phase, switches to BDF during the sharp transition.
- Linear ODE `y' = -1000y`: confirm immediate switch to BDF.

---

## Phase 8: Initial Step Size Selection

### 8.1 File: `src/initstep.ts`

If `h0` is not provided by the user, use the Hairer-Wanner algorithm:

1. Compute `d0 = ‖y0‖_WRMS` and `d1 = ‖f(t0, y0)‖_WRMS`.
2. First guess: `h0 = 0.01 · d0 / d1` (or `1e-6` if `d0` or `d1` is tiny).
3. Take an explicit Euler step: `y1 = y0 + h0·f0`.
4. Compute `f1 = f(t0 + h0, y1)`.
5. `d2 = ‖f1 - f0‖_WRMS / h0`.
6. `h1 = (0.01 / max(d1, d2))^(1/2)` (for order-1 start).
7. `h0 = min(100·h0, h1)`.
8. Ensure `h0 ≤ t_end - t0`.

---

## Phase 9: Main Integration Loop

### 9.1 File: `src/lsoda.ts`

```
function lsoda(f, tSpan, y0, options): ODEResult
```

**Pseudocode:**

```
1.  Parse options, set defaults.
2.  Initialize: t = t0, y = y0, compute f0 = f(t0, y0).
3.  Select initial step size h (Phase 8).
4.  Build Nordsieck array at order 1 for Adams method.
5.  While t < t_end and nSteps < maxSteps:
      a. PREDICT: zPred = nordsieckPredict(z, q)
      b. CORRECT: result = correct(state, zPred, f, jac, options)
      c. If Newton did not converge:
           - Reduce h by factor 0.25
           - Mark Jacobian stale
           - Retry from (a)
      d. If errNorm > 1.0 (step rejected):
           - Compute eta = computeEta(errNorm, q) (clamped to [0.2, 0.9])
           - h *= eta
           - If order > 1 and repeated failures: reduce order by 1
           - Retry from (a)
      e. Step accepted:
           - Update Nordsieck array with correction
           - Record solution if t_new is past a requested output point (interpolate)
           - Run order selection → get (qNew, etaNew)
           - Run stiffness detection (every N steps)
           - If method switch triggered: reinitialize (Phase 7.1.3)
           - Else: apply new order and step size
           - t = t_new, nSteps++
6.  Return collected solution points.
```

### 9.2 Dense Output / Interpolation

For output at arbitrary `tEval` points, use Nordsieck polynomial interpolation:

```typescript
function interpolate(z: number[][], t_n: number, h: number, t_out: number, q: number): number[] {
  const s = (t_out - t_n) / h;  // s ∈ [0, 1] within current step
  // Evaluate Nordsieck polynomial: y(t_out) ≈ Σ_{i=0}^{q} z[i] · s^i / (normalizations cancel)
  // Use Horner's method for stability
}
```

### 9.3 Testing Checkpoint
- Solve `y' = -y`, y(0)=1 on [0, 10]. Compare with exp(-t) — error should be within tolerance.
- Solve the Robertson chemical kinetics problem (classic stiff test):
  ```
  y1' = -0.04·y1 + 1e4·y2·y3
  y2' =  0.04·y1 - 1e4·y2·y3 - 3e7·y2²
  y3' =  3e7·y2²
  ```
  Confirm solver switches to BDF and solution stays physical (y_i ≥ 0).

---

## Phase 10: Edge Cases and Robustness

### 10.1 Minimum Step Size
- `h_min = 10 · ε · |t|` where `ε = Number.EPSILON`.
- If solver hits `h_min`, return with `success: false` and a diagnostic message.

### 10.2 Max Steps Guard
- Return partial solution with warning if `maxSteps` exceeded.

### 10.3 Repeated Failures
- If 3 consecutive rejected steps: force Jacobian recomputation.
- If 7 consecutive rejected steps: drop order to 1.
- If 10 consecutive rejected steps: return failure.

### 10.4 Overflow / NaN Protection
- After each function evaluation, check for NaN/Infinity in `f` output.
- If detected, halve step size and retry. After 3 NaN failures, return with error.

---

## Phase 11: Testing and Validation

### 11.1 File: `tests/`

**Unit tests** (per phase):
| Test | Validates |
|---|---|
| `linalg.test.ts` | LU decomposition, solve, edge cases |
| `nordsieck.test.ts` | Predict, correct, rescale operations |
| `coefficients.test.ts` | Known BDF/Adams coefficient values |
| `step.test.ts` | Single step on y' = -y |
| `jacobian.test.ts` | Finite diff Jacobian vs analytic for known systems |

**Integration tests** (end-to-end):
| Test Problem | Purpose |
|---|---|
| `y' = -y` | Basic exponential decay, non-stiff |
| `y' = -1000y + 3000 - 2000·exp(-t)` | Mildly stiff, verify BDF kicks in |
| Van der Pol (μ = 1000) | Classic stiff, verify method switching |
| Robertson kinetics | Very stiff 3-component system |
| Lotka-Volterra | Non-stiff oscillatory system, verify Adams stays active |
| `y' = cos(t)` with tEval | Verify dense output interpolation accuracy |

**Regression tests:**
- Compare against `scipy.integrate.solve_ivp(method='LSODA')` on all above problems.
- Error at final time should be within 10× of requested tolerance.

### 11.2 Performance Benchmark
- Robertson problem on `[0, 1e11]`: should complete in < 1 second, < 2000 steps.
- Van der Pol μ=1000 on `[0, 3000]`: verify reasonable step count vs SciPy.

---

## File Structure

```
src/
├── linalg.ts          # Phase 1: LU, matrix/vector ops
├── types.ts           # Phase 2: interfaces, enums
├── nordsieck.ts       # Phase 2: Nordsieck array operations
├── coefficients.ts    # Phase 3: Adams + BDF coefficients
├── jacobian.ts        # Phase 5: Jacobian computation
├── step.ts            # Phase 4: predict-correct cycle
├── adapt.ts           # Phase 6: step size + order selection
├── stiffness.ts       # Phase 7: stiffness detection + switching
├── initstep.ts        # Phase 8: initial step size
├── lsoda.ts           # Phase 9: main integration loop + interpolation
└── index.ts           # Public API export
tests/
├── linalg.test.ts
├── nordsieck.test.ts
├── coefficients.test.ts
├── step.test.ts
├── jacobian.test.ts
├── integration.test.ts  # End-to-end tests
└── benchmark.ts
tsconfig.json
package.json
```

---

## Implementation Order

Execute phases sequentially. Each phase has a testing checkpoint — **do not proceed to the next phase until all tests pass.**

1. **Phase 1** → Linear algebra utilities + tests
2. **Phase 2** → Data structures + Nordsieck helpers + tests
3. **Phase 3** → Coefficients + validation against known values
4. **Phase 4** → Single step (predict-correct) + single-step tests
5. **Phase 5** → Jacobian + tests
6. **Phase 6** → Step size / order adaptation + tests
7. **Phase 8** → Initial step size
8. **Phase 9** → Main loop with Adams only (no switching yet) → test on non-stiff problems
9. **Phase 7** → Stiffness detection + method switching → test on stiff problems
10. **Phase 10** → Edge cases + robustness
11. **Phase 11** → Full test suite + benchmarks

---

## Key Numerical Pitfalls to Watch For

- **Catastrophic cancellation** in `y_corrected - y_predicted` when correction is small relative to y. Use compensated summation if needed.
- **Nordsieck rescaling** after step size change must use `eta^i` for the i-th component — get the exponents right.
- **BDF order ≥ 3** are not A-stable. If solving a problem with eigenvalues near the imaginary axis, the solver may oscillate between Adams and BDF. The stiffness detection thresholds need tuning for this.
- **Jacobian perturbation size** `δ_j` must scale with both `|y_j|` and the tolerance. Too small → roundoff noise dominates; too large → truncation error dominates.
- **Pascal matrix** must be applied from high index to low (or use a temporary array) to avoid overwriting values that are still needed.
