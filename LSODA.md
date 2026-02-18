# Plan: Implement LSODA Method

## Context

LSODA (Livermore Solver for Ordinary Differential equations with Automatic method switching) automatically detects stiffness and switches between Adams methods (non-stiff) and BDF methods (stiff). This gives users a single solver that handles both stiff and non-stiff problems without requiring them to know the problem's character in advance.

The diff-grok library already has all the building blocks: Adams predictor-corrector methods (AB4/AB5), implicit solver infrastructure (Jacobian, LU decomposition), and adaptive step control. LSODA unifies these under automatic stiffness detection.

**References:**
- Hindmarsh, A.C. and Petzold, L.R., "LSODA, Ordinary Differential Equation Solver for Stiff or Non-Stiff System", ODEPACK, 1983.
- Hairer, E. and Wanner, G., "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems", Springer, 1996.

## Design Decisions

- **Fixed-order** approach (Adams order 4 + BDF order 2) rather than full variable-order Nordsieck arrays — keeps implementation consistent with existing codebase patterns (~400 lines vs ~1500+)
- **Non-stiff mode**: Adams-Bashforth-Moulton order 4 (identical to existing AB4 method), bootstrapped with RKF45
- **Stiff mode**: BDF2 with Modified Newton iteration, using existing `jacobian()`, `luDecomp()`, `luSolve()`, `solve1d2d()`
- **Stiffness detection**: consecutive step rejection counter (Adams→BDF), periodic trial explicit step (BDF→Adams)
- **Classification**: hybrid method — added to both `methods` (all correctness tests) and `implicitMethods` (stiff performance tests)

---

## Step 1: Create `src/solver-tools/lsoda-method.ts`

### Constants

**Adams-Bashforth-Moulton order 4** (same as [ab4-method.ts](src/solver-tools/ab4-method.ts)):
```
AB: 55/24, -59/24, 37/24, -9/24
AM: 9/24, 19/24, -5/24, 1/24
MILNE = 19/270
ADAMS_HIST_LEN = 4
```

**BDF2**:
```
y_{n+1} = (4/3)*y_n - (1/3)*y_{n-1} + (2/3)*h*f(t_{n+1}, y_{n+1})
BDF2_A0 = 4/3, BDF2_A1 = -1/3, BDF2_GAMMA = 2/3
```

**BDF2 adaptive step exponents** (order 2):
```
BDF_PSHRNK = -0.5       (= -1/p)
BDF_PSGROW = -1/3        (= -1/(p+1))
```

**Newton iteration**:
```
MAX_NEWTON_ITER = 4
NEWTON_TOL = 0.33         (convergence rate threshold)
```

**Stiffness detection**:
```
STIFF_DETECT_THRESHOLD = 5    (consecutive Adams rejects -> switch to BDF)
NONSTIFF_TEST_INTERVAL = 30   (accepted BDF steps -> trial Adams step)
```

**RKF45 Butcher tableau** (same as [ab4-method.ts](src/solver-tools/ab4-method.ts) lines 40-70).

### Buffer Allocation

Pre-allocate all `Float64Array` buffers (zero GC pressure):

| Category | Buffers |
|----------|---------|
| Common | `y`, `yPrev`, `dydt`, `yScale`, `yTemp`, `yErr` |
| Adams | `fBuf[0..3]` (f-history), `yCorr`, `fPred` |
| RKF45 bootstrap | `k1`..`k6` |
| BDF | `yBDF_prev` (y_{n-1}), `yNewton`, `fNewton`, `delta`, `G`, `yBDF1` |
| Linear algebra | `Ident`, `W`, `Jac`, `L`, `U`, `luBuf`, `f0Buf`, `f1Buf` |

### State Variables

```
isStiff: boolean              // current mode
adamsHistCount: number        // valid f-history points (0..4)
consecutiveAdamsRejects: number
acceptedBdfSteps: number
bdfHistCount: number          // 0=none, 1=BDF1 ready, 2=BDF2 ready
jacobianValid: boolean
jacobianAge: number           // steps since last J recomputation
lastAdamsH, lastBdfH: number  // for detecting step size changes
```

### Algorithm Structure

```
lsoda(odes, callback):
  Setup (identical to AB4 pattern)
  Initialize: f(t0, y0, fBuf[0]), adamsHistCount = 1

  Main loop (while flag):
    f(t, y, dydt)
    callback?.onIterationStart()
    yScale = |y| + h*|dydt| + TINY
    End-point check

    Adaptive step loop (while true):
      if NOT isStiff:
        Adams step (see Section below)
      else:
        BDF step (see Section below)

    Linear interpolation to dataframe grid (identical to all other methods)
    h = hNext, save yPrev

  Finalization (identical to all other methods)
```

### Adams Mode (Non-Stiff)

Follows [ab4-method.ts](src/solver-tools/ab4-method.ts) exactly, with added stiffness monitoring:

1. **History invalidation**: if `h !== lastAdamsH` and `adamsHistCount > 1`, reset to `adamsHistCount = 1`
2. **Bootstrap** (`adamsHistCount < 4`): RKF45 step, on accept shift fBuf and increment count
3. **PECE step** (`adamsHistCount >= 4`): AB4 predict -> evaluate -> AM3 correct -> Milne error
4. **On reject**: increment `consecutiveAdamsRejects`, if `>= 5` -> switch to BDF mode
5. **On accept**: reset `consecutiveAdamsRejects = 0`, compute hNext with `PSHRNK=-0.25`, `PSGROW=-0.2`, cap at `hMax`

### BDF Mode (Stiff)

Follows implicit method pattern from [mrt-method.ts](src/solver-tools/mrt-method.ts) for Jacobian/LU:

1. **Jacobian update**: recompute via `jacobian()` when `!jacobianValid || jacobianAge > 20 || h !== lastBdfH`
2. **Form W**: `W = I - gamma*h*J` where gamma = 2/3 (BDF2) or 1.0 (BDF1 fallback)
3. **LU decompose**: `luDecomp(W, L, U, dim)` for dim > 2, otherwise use `solve1d2d()`
4. **Newton predictor**: BDF2 extrapolation `y^(0) = 2*y_n - y_{n-1}`, or explicit Euler for BDF1
5. **Newton iteration** (up to 4 iterations):
   - Evaluate `f(t+h, y^(k), fNewton)`
   - Compute residual `G = -(y^(k) - BDF_rhs - gamma*h*fNewton)`
   - Solve `W * delta = G` via LU (or solve1d2d)
   - Update `y^(k+1) = y^(k) + delta`
   - Check convergence: `||delta/yScale|| < tolerance`
6. **Error estimation**: compare BDF2 vs BDF1 approximation:
   ```
   yBDF1[i] = y[i] + h * fNewton[i]
   yErr[i] = (yNewton[i] - yBDF1[i]) / 3.0
   ```
7. **Step control**: use `BDF_PSHRNK=-0.5`, `BDF_PSGROW=-1/3`; no hMax cap (stiff problems need large steps)
8. **On accept**: update BDF history (`yBDF_prev <- y`, `y <- yNewton`), increment `acceptedBdfSteps`
9. **On reject**: invalidate Jacobian, degrade `bdfHistCount` to 1 (use BDF1 next)
10. **Newton failure**: treat as step rejection, halve step, invalidate Jacobian

### Stiffness Detection

**Adams -> BDF** (non-stiff -> stiff):
- Track `consecutiveAdamsRejects` — incremented on every Adams step rejection
- When `>= STIFF_DETECT_THRESHOLD (5)`: call `switchToStiff()`
  - Set `isStiff = true`, `acceptedBdfSteps = 0`, `jacobianValid = false`
  - Initialize `yBDF_prev` from `yPrev` (already available from interpolation buffer)
  - Set `bdfHistCount = 2` (both y_n and y_{n-1} available)

**BDF -> Adams** (stiff -> non-stiff):
- After `NONSTIFF_TEST_INTERVAL (30)` accepted BDF steps, attempt a trial RKF45 step
- If trial error <= 1.0 (would be accepted): call `switchToAdams()`
  - Set `isStiff = false`, `consecutiveAdamsRejects = 0`
  - Reset `adamsHistCount = 1`, `f(t, y, fBuf[0])` — re-bootstrap needed
- Trial step does NOT advance the solution (diagnostic only)

---

## Step 2: Add Error Message

**File**: [src/solver-tools/solver-defs.ts](src/solver-tools/solver-defs.ts) (line ~161)

Add to `ERROR_MSG` enum:
```typescript
LSODA_FAILS = 'The LSODA method fails',
```

## Step 3: Export the Solver

**File**: [src/solver-tools/index.ts](src/solver-tools/index.ts) (after line 9)

```typescript
export {lsoda} from './lsoda-method';
```

**File**: [index.ts](index.ts) (line 1)

Add `lsoda` to the re-export list.

## Step 4: Add to Web Worker

**File**: [src/worker-tools/solving.ts](src/worker-tools/solving.ts)

- Add `lsoda` to the import (line 3)
- Add case `'lsoda': return lsoda;` to `getMethod()` switch

## Step 5: Add to Test Suite

**File**: [src/tests/test-defs.ts](src/tests/test-defs.ts)

- Import `lsoda` from `'../../index'`
- Add `['LSODA', lsoda]` to `methods` map (correctness tests on all 6 problems)
- Add `['LSODA', lsoda]` to `implicitMethods` map (performance tests on stiff benchmarks)

## Step 6: Update Documentation

**File**: [README.MD](README.MD) — Add LSODA to features and solving sections as automatic stiffness-detecting method.

**File**: [CLAUDE.md](CLAUDE.md) — Add `lsoda` to method lists (as a new **Automatic methods** category).

---

## Verification

```bash
# Clean compiled JS files
find src -name "*.js" -delete && find src -name "*.js.map" -delete && rm -f index.js index.js.map

# Build
npm run build

# Lint
npm run lint

# Run all tests
npm test
```

**Expected results:**
- **Correctness tests**: 6 problems x 9 methods = 54 tests passing (all errors < MAX_MAD = 0.1)
- **Performance tests**: 6 stiff benchmarks x 4 methods = 24 tests passing (all < 10s timeout)
- **Pipeline tests**: 3 tests passing (method-agnostic)
