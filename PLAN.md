# Plan Template: Adding New ODE Solver Methods

This document provides a step-by-step guide for implementing new numerical solver methods in the diff-grok library.

## Overview

The diff-grok library currently supports:
- **Implicit methods** (for stiff ODEs): `mrt`, `ros3prw`, `ros34prw`
- **Explicit methods** (for non-stiff ODEs): `rk4`

All solvers follow a common architecture with adaptive step size control.

## Implementation Checklist

### 1. Create the Solver Implementation File

**Location**: `src/solver-tools/<method-name>-method.ts`

**Structure**:
```typescript
/* Method description and references */

import {ODEs, max, abs, SAFETY, PSHRNK, PSGROW, REDUCE_COEF, GROW_COEF,
  ERR_CONTR, TINY, EPS, tDerivative, jacobian, ERROR_MSG} from './solver-defs';
import {Callback} from './callbacks/callback-base';
import {luDecomp, luSolve, solve1d2d} from './lin-alg-tools'; // for implicit methods

// Method-specific constants (Butcher tableau coefficients, gamma values, etc.)
const COEF_1 = ...;
const COEF_2 = ...;

/** Solve initial value problem using the <Method Name> method
 * @param odes initial value problem for ordinary differential equations
 * @param callback computations control callback
 * @returns solution of the problem
*/
export function <methodName>(odes: ODEs, callback?: Callback): Float64Array[] {
  // 1. Extract problem parameters
  const f = odes.func;
  const t0 = odes.arg.start;
  const t1 = odes.arg.finish;
  let h = odes.arg.step;
  const hDataframe = h;
  const tolerance = odes.tolerance;

  // For explicit methods: add maximum step size constraint
  const hMax = hDataframe * 10;  // prevents poor interpolation

  // 2. Pre-allocate result arrays
  const rowCount = Math.trunc((t1 - t0) / h) + 1;
  const dim = odes.initial.length;
  const tArr = new Float64Array(rowCount);
  const yArrs = Array<Float64Array>(dim);
  for (let i = 0; i < dim; ++i)
    yArrs[i] = new Float64Array(rowCount);

  // 3. Initialize working buffers
  const y = new Float64Array(odes.initial);
  const yPrev = new Float64Array(odes.initial);
  const dydt = new Float64Array(dim);
  const yScale = new Float64Array(dim);
  const yTemp = new Float64Array(dim);
  const yErr = new Float64Array(dim);

  // Method-specific buffers (k-vectors, matrices, etc.)
  // For implicit: W matrix, L/U decomposition, Jacobian
  // For explicit: k1, k2, k3, ..., kN stage vectors

  // 4. Set initial values
  tArr[0] = t0;
  for (let i = 0; i < dim; ++i)
    yArrs[i][0] = y[i];

  // 5. Main adaptive stepping loop
  let timeDataframe = t0 + hDataframe;
  let t = t0;
  let tPrev = t0;
  let hNext = 0.0;
  let flag = true;
  let index = 1;

  while (flag) {
    // Compute derivative for scale vector
    f(t, y, dydt);

    // Callback hook
    if (callback)
      callback.onIterationStart();

    // Compute scale vector for error control
    for (let i = 0; i < dim; ++i)
      yScale[i] = abs(y[i]) + h * abs(dydt[i]) + TINY;

    // Check if we're approaching the endpoint
    if (t + h > t1) {
      h = t1 - t;
      flag = false;
    }

    // 6. Adaptive step inner loop
    while (true) {
      // === METHOD-SPECIFIC COMPUTATION ===
      //
      // For implicit methods (Rosenbrock-type):
      //   - Compute Jacobian: jacobian(t, y, f, EPS, f0Buf, f1Buf, W)
      //   - Compute time derivative: tDerivative(t, y, f, EPS, f0Buf, f1Buf, hdT)
      //   - Form W matrix: W = I - h*gamma*J
      //   - LU decomposition: luDecomp(W, L, U, dim)
      //   - Solve stages: luSolve(L, U, b, luBuf, k, dim)
      //
      // For explicit methods (Runge-Kutta):
      //   - Stage 1: f(t, y, k1)
      //   - Stage 2: f(t + c2*h, y + h*a21*k1, k2)
      //   - ... (N stages total)
      //   - Solution: yTemp = y + h*(b1*k1 + b2*k2 + ...)
      //   - Error: yErr = h*(e1*k1 + e2*k2 + ...) where ei = bi - bHat_i

      // === ERROR ESTIMATION ===
      let errmax = 0;
      for (let i = 0; i < dim; ++i)
        errmax = max(errmax, abs(yErr[i] / yScale[i]));
      errmax /= tolerance;

      // === STEP SIZE CONTROL ===
      if (errmax > 1) {
        // Step rejected: shrink step
        const hTemp = SAFETY * h * errmax**PSHRNK;
        h = max(hTemp, REDUCE_COEF * h);
        const tNew = t + h;
        if (tNew == t)
          throw new Error(ERROR_MSG.<METHOD>_FAILS);
      } else {
        // Step accepted: compute next step size
        if (errmax > ERR_CONTR)
          hNext = SAFETY * h * errmax**PSGROW;
        else
          hNext = GROW_COEF * h;

        // For explicit methods: cap maximum step
        if (hNext > hMax)
          hNext = hMax;

        t = t + h;
        for (let i = 0; i < dim; ++i)
          y[i] = yTemp[i];

        break;
      }
    } // adaptive step loop

    // 7. Linear interpolation to output grid
    while (timeDataframe < t) {
      const cLeft = (t - timeDataframe) / (t - tPrev);
      const cRight = 1.0 - cLeft;

      tArr[index] = timeDataframe;
      for (let j = 0; j < dim; ++j)
        yArrs[j][index] = cRight * y[j] + cLeft * yPrev[j];

      timeDataframe += hDataframe;
      ++index;
    }

    h = hNext;
    tPrev = t;
    for (let i = 0; i < dim; ++i)
      yPrev[i] = y[i];
  } // main loop

  // 8. Finalization
  if (callback)
    callback.onComputationsCompleted();

  tArr[rowCount - 1] = t1;
  for (let i = 0; i < dim; ++i)
    yArrs[i][rowCount - 1] = y[i];

  // 9. Return solution
  const solution = Array<Float64Array>(dim + 1);
  solution[0] = tArr;
  for (let i = 0; i < dim; ++i)
    solution[i + 1] = yArrs[i];

  return solution;
}
```

### 2. Add Error Message Enum

**File**: `src/solver-tools/solver-defs.ts`

Add to the `ERROR_MSG` enum:
```typescript
export enum ERROR_MSG {
  // ... existing entries
  <METHOD>_FAILS = 'The <Method Name> method fails',
};
```

### 3. Export the Solver

**File**: `src/solver-tools/index.ts`

Add export:
```typescript
export {<methodName>} from './<method-name>-method';
```

**File**: `index.ts` (root)

Add to the solver-tools export line:
```typescript
export {Func, ODEs, mrt, ros3prw, ros34prw, rk4, <methodName>, ...} from './src/solver-tools';
```

### 4. Add to Test Suite

**File**: `src/tests/test-defs.ts`

Import the method:
```typescript
import {mrt, ros3prw, ros34prw, rk4, <methodName>} from '../../index';
```

Add to appropriate methods map:

**For methods suitable for all problems** (correctness tests only):
```typescript
export const methods = new Map([
  // ... existing entries
  ['<MethodName>', <methodName>],
]);
```

**For implicit methods** (performance tests on stiff benchmarks):
```typescript
export const implicitMethods = new Map([
  // ... existing entries
  ['<MethodName>', <methodName>],
]);
```

### 5. Build, Lint, and Test

```bash
# Clean any compiled JS files (they interfere with ts-jest)
find src -name "*.js" -delete
find src -name "*.js.map" -delete
rm -f index.js index.js.map

# Build
npm run build

# Lint
npm run lint

# Test
npm test
```

**Expected test results**:
- **Correctness tests**: 6 problems × N methods (all methods in `methods` map)
- **Performance tests**: 6 stiff benchmarks × N methods (only `implicitMethods`)
- **Pipeline tests**: 3 tests (method-agnostic)

## Method Classification Guidelines

### Implicit Methods
- Designed for **stiff ODEs**
- Require Jacobian computation and linear system solves
- Examples: Rosenbrock-Wanner methods (MRT, ROS3PRw, ROS34PRw)
- **Include in**: `methods` + `implicitMethods` maps
- **Test with**: All correctness tests + all performance tests

### Explicit Methods
- Designed for **non-stiff ODEs**
- Only require function evaluations (no Jacobian)
- Examples: Runge-Kutta methods (RK4)
- **Include in**: `methods` map only
- **Test with**: All correctness tests only
- **Important**: Add `hMax` constraint to prevent poor interpolation

## Performance Considerations

1. **Pre-allocate all buffers** (Float64Array) to minimize GC pressure
2. **Use `solve1d2d()` for small systems** (dim ≤ 2) to avoid LU overhead
3. **Finite difference step**: Use `EPS = 1.0e-10` for Jacobian/derivative
4. **Adaptive step constants**: Reuse existing `SAFETY`, `PSHRNK`, `PSGROW`, etc.
5. **Linear interpolation**: Required because adaptive steps don't align with output grid

## Code Style

- **Indentation**: 2 spaces
- **Line length**: 120 characters max
- **Curly braces**: `"multi-or-nest"` (single-statement bodies don't need braces)
- **Comments**: Explain algorithm steps, reference equations/papers
- **Naming**: Use descriptive variable names (e.g., `yScale`, `errmax`, `hNext`)

## References

When implementing a method from literature:
1. Cite the paper/book in the file header
2. Document the Butcher tableau or method coefficients
3. Note any deviations from the reference (e.g., modified error estimation)

## Example: RK4 Implementation

See [src/solver-tools/rk4-method.ts](src/solver-tools/rk4-method.ts) for a complete example of an explicit Runge-Kutta method implementation.
