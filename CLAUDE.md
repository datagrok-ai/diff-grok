# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**diff-grok** is a zero-dependency TypeScript library for solving initial value problems (IVPs) for ordinary differential equations (ODEs), focused on stiff equations. It implements Rosenbrock-Wanner type numerical methods and includes a declarative model scripting system.

## Build & Development Commands

```bash
npm install          # Install dev dependencies
npm run build        # Compile TypeScript (tsc)
npm run lint         # ESLint check: eslint "./src/**/*.ts"
npm run lint-fix     # Auto-fix lint issues
npm test             # Run Jest test suite
```

Run a single test file:
```bash
npx jest src/tests/correctness.test.ts
```

## Architecture

### Module Structure

- **`index.ts`** — Single entry point; re-exports everything from the four modules below.

- **`src/solver-tools/`** — Core numerical solvers
  - **Implicit methods** (for stiff ODEs): `mrt` (Modified Rosenbrock Triple), `ros3prw`, `ros34prw`
    - Require Jacobian computation and linear system solves
    - Use LU decomposition for solving W*k = b at each stage
  - **Explicit methods** (for non-stiff ODEs): `rk4` (Runge-Kutta-Fehlberg 4(5)), `ab5` (Adams-Bashforth-Moulton 5), `ab4` (Adams-Bashforth-Moulton 4), `rkdp` (Dormand-Prince 5(4))
    - Only require function evaluations (no Jacobian)
    - Include `hMax` constraint to prevent poor interpolation
    - `ab5` and `ab4` are multistep predictor-corrector methods bootstrapped with RKF45
  - Key types: `ODEs` (problem definition), `Func` (RHS function signature), `SolverMethod`
  - `solver-defs.ts` — Type definitions, adaptive step constants, Jacobian/derivative computation
  - `lin-alg-tools.ts` — Linear algebra routines (LU decomposition, matrix ops)
  - `callbacks/` — Flow control during solving (iteration limits, time limits)

- **`src/scripting-tools/`** — Declarative model parsing and JS code generation
  - `scripting-tools.ts` — Main parser: converts model DSL text → `IVP` object → JavaScript code
  - Model DSL uses `#name`, `#equations`, `#argument`, `#inits`, `#parameters`, `#constants`, `#expressions` blocks
  - Entry functions: `getIVP()` (parse model string), `getJScode()` (generate JS from IVP)

- **`src/pipeline/`** — Multi-stage simulation orchestration
  - Three pipeline creators for different model types:
    - `BasicModelPipelineCreator` — Single-stage models
    - `CyclicModelPipelineCreator` — Loop/dosing models (PK-PD)
    - `UpdatesModelPipelineCreator` — Multi-stage models with state transitions
  - `applyPipeline()` runs the full computation pipeline

- **`src/worker-tools/`** — WebWorker integration for browser parallel computation
  - `getIvp2WebWorker()` serializes IVP for worker transfer; `solveIvp()` runs in worker

### Data Flow

1. Define problem as `ODEs` object (programmatic) or parse model string via `getIVP()`
2. Call solver method (`mrt`/`ros3prw`/`ros34prw`/`rk4`/`ab5`/`ab4`/`rkdp`) → returns `Float64Array[]` (argument values + solution columns)
3. For complex models: create pipeline via `getPipelineCreator()` → `applyPipeline()`

### Performance Patterns

- Uses `Float64Array` throughout for memory efficiency
- Pre-allocated buffers to minimize GC pressure
- RHS function signature `(t, y, output) => void` writes to pre-allocated output array

## Code Style

- **ESLint config**: Google style base, 2-space indent, 120 char max line length
- **Curly braces**: `"multi-or-nest"` — single-statement bodies don't need braces, multi-line do
- **TypeScript strict mode** enabled
- Commit messages: `feat:`, `fix:`, `docs:` prefixes (conventional style)

## Testing

Tests are in `src/tests/` using Jest with ts-jest:
- `correctness.test.ts` — Validates solver accuracy against reference solutions (threshold: MAX_MAD = 0.1)
  - Tests all methods (implicit + explicit) against 6 problems with known exact solutions
  - 3 non-stiff problems (1D, 2D, 3D) + 3 stiff problems (1D, 2D, 3D)
- `performance.test.ts` — Benchmarks solver speed on stiff problems (timeout: 10,000ms)
  - Tests only implicit methods (MRT, ROS3PRw, ROS34PRw) — explicit methods not suitable for stiff benchmarks
  - Test problems: Robertson, HIRES, VDPOL, OREGO, E5, Pollution
- `pipeline.test.ts` — Pipeline integration tests (3 model types)
- Method definitions in `test-defs.ts`: `methods` map (all 7 solvers), `implicitMethods` map (stiff-capable only)

## Key Types

```typescript
type Func = (t: number, y: Float64Array, output: Float64Array) => void;
type ODEs = { name, arg: {name, start, finish, step}, initial, func: Func, tolerance, solutionColNames };
type SolverMethod = (odes: ODEs, callback?: Callback) => Float64Array[];
```
