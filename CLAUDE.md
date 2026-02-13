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
  - Three methods: `mrt` (Modified Rosenbrock Triple), `ros3prw`, `ros34prw`
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
2. Call solver method (`mrt`/`ros3prw`/`ros34prw`) → returns `Float64Array[]` (argument values + solution columns)
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
- `performance.test.ts` — Benchmarks solver speed (timeout: 10,000ms)
- `pipeline.test.ts` — Pipeline integration tests
- Test problems defined in `test-defs.ts` (Robertson, HIRES, VDPOL, OREGO, E5, Pollution)

## Key Types

```typescript
type Func = (t: number, y: Float64Array, output: Float64Array) => void;
type ODEs = { name, arg: {name, start, finish, step}, initial, func: Func, tolerance, solutionColNames };
type SolverMethod = (odes: ODEs, callback?: Callback) => Float64Array[];
```
