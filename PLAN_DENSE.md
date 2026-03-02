# Dense Output for LSODA — Implementation Plan

## Goal

Enable retrieval of ODE solution values at arbitrary points **after** integration
completes, by saving Nordsieck array snapshots at every successful internal step
and using Taylor/Horner interpolation on them.

Two public use-cases:

1. **Array of points**: user provides a sorted array `tArray: Float64Array` — return
   solution at each point.
2. **Uniform grid**: user provides `tStart`, `tEnd`, `step` — return solution on
   the grid `tStart, tStart+step, tStart+2*step, …, tEnd`.

---

## 1. Data structures

### 1.1 `NordsieckSnapshot`

New interface in `src/common.ts`:

```ts
export interface NordsieckSnapshot {
  tn: number;           // right edge of the step
  h: number;            // step size used
  hu: number;           // last successful step size (for range check)
  nq: number;           // method order
  yh: Float64Array[];   // yh[0..nq], each of length neq (0-based copies)
}
```

All arrays are **0-based** (converted from internal 1-based at capture time) so
that downstream code does not need to know about the 1-based convention.

### 1.2 `DenseOutput`

New class in a new file `src/dense.ts` that owns a `NordsieckSnapshot[]` array
and exposes the two interpolation APIs.

---

## 2. Capturing snapshots

### 2.1 Low-level: modify `lsoda()` in `src/lsoda.ts`

Inside the main `while(true)` loop, right after `kflag === 0` (line 427), before
the itask-based return logic, capture a snapshot:

```ts
if (kflag === 0) {
  // --- dense output: capture snapshot ---
  if (ctx.snapshots) {
    ctx.snapshots.push(captureSnapshot(c, neq));
  }
  // ... existing code unchanged ...
}
```

Add a helper `captureSnapshot(c, neq)` that deep-copies the needed fields from
`LsodaCommon` into a `NordsieckSnapshot`.

### 2.2 Enable via `LsodaContext`

Add an optional field to `LsodaContext`:

```ts
snapshots: NordsieckSnapshot[] | null = null;
```

When `null` (default), no snapshots are collected — zero overhead for users who
don't need dense output.

### 2.3 High-level: `Lsoda` class

Add a method (or constructor option) to enable snapshot collection, and expose
a `getDenseOutput(): DenseOutput` method that wraps the collected snapshots.

---

## 3. Interpolation logic (`src/dense.ts`)

### 3.1 Core: `evaluateAt(t: number): number[]`

Given a single time `t`:

1. Binary-search `snapshots` by `tn` to find the snapshot whose interval
   `[tn - hu, tn]` contains `t`.
   Because snapshots are ordered by `tn`, use `upperBound(t)` on `tn` values.
2. Compute `s = (t - snap.tn) / snap.h`.
3. Horner evaluation (0-based version of `intdy`):
   ```
   dky[i] = yh[nq][i]
   for j = nq-1 down to 0:
       dky[i] = yh[j][i] + s * dky[i]
   ```
4. Return `dky` as `number[]`.

### 3.2 Case 1: `solveAtTimes(tArray: number[]): number[][]`

```ts
solveAtTimes(tArray: number[]): number[][]
```

- Input: sorted ascending array of query points.
- Validate: all points within `[snapshots[0].tn - snapshots[0].hu, snapshots[last].tn]`.
- Optimization: since both `tArray` and `snapshots` are sorted, walk through
  them with a single pointer (no binary search per point). For each query point,
  advance the snapshot pointer until the right interval is found.
- Return: array of solution vectors (0-based, `number[]`), one per query point.

### 3.3 Case 2: `solveOnGrid(tStart, tEnd, step): { t: number[], y: number[][] }`

```ts
solveOnGrid(tStart: number, tEnd: number, step: number): { t: number[], y: number[][] }
```

- Generate the uniform grid `tStart, tStart+step, …` (clamp last point to `tEnd`).
- Delegate to `solveAtTimes` internally.
- Return both the `t` array and the corresponding `y` matrix.

---

## 4. File changes summary

| File | Change |
|------|--------|
| `src/common.ts` | Add `NordsieckSnapshot` interface; add `snapshots` field to `LsodaContext` |
| `src/lsoda.ts` | Add `captureSnapshot()` helper; insert snapshot capture in main loop; add dense output methods to `Lsoda` class |
| `src/dense.ts` | **New file**. `DenseOutput` class with `evaluateAt`, `solveAtTimes`, `solveOnGrid` |
| `src/index.ts` | Re-export `DenseOutput`, `NordsieckSnapshot`, and new `Lsoda` methods |
| `test/test.ts` | Add tests for both use-cases |

---

## 5. Testing strategy

### 5.1 Consistency test

Run the Robertson problem. Collect dense output. Query at the same `tout` points
used in the existing test. Verify that dense-output values match the direct
`solver.solve()` values to within a tight tolerance (e.g., `1e-10` relative).

### 5.2 Monotonicity / ordering test

Query at 1000 uniform points. Verify that returned `t` values are strictly
increasing and all `y` values are finite.

### 5.3 Boundary test

Query at `t = tStart` (left boundary) and `t = tEnd` (right boundary) —
verify correct results and no out-of-range errors.

### 5.4 Grid test

Use `solveOnGrid(0, 4e10, 1e9)`, verify length = 41, first and last `t` values,
and solution quality at sampled points.

---

## 6. Implementation order

1. Add `NordsieckSnapshot` interface and `snapshots` field to `LsodaContext`.
2. Implement `captureSnapshot()` and insert it into the `lsoda()` main loop.
3. Implement `DenseOutput` class with `evaluateAt`.
4. Implement `solveAtTimes` (linear scan).
5. Implement `solveOnGrid` (delegates to `solveAtTimes`).
6. Wire up `Lsoda` class: add `enableDenseOutput()` / constructor option and
   `getDenseOutput()`.
7. Update exports in `src/index.ts`.
8. Write tests.
