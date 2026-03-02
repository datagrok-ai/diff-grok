# LSODA C-to-TypeScript Transformation Patterns

This document captures the systematic patterns used to transform the C library
`liblsoda` (28 source files) into the TypeScript module in `src/` (15 files).

---

## 1. File & Module Structure

### Pattern: Many C files → Fewer TS modules with ES exports

| C (liblsoda/src/)            | TS (src/)            | Notes                                      |
|------------------------------|----------------------|--------------------------------------------|
| `lsoda.h`, `lsoda.c`        | `lsoda.ts`           | Public API + main solver                   |
| `common.h`, `common.c`      | `common.ts`          | Constants, types, context class            |
| `blas.h`, `dgefa.c`, `dgesl.c`, `daxpy.c`, `ddot.c`, `dscal.c`, `idamax.c`, `fnorm.c`, `vmnorm.c` | `blas.ts` | All BLAS routines consolidated into one module |
| `stoda.c`                    | `stoda.ts`           | 1-to-1                                     |
| `correction.c`               | `correction.ts`      | 1-to-1                                     |
| `prja.c`                     | `prja.ts`            | 1-to-1                                     |
| `methodswitch.c`             | `methodswitch.ts`    | 1-to-1 (static data moved to `common.ts`)  |
| `cfode.c`, `cfode_static.c`  | `cfode.ts`          | Merged                                     |
| `intdy.c`                    | `intdy.ts`           | 1-to-1                                     |
| `corfailure.c`               | `corfailure.ts`      | 1-to-1                                     |
| `orderswitch.c`              | `orderswitch.ts`     | 1-to-1                                     |
| `scaleh.c`                   | `scaleh.ts`          | 1-to-1                                     |
| `solsy.c`                    | `solsy.ts`           | 1-to-1                                     |
| `strdup_printf.c`, `printcf.c` | *(removed)*        | Debug/print utilities not needed           |
| *(none)*                     | `dense.ts`           | New: dense output interpolator (TS-only)   |
| *(none)*                     | `index.ts`           | New: barrel re-exports for public API      |

**Rule:** Header files (`.h`) are eliminated — their declarations become TS
`export`/`import` statements. One header can become multiple TS imports.
Small related C files (all BLAS routines) are consolidated into a single TS module.

---

## 2. Data Structures

### Pattern: C structs → TS classes/interfaces

| C                            | TypeScript                     |
|------------------------------|--------------------------------|
| `struct lsoda_context_t`     | `class LsodaContext`           |
| `struct lsoda_opt_t`         | `interface LsodaOpt`           |
| `struct lsoda_common_t`      | `class LsodaCommon`            |

**Mutable bags of state** (like `lsoda_common_t`) become **classes** with
all fields initialized to zero/empty defaults. Immutable option records become
**interfaces**.

```c
// C
struct lsoda_common_t {
    double **yh, **wm, *ewt, *savf, *acor;
    int    *ipvt;
    void   *memory;
    double  h, hu, rc, tn;
    int     ialth, ipup, nslp;
    // ...
};
```

```typescript
// TypeScript
export class LsodaCommon {
  yh: Float64Array[] = [];
  wm: Float64Array[] = [];
  ewt: Float64Array = new Float64Array(0);
  savf: Float64Array = new Float64Array(0);
  acor: Float64Array = new Float64Array(0);
  ipvt: Int32Array = new Int32Array(0);
  h: number = 0;
  hu: number = 0;
  rc: number = 0;
  tn: number = 0;
  ialth: number = 0;
  ipup: number = 0;
  nslp: number = 0;
  // ...
}
```

**Rules:**
- Every field gets an explicit default (`= 0`, `= new Float64Array(0)`, etc.).
- The `void *memory` field is dropped entirely (no manual memory management).
- `char *error` becomes `string | null = null`.

---

## 3. Function Pointer Types

### Pattern: C typedef for function pointer → TS type alias

```c
// C
typedef int (*_lsoda_f)(double, double *, double *, void *);
```

```typescript
// TypeScript
export type LsodaFunc = (t: number, y: Float64Array, ydot: Float64Array, data: any) => number;
```

**Rule:** `void *` user data becomes `any`. Return type stays `number`.

---

## 4. The `_C(x)` Macro

### Pattern: Macro accessor → Local alias `const c = ctx.common!`

The C code uses `#define _C(x) (ctx->common->x)` pervasively to access fields
of the internal state through the context pointer. In TypeScript, a local
variable with a non-null assertion replaces this:

```c
// C
_C(nje)++;
_C(h) * _C(el)[1]
_C(wm)[i][j] = (_C(acor)[i] - _C(savf)[i]) * fac;
```

```typescript
// TypeScript
const c = ctx.common!;
c.nje++;
c.h * c.el[1]
c.wm[i][j] = (c.acor[i] - c.savf[i]) * fac;
```

**Rule:** Every function that accesses `ctx->common` in C starts with
`const c = ctx.common!;` in TS. All `_C(field)` → `c.field`.

---

## 5. Arrays & Memory

### 5a. Pattern: `double *` → `Float64Array`, `int *` → `Int32Array`

```c
double *ewt, *savf, *acor;
int    *ipvt;
```

```typescript
ewt: Float64Array = new Float64Array(0);
savf: Float64Array = new Float64Array(0);
acor: Float64Array = new Float64Array(0);
ipvt: Int32Array = new Int32Array(0);
```

### 5b. Pattern: `double **` (2D) → `Float64Array[]`

```c
double **yh;    // allocated as array of pointers to rows
double **wm;
```

```typescript
yh: Float64Array[] = [];   // array of Float64Array rows
wm: Float64Array[] = [];
```

### 5c. Pattern: Fixed-size C arrays → `Float64Array` of same size

```c
double el[14];
double elco[13][14], tesco[13][4];
```

```typescript
el: Float64Array = new Float64Array(14);
elco: Float64Array[] = [];   // 13 elements, each Float64Array(14)
tesco: Float64Array[] = [];  // 13 elements, each Float64Array(4)
```

### 5d. Pattern: Single `malloc` block → Individual typed-array allocations

The C code allocates one contiguous block and manually carves it into
sub-arrays via pointer arithmetic. TypeScript allocates each array separately:

```c
// C: single malloc, manual offset arithmetic
_C(memory) = malloc(offset);
_C(yh) = (double **)((char *)_C(memory) + yhoff);
_C(ewt) = (double *)((char *)_C(memory) + ewtoff);
for(i = 0; i <= lenyh; i++)
    _C(yh)[i] = (double *)((char *)_C(memory) + yh0off + i * (1+nyh) * sizeof(double));
```

```typescript
// TypeScript: individual allocations, GC handles deallocation
c.yh = new Array(lenyh + 1);
for (let i = 0; i <= lenyh; i++)
  c.yh[i] = new Float64Array(nyh + 1);
c.ewt = new Float64Array(nyh + 1);
c.savf = new Float64Array(nyh + 1);
c.acor = new Float64Array(nyh + 1);
c.ipvt = new Int32Array(nyh + 1);
```

**Rule:** Drop all `malloc`/`free`/`memset`/pointer arithmetic. Replace with
typed-array constructors. Deallocation is handled by GC (set to `null`).

---

## 6. Array Indexing: 1-Based Convention Preserved

### Pattern: Keep Fortran-style 1-based indexing internally

The C code inherits Fortran's 1-based indexing by decrementing pointers:
```c
y--;  // now y[1] refers to what was y[0]
const double *rtol = opt->rtol - 1;
```

The TypeScript code **preserves the 1-based convention** for all internal
arrays by allocating `size + 1` elements and leaving index 0 unused:

```typescript
this.internalY = new Float64Array(neq + 1); // 1-indexed
c.ewt = new Float64Array(nyh + 1);          // ewt[1..nyh] used
```

The user-facing `Lsoda` class converts between 0-based (user) and 1-based (internal):

```typescript
// 0-based user → 1-based internal
for (let i = 0; i < neq; i++)
  this.internalY[i + 1] = y[i];
// 1-based internal → 0-based user
for (let i = 0; i < neq; i++)
  out[i] = this.internalY[i + 1];
```

---

## 7. Pointer Arithmetic → `Float64Array.subarray()`

### Pattern: `ptr + offset` → `.subarray(offset)`

C passes shifted pointers to BLAS and the user function. TypeScript uses
`subarray()` which creates a view without copying:

```c
// C: pointer arithmetic
(*ctx->function)(_C(tn), y + 1, _C(savf) + 1, ctx->data);
dscal(n - k, t, a[k] + k, 1);
daxpy(n - k, t, a[k] + k, 1, a[i] + k, 1);
```

```typescript
// TypeScript: subarray views
ctx.func(c.tn, y.subarray(1), c.savf.subarray(1), ctx.data);
dscal(n - k, t, a[k].subarray(k), 1);
daxpy(n - k, t, a[k].subarray(k), 1, a[i].subarray(k), 1);
```

**Rule:** Every `array + offset` in C becomes `array.subarray(offset)` in TS.
For 2D arrays `a[k] + k` → `a[k].subarray(k)`.

---

## 8. Output Parameters (Pointers)

### Pattern: `double *out` / `int *out` → Interface object or return value

C uses pointer parameters for multiple outputs. TypeScript uses either:

**(a) An interface for mutable output shared between functions:**

```c
// C: pointer parameters modified in-place
int correction(..., double *del, double *delp, ..., int *m) {
    *m = 0;
    *del = 0.;
    // ...
    (*m)++;
}
```

```typescript
// TypeScript: mutable state object
export interface CorrectionState {
  del: number;
  delp: number;
  m: number;
}

export function correction(..., cs: CorrectionState, ...): number {
  cs.m = 0;
  cs.del = 0.;
  // ...
  cs.m++;
}
```

**(b) A return value (when only one output is needed):**

```c
// C: output via pointer
void dgefa(double **a, int n, int *ipvt, int *info) {
    *info = 0;
    // ...
}
```

```typescript
// TypeScript: return value
export function dgefa(a: Float64Array[], n: number, ipvt: Int32Array): number {
  let info = 0;
  // ...
  return info;
}
```

---

## 9. Preprocessor Macros

### 9a. Pattern: `#define CONST` → `export const`

```c
#define ETA 2.2204460492503131e-16
#define SQRTETA 1.4901161193847656e-08
#define CCMAX 0.3
```

```typescript
export const ETA = 2.2204460492503131e-16;
export const SQRTETA = 1.4901161193847656e-08;
export const CCMAX = 0.3;
```

### 9b. Pattern: `#define max/min` → `Math.max/Math.min`

```c
#define max(a, b) ((a) > (b) ? (a) : (b))
min(opt->mxordn, mord[1])
```

```typescript
Math.min(opt.mxordn, mord[1])
```

### 9c. Pattern: Multi-statement macros → Local helper functions

The C code uses macros for inline code blocks with `return` statements
(`hardfailure`, `softfailure`, `successreturn`, `intdyreturn`, `ewset`).
These become local functions inside the enclosing TS function:

```c
// C: macro that returns from the enclosing function
#define hardfailure(fmt,...) { \
    ERROR(fmt, ## __VA_ARGS__); \
    ctx->state = -3; \
    return ctx->state; \
}
// Used as:  hardfailure("[lsoda] error %s", msg);
```

```typescript
// TypeScript: local function with explicit return
function hardfailure(msg: string): { t: number; state: number } {
  ctx.error = msg;
  ctx.state = -3;
  return { t, state: ctx.state };
}
// Used as:  return hardfailure(`[lsoda] error ${msg}`);
```

**Key difference:** C macros expand inline and their `return` exits the
enclosing function automatically. TS local functions require the caller to
explicitly `return` the result.

### 9d. Pattern: `ewset` macro → Named function

```c
#define ewset(ycur) { \
    for (i = 1; i <= neq; i++) \
        _C(ewt)[i] = rtol[i] * fabs((ycur)[i]) + atol[i]; \
    for (i = 1; i <= neq; i++) \
        _C(ewt)[i] = 1. / _C(ewt)[i]; \
}
```

```typescript
function ewset(ctx: LsodaContext, ycur: Float64Array): void {
  const c = ctx.common!;
  const neq = ctx.neq;
  const rtol = ctx.opt!.rtol;
  const atol = ctx.opt!.atol;
  for (let i = 1; i <= neq; i++)
    c.ewt[i] = rtol[i] * Math.abs(ycur[i]) + atol[i];
  for (let i = 1; i <= neq; i++)
    c.ewt[i] = 1. / c.ewt[i];
}
```

---

## 10. Standard Library Mapping

| C                           | TypeScript                    |
|-----------------------------|-------------------------------|
| `fabs(x)`                  | `Math.abs(x)`                 |
| `fmax(a, b)`               | `Math.max(a, b)`              |
| `fmin(a, b)`               | `Math.min(a, b)`              |
| `sqrt(x)`                  | `Math.sqrt(x)`                |
| `pow(x, y)`                | `Math.pow(x, y)`              |
| `fprintf(stderr, ...)`     | `console.error(...)`          |
| `printf(...)` / `ERROR(…)` | `ctx.error = \`template\``    |
| `malloc(n)`                | `new Float64Array(n)` / `new Array(n)` |
| `calloc(1, sizeof(T))`     | `new LsodaCommon()`          |
| `free(ptr)`                | `ctx.common = null` (GC)      |
| `memset(..., 0, ...)`      | `.fill(0)` on each typed array |

---

## 11. Variable Declarations

### Pattern: Hoist → Block-scoped `let`/`const`

C declares all variables at the top of a function. TypeScript declares them
at point of first use with appropriate scope:

```c
// C: all at top
int i, ier, j;
double fac, hl0, r, r0, yj;
const int neq = ctx->neq;
```

```typescript
// TypeScript: at point of use
const neq = ctx.neq;
const hl0 = c.h * c.el[1];
let fac = vmnorm(neq, c.savf, c.ewt);
let r0 = 1000. * Math.abs(c.h) * ETA * neq * fac;
```

**Rules:**
- Loop variables: `for (let i = 1; ...)` (not hoisted)
- Immutable intermediates: `const`
- Mutated values: `let`
- C `int` → TS `number` (no separate integer type)

---

## 12. `static` Data

### Pattern: File-scoped `static` arrays → Module-level `export const`

```c
// C: methodswitch.c
static double cm1[13] = {0x0p+0, 0x1p+1, ...};
static double cm2[13] = {0x0p+0, 0x1p+1, ...};
```

```typescript
// TypeScript: common.ts (centralized)
export const cm1: number[] = [0, 2.0, 5.999999999999998, ...];
export const cm2: number[] = [0, 2.0, 1.5, ...];
```

**Rules:**
- Hex float literals (`0x1p+1`) are converted to decimal (`2.0`).
- Static data scattered across C files is centralized in `common.ts`.
- TS uses `number[]` (not `Float64Array`) for compile-time constant arrays.

---

## 13. Error Handling

### Pattern: `fprintf(stderr, ...)` + state code → `ctx.error` string + return

```c
// C: variadic format strings
fprintf(stderr, "[lsoda] mxordn = %d is less than 0\n", opt->mxordn);
return 0;
```

```typescript
// TypeScript: template literals
ctx.error = `[lsoda] mxordn = ${opt.mxordn} is less than 0`;
return false;
```

**Rules:**
- Fatal errors: set `ctx.error` + set `ctx.state = -3` + return.
- Soft errors: set `ctx.error` + restore `y` from `yh[1]` + set state code + return.
- Format specifiers (`%d`, `%g`, `%s`) → template literal interpolation.
- `strdup_printf` helper eliminated entirely.

---

## 14. Control Flow

### Pattern: `while(1)` → `while(true)`

```c
while (1) { ... }
```
```typescript
while (true) { ... }
```

### Pattern: C-style cast → Implicit (unnecessary in TS)

```c
double conit = 0.5 / (double)(_C(nq) + 2);
```
```typescript
const conit = 0.5 / (c.nq + 2);
```

### Pattern: C truthiness `if (!ierpj)` → Same in TS

Integer truthiness tests (`0` = falsy) work identically in both languages, so
these are preserved as-is:

```c
if (!ierpj) { return corfailure(ctx, told); }
```
```typescript
if (!ierpj) { return corfailure(ctx, told); }
```

---

## 15. `#include` Guards → ES Module Imports

```c
// C
#ifndef _LSODA_H_
#define _LSODA_H_
#include "lsoda.h"
#include "blas.h"
#include "common.h"
```

```typescript
// TypeScript
import { LsodaContext, ETA, SQRTETA } from './common';
import { vmnorm, fnorm, dgefa } from './blas';
```

**Rule:** Named imports replace `#include`. Only the symbols actually used are
imported. No include guards needed (ES modules are singletons).

---

## 16. Public API Layering

### Pattern: Flat C API → Two-tier TS API (low-level + high-level class)

The C library exposes a flat procedural API:
```c
struct lsoda_context_t *ctx = lsoda_create_ctx();
lsoda_prepare(ctx, opt);
lsoda(ctx, y, &t, tout);
lsoda_free(ctx);
```

The TypeScript port preserves this as the **low-level API** and adds a
**high-level `Lsoda` class** on top:

```typescript
// Low-level (mirrors C)
export function lsodaPrepare(ctx, opt): boolean
export function lsoda(ctx, y, tIn, tout): { t, state }
export function lsodaReset(ctx): void
export function lsodaFree(ctx): void

// High-level (idiomatic TS)
export class Lsoda {
  constructor(f: OdeFunction, neq: number, opt?, data?)
  solve(y: ArrayLike<number>, t: number, tout: number): LsodaSolveResult
  getDenseOutput(): DenseOutput
  reset(): void
}
```

The high-level class handles:
- 0-based ↔ 1-based array conversion
- Option defaults and validation
- Dense output toggle (`{ dense: true }`)

---

## 17. New Capabilities (TS-Only)

### Dense Output (`dense.ts`)

The C library has no equivalent. The TS port adds a `NordsieckSnapshot`
mechanism: on each successful step, a snapshot of the Nordsieck array is
captured. The `DenseOutput` class then interpolates the solution at arbitrary
points using Horner's method.

This is an architectural addition enabled by the TS class system; it would
require significant refactoring to add to the C version.

### Barrel Exports (`index.ts`)

A standard TS pattern: re-exports all public symbols from a single entry point.

---

## Summary Checklist

When transforming a C function to TypeScript, apply these steps in order:

1. **File:** Map `.c` → `.ts` (consolidate small related files).
2. **Headers:** Replace `#include` with `import { ... } from '...'`.
3. **Signature:** `int func(ctx, ...)` → `export function func(ctx: LsodaContext, ...): number`.
4. **`_C(x)`:** Add `const c = ctx.common!;` then replace all `_C(x)` → `c.x`.
5. **Types:** `double` → `number`, `int` → `number`, `double *` → `Float64Array`.
6. **Vars:** Move declarations to point of use; use `const`/`let`.
7. **Pointer math:** `ptr + n` → `ptr.subarray(n)`.
8. **Output params:** `int *info` → return value; multiple outputs → interface object.
9. **Math:** `fabs` → `Math.abs`, `fmax` → `Math.max`, `pow` → `Math.pow`.
10. **Macros:** Constants → `export const`; code macros → local functions.
11. **Errors:** `fprintf(stderr, fmt, ...)` → `ctx.error = \`...\``; `console.error(...)`.
12. **Memory:** Drop `malloc`/`free`; use typed-array constructors; GC handles cleanup.
13. **Indexing:** Preserve 1-based internal convention; add 0↔1 conversion at API boundary.
