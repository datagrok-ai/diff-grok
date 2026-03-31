# IVP File Syntax Reference

This document describes the IVP file format used by Datagrok Diff Studio. It serves as the input specification for the IVP-to-LaTeX converter.

## File structure

An IVP file is plain text composed of **blocks**. Each block begins with a line starting with `#keyword`. Content lines within a block are indented with spaces.

## Comments

- `//` starts a single-line comment. Everything from `//` to the end of the line is ignored.
- `#comment:` block — multi-line free-form text (not parsed as math).

```
// This is a comment
dx/dt = x + y  // inline comment
```

**Important:** `//` is always a comment. IVP has no string literals. A single `/` is division.

## Block reference

### `#name`
Single-line model identifier.
```
#name: Robertson
```

### `#description`
Single-line model description.
```
#description: Robertson's chemical reaction model
```

### `#comment`
Multi-line free-form text.
```
#comment:
  Source: https://example.com
  This model describes...
```

### `#equations`
System of ODEs. Each equation has the form `d(var)/dt = expression` or `dvar/dt = expression`.

```
#equations:
  dx/dt = -k1 * x + k2 * y**2
  d(FFox)/dt = -E11 + E12
```

Formulas can span multiple lines — continuation lines are indented and do NOT start with `d.../dt` or contain a top-level `=`:
```
#equations:
  d(MEAthiol)/dt = 2 * (-E11 + E12 - E21 + E22 + E31 + E41 - E32 - E42 - E62 - ktox * E71 * E72)
                   - (MEAthiol + MA) * (Fin + Fper) / VL
```

### `#expressions`
Auxiliary formulas: `name = expression`.

```
#expressions:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2
```

Expressions may also contain:
- **Aliases:** `ceil = Math.ceil`
- **Arrow functions:** `func1 = (p, t) => expr`
- **Ternary operators:** `control = (P1 < 3) ? f(t) : g(t)`

### `#argument`
Independent variable definition with optional stage name.

```
#argument: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10}
  finish = 10 {caption: Final; category: Time}
  step = 0.01 {caption: Step; category: Time}
```

With stage name:
```
#argument: t, 1-st stage
  t0 = 0.01
  t1 = 15
  h = 0.01
```

### `#inits`
Initial values for ODE variables.

```
#inits:
  x = 2 {units: mol/L; category: Initial values; min: 0; max: 5} [Initial x]
  y = 0 {units: mol/L; category: Initial values; min: -2; max: 2} [Initial y]
```

### `#parameters`
Model parameters (generate UI controls).

```
#parameters:
  k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}
  k2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}
```

### `#constants`
Fixed values (no UI controls).

```
#constants:
  C1 = 1
  C2 = 3
```

### `#output`
Specifies which variables appear in output, with optional captions.

```
#output:
  t {caption: Time, h}
  A1 {caption: Central}
```

### `#loop`
Cyclic process. `count` sets the number of cycles.

```
#loop:
  count = 10 {caption: count; category: Dosing; min: 1; max: 20}
  depot += dose
```

### `#update`
Multi-stage process. `duration` sets stage length.

```
#update: 2-nd stage
  duration = overall - _t1
  S += 70
```

### `#tolerance`
Solver tolerance (single numeric value).

```
#tolerance: 1e-7
```

### `#meta.solver`
Solver configuration in JSON-like syntax.

```
#meta.solver: {method: 'mrt'; maxTimeMs: 50}
```

### `#meta.inputs`
Lookup table reference.

```
#meta.inputs: mode {caption: Process mode; choices: OpenFile("path.csv")}
```

## Annotation syntax

Values in `#inits`, `#parameters`, `#argument` can have annotations in `{...}` and tooltips in `[...]`:

```
name = value {key1: val1; key2: val2; ...} [Tooltip text]
```

Supported annotation keys:
- `units` — measurement units
- `caption` — display label
- `category` — UI grouping
- `min`, `max` — range for sliders
- `step` — slider step size

## Expression syntax

Expressions use JavaScript-like math syntax:

| Feature | Syntax | Example |
|---|---|---|
| Addition | `+` | `a + b` |
| Subtraction | `-` | `a - b` |
| Multiplication | `*` | `a * b` |
| Division | `/` | `a / b` |
| Exponentiation | `**` | `x**2`, `(a+b)**n` |
| Function call | `f(x)` | `sin(t)`, `pow(x, 2)` |
| Ternary | `c ? a : b` | `(x < 0) ? 0 : x` |
| Arrow function | `(args) => expr` | `(p, t) => p + t` |
| Comparison | `< > <= >= == !=` | `x >= 0` |
| Scientific notation | `1e4`, `1.2E-3` | `9.2E-2` |
| Grouping | `(...)` | `(a + b) * c` |

Available math functions: `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `pow`, `abs`, `Math.ceil`, `Math.floor`, `Math.ceil` (can be aliased).
