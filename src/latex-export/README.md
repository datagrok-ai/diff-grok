# latex-export

Converts [Diff Studio](https://datagrok.ai/help/compute/diff-studio) IVP model files into **LaTeX** and **Markdown** with LaTeX math. Part of the [diff-grok](../../README.md) library.

## Quick start

```typescript
import {convertIvpToLatex} from 'diff-grok';

const ivp = `#name: Robertson
#equations:
  dA/dt = -0.04 * A + 1e4 * B * C
  dB/dt = 0.04 * A - 1e4 * B * C - 3e7 * B**2
  dC/dt = 3e7 * B**2

#inits:
  A = 1 {units: mol/L}
  B = 0 {units: mol/L}
  C = 0 {units: mol/L}

#argument: t
  start = 0
  finish = 40
  step = 0.01`;

// LaTeX output (default)
const latex = convertIvpToLatex(ivp);

// Markdown output
const markdown = convertIvpToLatex(ivp, {format: 'markdown'});
```

## API

### `convertIvpToLatex(ivpText, options?)`

Converts the full text of an `.ivp` file to LaTeX or Markdown.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `ivpText` | `string` | Raw text content of an `.ivp` file |
| `options` | `Partial<ConvertOptions>` | Conversion options (all optional) |

**Returns:** `string` — LaTeX or Markdown markup.

### `ConvertOptions`

| Option | Type | Default | Description |
|---|---|---|---|
| `format` | `'latex' \| 'markdown'` | `'latex'` | Output format |
| `includeMetadata` | `boolean` | `true` | Include name, description, and comment |
| `includeInits` | `boolean` | `true` | Include initial conditions table |
| `includeParameters` | `boolean` | `true` | Include parameters table |
| `includeConstants` | `boolean` | `true` | Include constants table |
| `useCdot` | `boolean` | `true` | Use `\cdot` for multiplication; if `false`, use juxtaposition |
| `compact` | `boolean` | `false` | Compact mode: no section headings, inline initial conditions |

### Compact mode

When `compact: true`, the output uses text labels instead of section headings — suitable for small models:

- **"Equation:"** / **"Equations:"** followed by the align block and argument range
- **"initial condition:"** / **"initial conditions:"** with inline `$x(0) = 2$, $y(0) = 0$`
- **"where"** for expressions, listed inline as `$E_1 = ...$, $E_2 = ...$`
- Parameters and constants remain as tables

```typescript
const compact = convertIvpToLatex(ivp, {compact: true});
```

## Examples

The [examples/](examples/) directory contains 13 `.ivp` model files and 13 `.ts` scripts demonstrating `convertIvpToLatex` with different option combinations. Models are based on the [Datagrok model library](https://datagrok.ai/help/compute/models).

Models with up to 3 equations use `compact: true` for concise output; larger models use section headings.

### Running examples

```bash
# Run a single example
npx tsx src/latex-export/examples/basic-to-latex.ts

# Run all examples
for f in src/latex-export/examples/*.ts; do echo "=== $f ==="; npx tsx "$f"; done
```

### Example catalog

| Example | Model | Format | Options |
|---|---|---|---|
| [basic-to-latex.ts](examples/basic-to-latex.ts) | Simple 1D ODE template | LaTeX | `compact` |
| [robertson-to-latex-no-meta.ts](examples/robertson-to-latex-no-meta.ts) | Robertson's stiff chemical kinetics (3 species) | LaTeX | `compact`, `includeMetadata: false` |
| [advanced-to-latex.ts](examples/advanced-to-latex.ts) | 2D ODE template with expressions, constants, and parameters | LaTeX | `compact`, `includeMetadata: false`, `useCdot: false` |
| [fermentation-to-latex.ts](examples/fermentation-to-latex.ts) | Ethanol fermentation process kinetics | LaTeX | `compact`, `includeParameters: false` |
| [nimotuzumab-to-latex.ts](examples/nimotuzumab-to-latex.ts) | Population PK of the monoclonal antibody nimotuzumab | LaTeX | `compact`, `includeInits: false`, `useCdot: false` |
| [chem-react-to-latex-no-cdot.ts](examples/chem-react-to-latex-no-cdot.ts) | Deterministic mass-action kinetics (4 species) | LaTeX | `useCdot: false` |
| [pollution-to-latex.ts](examples/pollution-to-latex.ts) | Air pollution chemistry (20 species, 25 reactions) | LaTeX | `includeMetadata: false`, `includeInits: false` |
| [bioreactor-to-latex-selective.ts](examples/bioreactor-to-latex-selective.ts) | Fab-arm exchange bioreactor simulation (13 species) | LaTeX | `includeConstants: false`, `includeParameters: false` |
| [extended-to-markdown.ts](examples/extended-to-markdown.ts) | 2D ODE system with expressions and annotations | Markdown | `compact` |
| [pk-to-markdown.ts](examples/pk-to-markdown.ts) | One-compartment pharmacokinetic model with dosing loop | Markdown | `compact` |
| [energy-n-control-to-markdown.ts](examples/energy-n-control-to-markdown.ts) | Custom JS functions, ternary operators, arrow functions | Markdown | `compact`, `includeConstants: false` |
| [ga-production-to-markdown.ts](examples/ga-production-to-markdown.ts) | Gluconic acid production by Aspergillus niger (multi-stage) | Markdown | `useCdot: false` |
| [pk-pd-to-markdown-minimal.ts](examples/pk-pd-to-markdown-minimal.ts) | Two-compartment PK-PD model with dosing loop | Markdown | `includeInits: false`, `includeParameters: false`, `includeConstants: false` |

## What gets converted

### Equations

ODE system from `#equations` block. Derivatives become fractions:

| IVP | LaTeX |
|---|---|
| `dy/dt = -k * y` | `\frac{dy}{dt} = -k \cdot y` |
| `d(FFox)/dt = -E11 + E12` | `\frac{d\,\mathrm{FFox}}{dt} = -\mathrm{E11} + \mathrm{E12}` |
| `dy1/dt = r1 + r2` | `\frac{dy_{1}}{dt} = r_{1} + r_{2}` |

Equations are rendered inside `align` (LaTeX) or `aligned` (Markdown) environments.

### Expressions

Auxiliary formulas from `#expressions` block, including:

- Regular expressions: `E1 = C1 * exp(-t) + P1`
- Arrow functions: `func1 = (p, t) => (p < 8) ? ceil(t) : ceil(t**2 / 20)`
- Ternary operators: `control = (P1 < 3) ? cos(t) : sin(t)`
- Aliases: `ceil = Math.ceil`

### Identifiers

Automatic recognition and formatting:

| Pattern | Example | LaTeX |
|---|---|---|
| Greek letter | `alpha` | `\alpha` |
| Greek + subscript | `mu2` | `\mu_{2}` |
| Letter + digits | `x1`, `k25` | `x_{1}`, `k_{25}` |
| Chemical formula | `CO2`, `N2O5` | `\mathrm{CO}_{2}`, `\mathrm{N_{2}O_{5}}` |
| Multi-letter name | `FFox`, `MEAthiol` | `\mathrm{FFox}`, `\mathrm{MEAthiol}` |
| Special constants | `PI`, `Inf` | `\pi`, `\infty` |

### Operators

| IVP | LaTeX | Notes |
|---|---|---|
| `a * b` | `a \cdot b` | Configurable via `useCdot` |
| `a / b` | `\frac{a}{b}` | Inline `/` inside exponents |
| `x**2` | `x^{2}` | Parentheses added for compound bases |
| `>=` | `\geq` | |
| `<=` | `\leq` | |
| `!=` | `\neq` | |
| `==` | `=` | |

### Functions

| IVP | LaTeX |
|---|---|
| `sin(x)`, `cos(x)`, `tan(x)` | `\sin\!\left(x\right)`, `\cos\!\left(x\right)`, `\tan\!\left(x\right)` |
| `exp(x)` | `e^{x}` |
| `sqrt(x)` | `\sqrt{x}` |
| `pow(a, b)` | `a^{b}` (parentheses around compound base) |
| `log(x)` | `\ln\!\left(x\right)` |
| `log10(x)` | `\log_{10}\!\left(x\right)` |
| `abs(x)` | `\left\lvert x \right\rvert` |
| `ceil(x)` / `Math.ceil(x)` | `\left\lceil x \right\rceil` |
| `Math.floor(x)` | `\left\lfloor x \right\rfloor` |

### Scientific notation

| IVP | LaTeX |
|---|---|
| `1e4` | `10^{4}` |
| `3e7` | `3 \times 10^{7}` |
| `9.2E-2` | `9.2 \times 10^{-2}` |

### Ternary operator

```
Fin = t < switchTime ? 0 : 0.025
```

Rendered as a `cases` environment:

```latex
\mathrm{Fin} = \begin{cases} 0, & \text{if } t < \mathrm{switchTime} \\ 0.025, & \text{otherwise} \end{cases}
```

### Tables

Initial conditions, parameters, and constants from `#inits`, `#parameters`, and `#constants` blocks are rendered as tables:

- **LaTeX:** `tabular` with `booktabs` rules (`\toprule`, `\midrule`, `\bottomrule`)
- **Markdown:** pipe-delimited tables with inline math `$...$`

Tables include name, value, and units (when available).

## Output formats

### LaTeX

Produces a standalone `.tex`-ready document with:

- `\documentclass{article}` preamble
- `amsmath`, `amssymb`, `booktabs` packages
- `\section` / `\subsection` structure
- `align` environments for equations
- `tabular` environments for data tables

### Markdown

Produces Markdown with:

- `##` / `###` heading structure
- `$$\begin{aligned}...\end{aligned}$$` display math blocks
- `$...$` inline math in tables
- Standard Markdown pipe tables

## Supported IVP blocks

| Block | Converted | Description |
|---|---|---|
| `#name` | Yes | Model title |
| `#description` | Yes | Short description |
| `#comment` | Yes | Free-form text |
| `#equations` | Yes | ODE system (primary content) |
| `#expressions` | Yes | Auxiliary formulas |
| `#inits` | Yes | Initial conditions |
| `#parameters` | Yes | Model parameters |
| `#constants` | Yes | Fixed constants |
| `#argument` | Yes | Independent variable and range |
| `#loop` | Yes | Cyclic process (dosing) |
| `#update` | Yes | Multi-stage process |
| `#output` | Parsed | Output columns |
| `#tolerance` | Parsed | Solver tolerance |
| `#meta.solver` | Parsed | Solver settings |

## Architecture

```
index.ts                  Public API: convertIvpToLatex()
  │
  ├── parser/
  │   ├── line-joiner.ts       Strip comments, join multi-line formulas
  │   ├── ivp-parser.ts        Parse .ivp text → ParsedModel
  │   ├── tokenizer.ts         Expression string → token stream
  │   └── ast-parser.ts        Tokens → Abstract Syntax Tree
  │
  ├── transformer/
  │   ├── identifier.ts        Identifiers → Greek, subscripts, \mathrm
  │   ├── functions.ts         Function calls → LaTeX
  │   └── operators.ts         Operators → LaTeX
  │
  ├── generator/
  │   ├── latex-generator.ts   AST → LaTeX string
  │   ├── bracket-manager.ts   Precedence-based parenthesis logic
  │   └── document-builder.ts  Assemble complete document
  │
  ├── output/
  │   ├── tex-formatter.ts     Wrap in \documentclass preamble
  │   └── md-formatter.ts      Wrap math in $...$ / $$...$$
  │
  └── types.ts                 All TypeScript interfaces
```

**Data flow:**

```
.ivp text → parseIvp() → ParsedModel → buildLatexDocument() / buildMarkdownDocument() → string
                              │
                    For each formula:
                    expressionToLatex(rhs)
                        │
                  tokenize → parse → AST → nodeToLatex
                                            ├── identifierToLatex
                                            ├── functionToLatex
                                            ├── operatorToLatex
                                            └── needsParentheses
```

## Testing

```bash
npx jest src/tests/latex-export/                   # all latex-export tests
npx jest src/tests/latex-export/parser/            # parser tests only
npx jest src/tests/latex-export/transformer/       # transformer tests only
npx jest src/tests/latex-export/generator/         # generator tests only
npx jest src/tests/latex-export/integration/       # end-to-end with real .ivp files
```
