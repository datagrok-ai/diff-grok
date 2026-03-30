# latex-export — Implementation Specification

> Module of the `diff-grok` library. Location: `src/latex-export/`

## Overview

Build a **TypeScript** module that converts IVP model files (a declarative ODE notation used by Datagrok Diff Studio) into **LaTeX** markup suitable for inclusion in `.tex` documents and Markdown files.

The converter receives the full text of an `.ivp` file and produces a LaTeX string.

## Input format

IVP files are plain-text files composed of **blocks** that start with `#keyword`. Full syntax reference is in `docs/latex-export/SYNTAX.md`. Example files are in the `examples/` directory.

Key blocks relevant to LaTeX conversion:

| Block | Purpose |
|---|---|
| `#name` | Model title |
| `#description` | Short description |
| `#comment` | Free-form comment (can be multi-line) |
| `#equations` | System of ODEs (the primary content) |
| `#expressions` | Auxiliary computed formulas |
| `#inits` | Initial conditions with optional `{annotations}` and `[tooltips]` |
| `#parameters` | Model parameters with optional annotations |
| `#constants` | Fixed constants |
| `#argument` | Independent variable and its range |
| `#loop` | Cyclic process definition |
| `#update` | Multi-stage process definition |
| `#output` | Output columns configuration |
| `#tolerance` | Solver tolerance |
| `#meta.solver` | Solver settings |
| `#meta.inputs` | Lookup table reference |

## Output format

The converter must produce **two output formats**, controlled by an option:

1. **LaTeX** — standalone content for a `.tex` file (using `amsmath` environments: `align`, `cases`, `tabular`)
2. **Markdown** — LaTeX formulas wrapped in `$...$` and `$$...$$` delimiters

## Architecture

```
src/latex-export/
├── parser/
│   ├── ivp-parser.ts         // Parse IVP text → structured blocks
│   ├── line-joiner.ts        // Strip comments, join multi-line formulas
│   ├── tokenizer.ts          // Expression text → token stream
│   └── ast-parser.ts         // Token stream → AST
├── transformer/
│   ├── identifier.ts         // Identifier → LaTeX (Greek, subscripts, mathrm)
│   ├── functions.ts          // Function calls → LaTeX
│   └── operators.ts          // Operator precedence and rendering rules
├── generator/
│   ├── latex-generator.ts    // AST → LaTeX string
│   ├── bracket-manager.ts    // Parentheses logic
│   └── document-builder.ts   // Assemble full document from all blocks
├── output/
│   ├── tex-formatter.ts      // Format for .tex
│   └── md-formatter.ts       // Format for .md
├── types.ts                  // All TypeScript types/interfaces
└── index.ts                  // Public API: convertIvpToLatex(text, options)
```

Tests mirror this structure under `src/tests/latex-export/`.

## Public API

```typescript
interface ConvertOptions {
  /** Output format: 'latex' or 'markdown'. Default: 'latex' */
  format: 'latex' | 'markdown';

  /** Include metadata sections (name, description, comment). Default: true */
  includeMetadata: boolean;

  /** Include initial conditions table. Default: true */
  includeInits: boolean;

  /** Include parameters table. Default: true */
  includeParameters: boolean;

  /** Include constants table. Default: true */
  includeConstants: boolean;

  /** Use \cdot for multiplication. Default: true. If false, use juxtaposition. */
  useCdot: boolean;

  /** Compact mode: no section headings, inline initial conditions. Default: false */
  compact: boolean;
}

function convertIvpToLatex(ivpText: string, options?: Partial<ConvertOptions>): string;
```

---

## Detailed requirements

### 1. IVP parser (`parser/ivp-parser.ts`)

Parse the IVP text into structured blocks.

**Input:** raw string content of an `.ivp` file.

**Output:** a `ParsedModel` object:

```typescript
interface ParsedModel {
  name?: string;
  description?: string;
  comment?: string;
  equations: FormulaLine[];
  expressions: FormulaLine[];
  inits: AnnotatedLine[];
  parameters: AnnotatedLine[];
  constants: AnnotatedLine[];
  argument: { name: string; entries: AnnotatedLine[] };
  loops: LoopBlock[];
  updates: UpdateBlock[];
  output: OutputLine[];
  tolerance?: string;
  solverMeta?: string;
}

interface FormulaLine {
  lhs: string;        // e.g. "dx/dt", "E1", "func1"
  rhs: string;        // full right-hand side expression
  isDerivative: boolean;
  isArrowFunction: boolean;
  arrowParams?: string[];  // e.g. ["p", "t"]
}

interface AnnotatedLine {
  name: string;
  value: string;
  units?: string;
  caption?: string;
  category?: string;
  tooltip?: string;
  min?: string;
  max?: string;
  step?: string;
}
```

### 2. Line preprocessing (`parser/line-joiner.ts`)

Process raw lines before tokenization:

1. **Strip comments:** remove everything from `//` to end of line. Note: `//` is always a comment (IVP has no string literals). Do NOT confuse with the division operator `/` — a single slash is division, only double-slash `//` is a comment.

2. **Strip annotations:** remove `{...}` blocks and `[...]` tooltip blocks from lines in `#inits`, `#parameters`, `#argument` (but parse them first for metadata).

3. **Join multi-line formulas:** if a line is indented and does NOT contain `=` at the top level (outside parentheses), it is a continuation of the previous formula. Concatenate with a space.

4. **Handle block headers:** lines starting with `#keyword` begin a new block.

### 3. Tokenizer (`parser/tokenizer.ts`)

Convert an expression string into a flat array of tokens.

**Token types:**

```typescript
type TokenType =
  | 'NUMBER'       // 0.04, 1e4, 1E-2, 3e7, 9.2E-2
  | 'IDENTIFIER'   // x1, FFox, MEAthiol, pKa2MEA, Math.ceil
  | 'OPERATOR'     // +  -  *  /  **
  | 'LPAREN'       // (
  | 'RPAREN'       // )
  | 'COMPARISON'   // >=  <=  >  <  ==  !=
  | 'QUESTION'     // ?
  | 'COLON'        // :
  | 'COMMA'        // ,
  | 'ARROW'        // =>
  | 'ASSIGN'       // =  +=
  | 'SEMICOLON';   // ;

interface Token {
  type: TokenType;
  value: string;
  position: number;
}
```

**Rules:**
- `**` is a single OPERATOR token (not two `*`).
- `=>` is a single ARROW token.
- `>=`, `<=`, `==`, `!=` are single COMPARISON tokens.
- `Math.ceil` should be tokenized as a single IDENTIFIER (detect `Math.` prefix).
- Scientific notation: `1.23e4`, `9.2E-2`, `1E4` — single NUMBER token.
- Negative exponents in scientific notation: `1E-2` — the `-` is part of the number. Heuristic: if `-` follows `e`/`E` immediately in a number context, it's part of the number.

### 4. AST parser (`parser/ast-parser.ts`)

Build an Abstract Syntax Tree from tokens using **recursive descent** with the following precedence (lowest to highest):

```
1. Ternary           ? :          (right-associative)
2. Comparison        < > <= >= == !=
3. Addition          + -
4. Multiplication    * /
5. Unary             -x  +x
6. Exponentiation    **           (right-associative)
7. Function call     f(...)
8. Atom              number, identifier, (group)
```

**AST node types:**

```typescript
type ASTNode =
  | NumberNode
  | IdentifierNode
  | BinaryNode
  | UnaryNode
  | CallNode
  | TernaryNode;

interface NumberNode { type: 'number'; value: string; }
interface IdentifierNode { type: 'identifier'; name: string; }
interface BinaryNode { type: 'binary'; op: string; left: ASTNode; right: ASTNode; }
interface UnaryNode { type: 'unary'; op: string; operand: ASTNode; }
interface CallNode { type: 'call'; name: string; args: ASTNode[]; }
interface TernaryNode { type: 'ternary'; condition: ASTNode; consequent: ASTNode; alternate: ASTNode; }
```

### 5. Identifier transformation (`transformer/identifier.ts`)

Convert identifiers to LaTeX representation.

**Step 1: Check Greek alphabet**

Full Greek alphabet mapping is in `docs/latex-export/GREEK.md`. Check if the identifier matches a Greek letter name. Match must be **greedy by length** (check longer names first to avoid `varepsilon` matching as `var` + `epsilon`).

Rules:
- **Exact match:** `alpha` → `\alpha`, `Gamma` → `\Gamma`
- **Greek prefix + digit suffix:** `alpha1` → `\alpha_{1}`, `mu2` → `\mu_{2}`
- **Do NOT match** if the remaining suffix starts with a letter: `muM` → NOT `\mu M`, treat as plain identifier `\mathrm{muM}`

**Step 2: Subscript extraction**

If the identifier is NOT Greek:
- Single letter + digits: `x1` → `x_{1}`, `y20` → `y_{20}`, `k25` → `k_{25}`
- Single letter + single digit: `V2` → `V_{2}`, `C3` → `C_{3}`

**Step 3: Multi-letter identifiers**

- If more than one letter and not Greek: wrap in `\mathrm{}`: `FFox` → `\mathrm{FFox}`, `MEAthiol` → `\mathrm{MEAthiol}`
- Chemical formulas heuristic: uppercase letter(s) followed by digits → separate: `CO2` → `\mathrm{CO}_{2}`, `N2O5` → `\mathrm{N_{2}O_{5}}`

**Special constants:**
- `PI` → `\pi`
- `Inf`, `inf` → `\infty`

### 6. Function transformation (`transformer/functions.ts`)

Map function calls to LaTeX:

| IVP function | LaTeX |
|---|---|
| `sin(x)` | `\sin\!\left(x\right)` |
| `cos(x)` | `\cos\!\left(x\right)` |
| `tan(x)` | `\tan\!\left(x\right)` |
| `exp(x)` | `e^{x}` — if argument is simple; `\exp\!\left(x\right)` if complex |
| `sqrt(x)` | `\sqrt{x}` |
| `pow(a, b)` | `a^{b}` — with brackets around `a` if compound |
| `log(x)` | `\ln(x)` |
| `log10(x)` | `\log_{10}(x)` |
| `abs(x)` | `\left\lvert x \right\rvert` |
| `Math.ceil(x)` | `\left\lceil x \right\rceil` |
| `ceil(x)` | `\left\lceil x \right\rceil` (alias) |
| `Math.floor(x)` | `\left\lfloor x \right\rfloor` |

**`exp()` heuristic:** if the argument is a single token or a simple negation (e.g., `-t`), render as `e^{-t}`. If the argument is complex, render as `e^{...}` with the full expression in the exponent.

### 7. Operator transformation (`transformer/operators.ts`)

**Division → fraction:**
- Default: render `/` as `\frac{numerator}{denominator}`
- Exception: if inside an exponent or subscript, use inline `/`

**Multiplication:**
- `*` between two identifiers/calls: `\cdot` (or juxtaposition if `useCdot: false`)
- `*` between number and identifier: `2 \cdot x` if `useCdot: true`
- `*` between number and number: `\cdot`

**Exponentiation:**
- `**` or `pow()`: render as `base^{exponent}`
- If base is compound, wrap in parentheses: `(x + y)**2` → `\left(x + y\right)^{2}`

### 8. Bracket management (`generator/bracket-manager.ts`)

Determine when parentheses are needed based on parent-child precedence.

**Rules:**
- Child has lower precedence than parent → add `\left(` ... `\right)`
- Inside `\frac{}{}` → no outer brackets needed
- Inside `\sqrt{}` → no outer brackets

### 9. LaTeX generator (`generator/latex-generator.ts`)

Recursively walk the AST and produce a LaTeX string.

**Derivative LHS patterns:**
- `dx/dt` → `\frac{dx}{dt}`
- `d(FFox)/dt` → `\frac{d\,\mathrm{FFox}}{dt}`
- `dy1/dt` → `\frac{dy_{1}}{dt}`

**Ternary → `cases`:**
```latex
\begin{cases}
  consequent, & \text{if } condition \\
  alternate, & \text{otherwise}
\end{cases}
```

**Arrow functions:**
`func1 = (p, t) => expr` → `\mathrm{func1}(p, t) = expr_latex`

### 10. Document builder (`generator/document-builder.ts`)

Assemble all converted blocks into a complete document.

**LaTeX output:** `\section`, `\subsection`, `\begin{align}`, `\begin{tabular}` with `\toprule/\midrule/\bottomrule`.

**Markdown output:** `##`, `###`, `$$\begin{aligned}...\end{aligned}$$`, markdown tables with `$...$` inline math.

---

## Implementation order

| Step | Module | Test file |
|------|--------|-----------|
| 1 | `types.ts` | — |
| 2 | `parser/line-joiner.ts` | `tests/latex-export/parser/line-joiner.test.ts` |
| 3 | `parser/ivp-parser.ts` | `tests/latex-export/parser/ivp-parser.test.ts` |
| 4 | `parser/tokenizer.ts` | `tests/latex-export/parser/tokenizer.test.ts` |
| 5 | `parser/ast-parser.ts` | `tests/latex-export/parser/ast-parser.test.ts` |
| 6 | `transformer/identifier.ts` | `tests/latex-export/transformer/identifier.test.ts` |
| 7 | `transformer/functions.ts` | `tests/latex-export/transformer/functions.test.ts` |
| 8 | `transformer/operators.ts` | `tests/latex-export/transformer/operators.test.ts` |
| 9 | `generator/bracket-manager.ts` | — |
| 10 | `generator/latex-generator.ts` | `tests/latex-export/generator/latex-generator.test.ts` |
| 11 | `generator/document-builder.ts` | `tests/latex-export/generator/document-builder.test.ts` |
| 12 | `output/tex-formatter.ts`, `md-formatter.ts` | — |
| 13 | `index.ts` | `tests/latex-export/integration/ivp-files.test.ts` |

## Testing

```bash
npx vitest run                                             # all tests
npx vitest run src/tests/latex-export/parser/              # parser only
npx vitest run src/tests/latex-export/integration/         # integration only
```
