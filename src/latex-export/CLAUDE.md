# CLAUDE.md — diff-grok / latex-export

## Project context

`diff-grok` is a library for working with Diff Studio ODE models. The `latex-export` module converts `.ivp` model files into LaTeX and Markdown with LaTeX math.

## Repository layout

```
diff-grok/
├── CLAUDE.md                          ← you are here
├── SPEC.md                            ← full implementation spec (read first!)
├── package.json
├── tsconfig.json
├── vitest.config.ts
├── docs/
│   └── latex-export/
│       ├── GREEK.md                   ← Greek alphabet mapping
│       └── SYNTAX.md                  ← IVP file format reference
├── examples/
│   └── *.ivp                          ← 12 real model files for testing
└── src/
    ├── latex-export/                  ← MODULE SOURCE CODE
    │   ├── index.ts                   ← public API
    │   ├── types.ts                   ← all interfaces
    │   ├── parser/
    │   │   ├── line-joiner.ts
    │   │   ├── ivp-parser.ts
    │   │   ├── tokenizer.ts
    │   │   └── ast-parser.ts
    │   ├── transformer/
    │   │   ├── identifier.ts
    │   │   ├── functions.ts
    │   │   └── operators.ts
    │   ├── generator/
    │   │   ├── latex-generator.ts
    │   │   ├── bracket-manager.ts
    │   │   └── document-builder.ts
    │   └── output/
    │       ├── tex-formatter.ts
    │       └── md-formatter.ts
    └── tests/
        └── latex-export/              ← MODULE TESTS
            ├── parser/
            │   ├── line-joiner.test.ts
            │   ├── ivp-parser.test.ts
            │   ├── tokenizer.test.ts
            │   └── ast-parser.test.ts
            ├── transformer/
            │   ├── identifier.test.ts
            │   ├── functions.test.ts
            │   └── operators.test.ts
            ├── generator/
            │   ├── latex-generator.test.ts
            │   └── document-builder.test.ts
            └── integration/
                └── ivp-files.test.ts
```

## Key documents

- **`SPEC.md`** — Complete implementation specification. **Read first.**
- **`docs/latex-export/SYNTAX.md`** — IVP file format reference.
- **`docs/latex-export/GREEK.md`** — Full Greek alphabet mapping with matching rules.
- **`examples/*.ivp`** — 12 real model files to test against.

## Implementation order

Implement modules in this order. Each step builds on the previous.

1. **`src/latex-export/parser/line-joiner.ts`** → run `src/tests/latex-export/parser/line-joiner.test.ts`
2. **`src/latex-export/parser/ivp-parser.ts`** → run `src/tests/latex-export/parser/ivp-parser.test.ts`
3. **`src/latex-export/parser/tokenizer.ts`** → run `src/tests/latex-export/parser/tokenizer.test.ts`
4. **`src/latex-export/parser/ast-parser.ts`** → run `src/tests/latex-export/parser/ast-parser.test.ts`
5. **`src/latex-export/transformer/identifier.ts`** → run `src/tests/latex-export/transformer/identifier.test.ts`
6. **`src/latex-export/transformer/functions.ts`** → run `src/tests/latex-export/transformer/functions.test.ts`
7. **`src/latex-export/transformer/operators.ts`** → run `src/tests/latex-export/transformer/operators.test.ts`
8. **`src/latex-export/generator/bracket-manager.ts`** — used by latex-generator
9. **`src/latex-export/generator/latex-generator.ts`** → run `src/tests/latex-export/generator/latex-generator.test.ts`
10. **`src/latex-export/generator/document-builder.ts`** → run `src/tests/latex-export/generator/document-builder.test.ts`
11. **`src/latex-export/output/tex-formatter.ts`** and **`md-formatter.ts`** — partially done
12. **`src/latex-export/index.ts`** — wire everything. Run `src/tests/latex-export/integration/ivp-files.test.ts`

## Critical conversion rules

### Comments
`//` is a single-line comment (like TypeScript). Everything after `//` to end of line is ignored. Do NOT confuse with `/` (division). Only `//` is a comment.

### Multi-line formulas
Continuation lines are indented and do NOT contain a top-level `=`. Join them with a space.

### Derivatives
- `dy/dt` → `\frac{dy}{dt}`
- `d(FFox)/dt` → `\frac{d\,\mathrm{FFox}}{dt}`
- `dy1/dt` → `\frac{dy_{1}}{dt}`

### Identifiers
- Greek letters: exact match or prefix+digit. See `docs/latex-export/GREEK.md`. Greedy by length.
- `muM` is NOT `\mu M` — suffix starts with letter → treat as plain identifier.
- Single letter + digits: `x1` → `x_{1}`
- Multi-letter: `FFox` → `\mathrm{FFox}`
- Chemical: `CO2` → `\mathrm{CO}_{2}`, `N2O5` → `\mathrm{N_{2}O_{5}}`
- Special: `PI` → `\pi`, `Inf` → `\infty`

### Operators
- `/` → `\frac{}{}`
- `*` → `\cdot` (configurable)
- `**` → `^{}`
- `>=` → `\geq`, `<=` → `\leq`, `!=` → `\neq`, `==` → `=`

### Functions
- `sin`, `cos`, `tan` → `\sin`, `\cos`, `\tan` with `\left( \right)`
- `exp(x)` → `e^{x}` (simple arg) or `e^{...}` (complex)
- `pow(a, b)` → `a^{b}` (parens around compound `a`)
- `sqrt(x)` → `\sqrt{x}`
- `Math.ceil(x)` / `ceil(x)` → `\lceil x \rceil`

### Ternary
`a ? b : c` → `\begin{cases} b, & \text{if } a \\ c, & \text{otherwise} \end{cases}`

### Arrow functions
`func1 = (p, t) => expr` → `\mathrm{func1}(p, t) = expr_latex`

### Scientific notation
`1e4` → `10^{4}`, `9.2E-2` → `9.2 \times 10^{-2}`, `3e7` → `3 \times 10^{7}`

## Running tests

```bash
npm install
npx vitest run                                       # all tests
npx vitest run src/tests/latex-export/parser/         # parser tests only
npx vitest run src/tests/latex-export/integration/    # integration tests
npx vitest run --reporter=verbose                     # detailed output
```

## Quality bar

- All test files must pass.
- All 12 `.ivp` example files must convert without errors.
- Generated LaTeX must be valid (compilable with amsmath).
- Brackets must be correctly balanced.
- No `//` comments in output.
