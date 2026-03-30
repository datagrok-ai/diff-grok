# latex-export — Test Plan

## Test structure

```
src/tests/latex-export/
├── parser/
│   ├── line-joiner.test.ts      # Comment stripping, multi-line joining
│   ├── ivp-parser.test.ts       # Block parsing
│   ├── tokenizer.test.ts        # Tokenization
│   └── ast-parser.test.ts       # AST construction
├── transformer/
│   ├── identifier.test.ts       # Greek letters, subscripts, mathrm
│   ├── functions.test.ts        # Math function → LaTeX
│   └── operators.test.ts        # Division, multiplication, power
├── generator/
│   ├── latex-generator.test.ts  # AST → LaTeX string
│   └── document-builder.test.ts # Full document assembly
└── integration/
    └── ivp-files.test.ts        # End-to-end tests with real .ivp files
```

Source modules under test: `src/latex-export/`

## Running tests

```bash
npx vitest run                                           # all project tests
npx vitest run src/tests/latex-export/                   # all latex-export tests
npx vitest run src/tests/latex-export/parser/            # parser tests only
npx vitest run src/tests/latex-export/transformer/       # transformer tests only
npx vitest run src/tests/latex-export/generator/         # generator tests only
npx vitest run src/tests/latex-export/integration/       # integration tests only
npx vitest run --reporter=verbose                        # detailed output
npm run test:latex                                       # shortcut (see package.json)
```

## Coverage targets

Each module should have >90% line coverage. Key areas:

1. **Comment handling:** `//` at various positions, not confused with `/`
2. **Multi-line formulas:** continuation detection and joining
3. **Scientific notation:** `1e4`, `9.2E-2`, `1E-2`, `3e7`
4. **Operator precedence:** nested expressions, right-associative `**`
5. **Greek letters:** exact match, prefix+digit, no false positives
6. **Subscripts:** single letter+digits, multi-letter names, chemical formulas
7. **Fractions:** simple division, nested fractions, division in exponents
8. **Ternary → cases:** simple and nested ternary operators
9. **Arrow functions:** parameter extraction, body parsing
10. **Derivatives:** `dx/dt`, `d(name)/dt` patterns
