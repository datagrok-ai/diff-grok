# Bug: `useCdot: false` is ignored by `convertIvpToLatex`

**Package:** `diff-grok`
**Version:** 1.3.0 (latest on npm at the time of filing)
**Repo:** `datagrok-ai/DiffGrok`
**Affected entry point:** `convertIvpToLatex(ivpText, options)`
**Affected option:** `ConvertOptions.useCdot`

## Summary

Setting `useCdot: false` has no effect on the rendered output. Every `*` in the source expressions is still rendered as `\cdot`. The option is documented in `src/latex-export/types.ts` as:

```ts
/** Use \cdot for multiplication. Default: true. If false, use juxtaposition. */
useCdot: boolean;
```

and the low-level operator renderer in `src/latex-export/transformer/operators.ts` **does** honour the context flag (see below). The bug lives in the intermediate AST-walking layer, which never threads the user option down to the point of use and instead hardcodes `useCdot: true` at three call sites.

## Repro

```ts
import {convertIvpToLatex} from 'diff-grok';

const MODEL = `
#name: tiny
#equations:
  dx/dt = 2 * a * b
#argument: t
  initial = 0
  final = 1
  step = 0.1
#inits:
  x = 0
`;

// Default → should render '2 \cdot a \cdot b'
console.log(convertIvpToLatex(MODEL));

// Request juxtaposition → should render '2 \, a \, b' (thin-space),
// matching what operators.ts emits for useCdot: false.
console.log(convertIvpToLatex(MODEL, {useCdot: false}));
```

### Expected

The second call emits

```
... 2 \, a \, b ...
```

### Actual

Both calls emit identical output:

```
... 2 \cdot a \cdot b ...
```

This is trivially observable by diffing the outputs of `examples/chem-react-to-latex-no-cdot.ts` against the same model with default options — the files should differ but don't.

## Root cause

The option is correctly typed in `types.ts` and merged into `opts` by `convertIvpToLatex` in `src/latex-export/index.ts`, but from there it is only plumbed into the *document-level* decisions (include sections, format). It is **never** passed into the expression renderer `expressionToLatex` / `nodeToLatex`, and the three call sites that emit `\cdot` hardcode the string directly.

### Call site 1 — `renderBinary` (multiplication path)

`src/latex-export/generator/latex-generator.ts`, inside `renderBinary`, for `op === '*'`:

```ts
if (op === '*')
    return operatorToLatex('*', l, r, { useCdot: true });   // ← hardcoded
```

`operatorToLatex` in `transformer/operators.ts` handles `useCdot: false` correctly — it returns `${left} \, ${right}` when the flag is off — but `renderBinary` always passes `true` regardless of user options.

### Call site 2 — `renderProductTerm` (multiple factors joined)

`src/latex-export/generator/latex-generator.ts`, end of `renderProductTerm`:

```ts
return nodes.map((n) => { ... }).join(' \\cdot ');   // ← hardcoded literal
```

This path is hit when `renderProductOfFractions` assembles numerator/denominator groups with more than one factor. The separator is a literal string, not computed from options.

### Call site 3 — `renderProductOfFractions` (groups joined)

`src/latex-export/generator/latex-generator.ts`, end of `renderProductOfFractions`:

```ts
return parts.join(' \\cdot ');   // ← hardcoded literal
```

Same issue — hit when a `*`/`/` chain produces multiple fraction groups.

## Proposed fix

Thread `useCdot` through the expression renderer. The minimal invasive change:

1. **Extend the public entry point** in `src/latex-export/generator/latex-generator.ts` to accept an options object:

   ```ts
   export interface RenderOptions {
     useCdot: boolean;
   }

   export function expressionToLatex(expr: string, opts: RenderOptions = {useCdot: true}): string {
     const ast = parseExpression(expr);
     return nodeToLatex(ast, opts);
   }
   ```

2. **Thread `opts`** through `nodeToLatex`, `renderBinary`, `renderProductTerm`, `renderProductOfFractions`. Each function takes `opts` as the last parameter.

3. **Replace the three hardcoded sites**:

   - `renderBinary`, `op === '*'` branch:
     ```ts
     if (op === '*')
       return operatorToLatex('*', l, r, { useCdot: opts.useCdot });
     ```
   - `renderProductTerm`:
     ```ts
     const sep = opts.useCdot ? ' \\cdot ' : ' \\, ';
     return nodes.map(/* ... */).join(sep);
     ```
   - `renderProductOfFractions`:
     ```ts
     const sep = opts.useCdot ? ' \\cdot ' : ' \\, ';
     return parts.join(sep);
     ```

4. **Update call sites in `document-builder.ts`** that invoke `expressionToLatex` to pass `{useCdot: opts.useCdot}`. There are several — every place that renders an equation, expression, stage transition, loop update, or arg-range expression.

5. **Add regression tests** (Jest) asserting that for a model containing `a*b`:
   - default → output contains `a \cdot b` and does NOT contain `a \, b`;
   - `useCdot: false` → output contains `a \, b` and does NOT contain `a \cdot b`.

   Also cover the product-of-fractions path: `a*b/c*d` or similar.

## Scope of impact

- All users of `convertIvpToLatex({useCdot: false})` get silently wrong output.
- `examples/chem-react-to-latex-no-cdot.ts`, `advanced-to-latex.ts`, `nimotuzumab-to-latex.ts`, `ga-production-to-markdown.ts` — all four example files that set `useCdot: false` produce output identical to the default, contradicting their comments (`"juxtaposition instead of \cdot"`).
- Downstream: the Diff Studio UI exposes this option as a user-facing toggle in the Markdown/LaTeX export dialog. Until the library ships the fix, the dialog carries a local post-process workaround (`replace(/ \\cdot /g, ' \\, ')`) that can be removed when `diff-grok >= 1.3.1` is released.

## Notes

- `derivativeToLatex` uses `\,` correctly for the `d` operator — no change needed there.
- `numberToLatex` emits `\times` for scientific notation. Outside the scope of this issue, but worth a follow-up: should `\times` also have a juxtaposition mode for consistency?
- Consider adding a `useCdot: false` Jest snapshot for the full pollution example, to guard against regressions in the product-of-fractions path on non-trivial models.
