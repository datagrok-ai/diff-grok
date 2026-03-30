// ============================================================
// Function Call → LaTeX Transformer
// ============================================================
// Maps math function calls to their LaTeX equivalents.
//
// Mapping table:
//   sin(x)       → \sin\!\left(x\right)
//   cos(x)       → \cos\!\left(x\right)
//   tan(x)       → \tan\!\left(x\right)
//   exp(x)       → e^{x}  (simple arg) or \exp\!\left(x\right) (complex)
//   sqrt(x)      → \sqrt{x}
//   pow(a, b)    → a^{b}  (with brackets if a is compound)
//   log(x)       → \ln(x)
//   log10(x)     → \log_{10}(x)
//   abs(x)       → \left\lvert x \right\rvert
//   Math.ceil(x) → \left\lceil x \right\rceil
//   ceil(x)      → \left\lceil x \right\rceil  (alias)
//   Math.floor(x)→ \left\lfloor x \right\rfloor
//   unknown(x)   → \mathrm{unknown}\!\left(x\right)

/**
 * Convert a function call to LaTeX.
 *
 * @param name - Function name (e.g., 'sin', 'pow', 'Math.ceil')
 * @param args - Array of already-rendered LaTeX argument strings
 * @returns LaTeX string for the function call
 */
export function functionToLatex(name: string, args: string[]): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Determine if an argument string is "simple" (single token, no binary ops).
 * Used by exp() to decide between e^{x} and \exp(...) notation.
 */
export function isSimpleArg(arg: string): boolean {
  // TODO: Implement
  throw new Error('Not implemented');
}
