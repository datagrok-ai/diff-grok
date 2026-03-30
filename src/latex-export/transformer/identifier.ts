// ============================================================
// Identifier → LaTeX Transformer
// ============================================================
// Converts identifier names to their LaTeX representations.
//
// Processing order:
//   1. Check special constants (PI, Inf)
//   2. Check Greek alphabet (exact match, then prefix + digit suffix)
//   3. Single letter + digits → subscript (x1 → x_{1})
//   4. Multi-letter + trailing digits → mathrm with subscript
//   5. Chemical formula heuristic (CO2, N2O5)
//   6. Plain multi-letter → \mathrm{name}
//
// See docs/GREEK.md for the complete Greek alphabet mapping.

/**
 * Convert an identifier to its LaTeX representation.
 *
 * @param name - Identifier string from the expression
 * @returns LaTeX string
 */
export function identifierToLatex(name: string): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Check if an identifier is a known Greek letter (exact match).
 * Returns the LaTeX command or undefined.
 */
export function matchGreekExact(name: string): string | undefined {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Check if an identifier starts with a Greek letter name followed by digits.
 * Returns { latex: string, suffix: string } or undefined.
 */
export function matchGreekPrefix(name: string): { latex: string; suffix: string } | undefined {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Detect and format chemical formula patterns like CO2, N2O5, SO4.
 * Returns formatted LaTeX or undefined if not a chemical formula.
 */
export function matchChemicalFormula(name: string): string | undefined {
  // TODO: Implement
  throw new Error('Not implemented');
}
