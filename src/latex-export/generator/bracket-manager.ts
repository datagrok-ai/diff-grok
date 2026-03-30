// ============================================================
// Bracket Manager
// ============================================================
// Determines when parentheses are needed based on AST context.
//
// Rules:
//   - Child has lower precedence than parent → wrap in \left( \right)
//   - Inside \frac{}{} → no outer brackets (fraction bar groups)
//   - Inside \sqrt{} → no outer brackets
//   - Unary minus before multiplication: no extra brackets
//   - Power of grouped expression: keep original parens

/**
 * Get the precedence level of an operator.
 * Higher number = higher precedence = binds tighter.
 */
export function getPrecedence(op: string): number {
  // TODO: Implement
  // Suggested levels:
  //   ternary: 1
  //   comparison: 2
  //   +, -: 3
  //   *, /: 4
  //   unary: 5
  //   **: 6
  throw new Error('Not implemented');
}

/**
 * Determine if a child expression needs parentheses
 * given the parent operator context.
 *
 * @param childOp - The operator of the child node (or undefined for atoms)
 * @param parentOp - The operator of the parent node
 * @param position - 'left' or 'right' child
 * @returns true if parentheses are needed
 */
export function needsParentheses(
  childOp: string | undefined,
  parentOp: string,
  position: 'left' | 'right',
): boolean {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Wrap a LaTeX string in \left( \right) if needed.
 */
export function wrapParens(latex: string, needed: boolean): string {
  if (needed) return `\\left(${latex}\\right)`;
  return latex;
}
