/* Bracket manager: determines when parentheses are needed based on operator precedence. */

/** Get the precedence level of an operator.
 *  Higher number = higher precedence = binds tighter.
 *  @param op  operator string
 *  @returns   numeric precedence level (0–6)
 */
export function getPrecedence(op: string): number {
  switch (op) {
  case '?': return 1;
  case '<': case '>': case '<=': case '>=': case '==': case '!=':
    return 2;
  case '+': case '-':
    return 3;
  case '*': case '/':
    return 4;
  case 'unary': return 5;
  case '**': return 6;
  default: return 0;
  }
}

/** Determine if a child expression needs parentheses given the parent operator context.
 *  @param childOp   operator of the child node (or undefined for atoms)
 *  @param parentOp  operator of the parent node
 *  @param position  'left' or 'right' child
 *  @returns         true if parentheses are needed
 */
export function needsParentheses(
  childOp: string | undefined,
  parentOp: string,
  position: 'left' | 'right',
): boolean {
  if (childOp === undefined) return false;

  const childPrec = getPrecedence(childOp);
  const parentPrec = getPrecedence(parentOp);

  if (childPrec < parentPrec) return true;
  if (childPrec > parentPrec) return false;

  // Same precedence: right child of non-commutative ops needs parens
  if (position === 'right' && (parentOp === '-' || parentOp === '/'))
    return true;

  return false;
}

/** Wrap a LaTeX string in \left( \right) if needed.
 *  @param latex   LaTeX string to wrap
 *  @param needed  whether wrapping is required
 *  @returns       wrapped or original string
 */
export function wrapParens(latex: string, needed: boolean): string {
  if (needed) return `\\left(${latex}\\right)`;
  return latex;
}
