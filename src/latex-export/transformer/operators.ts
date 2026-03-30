// ============================================================
// Operator → LaTeX Transformer
// ============================================================
// Renders binary and comparison operators as LaTeX.
//
// Rules:
//   +, -         → as-is
//   *            → \cdot (or juxtaposition with \,)
//   /            → \frac{num}{den}  (inline / if inside exponent)
//   **           → base^{exp}  (parens around compound base)
//   <, >         → as-is
//   <=           → \leq
//   >=           → \geq
//   ==           → =  (in math context)
//   !=           → \neq

import { OperatorContext } from '../types';

/**
 * Render a binary operator with its operands as LaTeX.
 *
 * @param op - Operator string: +, -, *, /, **, <, >, <=, >=, ==, !=
 * @param left - Already-rendered LaTeX for the left operand
 * @param right - Already-rendered LaTeX for the right operand
 * @param context - Optional context flags
 * @returns LaTeX string
 */
export function operatorToLatex(
  op: string,
  left: string,
  right: string,
  context?: OperatorContext,
): string {
  // TODO: Implement
  throw new Error('Not implemented');
}
