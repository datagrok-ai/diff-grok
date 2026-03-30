// ============================================================
// LaTeX Generator
// ============================================================
// Recursively walks the AST and produces LaTeX strings.
//
// Key responsibilities:
//   - Render each AST node type to LaTeX
//   - Apply identifier, function, and operator transformations
//   - Manage brackets via bracket-manager
//   - Handle derivative LHS patterns: d(X)/dt → \frac{d\,X}{dt}
//   - Render ternary as \begin{cases} ... \end{cases}
//   - Render scientific notation: 1e4 → 1 \times 10^{4}

import { ASTNode } from '../types';

/**
 * Convert an expression string to LaTeX by parsing and rendering the AST.
 *
 * @param expr - Expression string (RHS of a formula)
 * @returns LaTeX string
 */
export function expressionToLatex(expr: string): string {
  // TODO: Implement
  // 1. Parse expression to AST
  // 2. Recursively render AST to LaTeX
  throw new Error('Not implemented');
}

/**
 * Render an AST node to LaTeX.
 *
 * @param node - AST node
 * @param parentOp - Parent operator for bracket decisions (optional)
 * @param position - 'left' or 'right' in parent (optional)
 * @returns LaTeX string
 */
export function nodeToLatex(
  node: ASTNode,
  parentOp?: string,
  position?: 'left' | 'right',
): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Convert a derivative LHS to LaTeX.
 *
 * Patterns:
 *   dy/dt       → \frac{dy}{dt}
 *   dx1/dt      → \frac{dx_{1}}{dt}
 *   d(FFox)/dt  → \frac{d\,\mathrm{FFox}}{dt}
 *   d(depot)/dt → \frac{d\,\mathrm{depot}}{dt}
 *
 * @param lhs - The LHS string (e.g., "d(FFox)/dt")
 * @returns LaTeX string for the derivative
 */
export function derivativeToLatex(lhs: string): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Render a number, potentially converting scientific notation.
 *
 * Rules:
 *   42      → 42
 *   0.04    → 0.04
 *   1e4     → 10^{4}   or  1 \times 10^{4}
 *   9.2E-2  → 9.2 \times 10^{-2}
 *   3e7     → 3 \times 10^{7}
 *
 * @param value - Number string from token
 * @returns LaTeX string
 */
export function numberToLatex(value: string): string {
  // TODO: Implement
  throw new Error('Not implemented');
}
