// ============================================================
// Expression Tokenizer
// ============================================================
// Converts a mathematical expression string into a flat array of tokens.
//
// Key rules:
//   - `**` is a single OPERATOR token (exponentiation)
//   - `=>` is a single ARROW token
//   - `>=`, `<=`, `==`, `!=` are single COMPARISON tokens
//   - `Math.ceil`, `Math.floor` are single IDENTIFIER tokens
//   - Scientific notation: `1e4`, `9.2E-2` are single NUMBER tokens
//   - `-` after `e`/`E` in a number context is part of the number

import { Token } from '../types';

export type { Token } from '../types';

/**
 * Tokenize a mathematical expression string.
 *
 * @param expr - Expression string (RHS of a formula)
 * @returns Array of tokens
 */
export function tokenize(expr: string): Token[] {
  // TODO: Implement
  throw new Error('Not implemented');
}
