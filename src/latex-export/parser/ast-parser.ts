// ============================================================
// AST Parser (Recursive Descent)
// ============================================================
// Builds an Abstract Syntax Tree from a token stream.
//
// Operator precedence (lowest to highest):
//   1. Ternary           ? :          (right-associative)
//   2. Comparison        < > <= >= == !=
//   3. Addition          + -
//   4. Multiplication    * /
//   5. Unary             -x  +x
//   6. Exponentiation    **           (right-associative)
//   7. Function call     f(...)
//   8. Atom              number, identifier, (group)

import { ASTNode } from '../types';

export type { ASTNode } from '../types';

/**
 * Parse an expression string into an AST.
 * Internally tokenizes the input first, then applies recursive descent.
 *
 * @param expr - Expression string
 * @returns Root AST node
 */
export function parseExpression(expr: string): ASTNode {
  // TODO: Implement
  // 1. Tokenize the expression
  // 2. Apply recursive descent parsing
  throw new Error('Not implemented');
}

/**
 * Parse a token array into an AST (for internal use / testing).
 *
 * @param tokens - Token array from the tokenizer
 * @returns Root AST node
 */
// export function parseTokens(tokens: Token[]): ASTNode { ... }
