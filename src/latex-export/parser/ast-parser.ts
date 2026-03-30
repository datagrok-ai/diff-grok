/* AST parser (recursive descent): builds an Abstract Syntax Tree from a token stream. */

import {Token, ASTNode, ERROR_MSG} from '../types';
import {tokenize} from './tokenizer';

export type {ASTNode} from '../types';

/** Parse an expression string into an AST.
 *  @param expr  expression string
 *  @returns     root AST node
 */
export function parseExpression(expr: string): ASTNode {
  const tokens = tokenize(expr);
  const parser = new Parser(tokens);
  return parser.parse();
}

class Parser {
  private tokens: Token[];
  private pos: number;

  constructor(tokens: Token[]) {
    this.tokens = tokens;
    this.pos = 0;
  }

  parse(): ASTNode {
    const node = this.parseTernary();
    return node;
  }

  private peek(): Token | undefined {
    return this.tokens[this.pos];
  }

  private advance(): Token {
    return this.tokens[this.pos++];
  }

  private expect(type: string, value?: string): Token {
    const tok = this.advance();
    if (!tok || tok.type !== type || (value !== undefined && tok.value !== value))
      throw new Error(ERROR_MSG.expectedToken(type, value, tok));
    return tok;
  }

  // 1. Ternary: condition ? consequent : alternate (right-associative)
  private parseTernary(): ASTNode {
    const node = this.parseComparison();
    const tok = this.peek();
    if (tok && tok.type === 'QUESTION') {
      this.advance(); // skip ?
      const consequent = this.parseTernary();
      this.expect('COLON');
      const alternate = this.parseTernary();
      return {type: 'ternary', condition: node, consequent, alternate};
    }
    return node;
  }

  // 2. Comparison: < > <= >= == !=
  private parseComparison(): ASTNode {
    let node = this.parseAddition();
    while (true) {
      const tok = this.peek();
      if (tok && tok.type === 'COMPARISON') {
        const op = this.advance().value;
        const right = this.parseAddition();
        node = {type: 'binary', op, left: node, right};
      } else
        break;
    }
    return node;
  }

  // 3. Addition: + -
  private parseAddition(): ASTNode {
    let node = this.parseMultiplication();
    while (true) {
      const tok = this.peek();
      if (tok && tok.type === 'OPERATOR' && (tok.value === '+' || tok.value === '-')) {
        const op = this.advance().value;
        const right = this.parseMultiplication();
        node = {type: 'binary', op, left: node, right};
      } else
        break;
    }
    return node;
  }

  // 4. Multiplication: * /
  private parseMultiplication(): ASTNode {
    let node = this.parseUnary();
    while (true) {
      const tok = this.peek();
      if (tok && tok.type === 'OPERATOR' && (tok.value === '*' || tok.value === '/')) {
        const op = this.advance().value;
        const right = this.parseUnary();
        node = {type: 'binary', op, left: node, right};
      } else
        break;
    }
    return node;
  }

  // 5. Unary: -x +x
  private parseUnary(): ASTNode {
    const tok = this.peek();
    if (tok && tok.type === 'OPERATOR' && (tok.value === '-' || tok.value === '+')) {
      const op = this.advance().value;
      const operand = this.parseUnary();
      return {type: 'unary', op, operand};
    }
    return this.parseExponentiation();
  }

  // 6. Exponentiation: ** (right-associative)
  private parseExponentiation(): ASTNode {
    const base = this.parseCallOrAtom();
    const tok = this.peek();
    if (tok && tok.type === 'OPERATOR' && tok.value === '**') {
      this.advance();
      const exp = this.parseUnary(); // right-associative: recurse through unary to handle -x**2
      return {type: 'binary', op: '**', left: base, right: exp};
    }
    return base;
  }

  // 7. Function call / 8. Atom
  private parseCallOrAtom(): ASTNode {
    const tok = this.peek();

    // Parenthesized expression
    if (tok && tok.type === 'LPAREN') {
      this.advance(); // skip (
      const node = this.parseTernary();
      this.expect('RPAREN');
      return node;
    }

    // Number
    if (tok && tok.type === 'NUMBER') {
      this.advance();
      return {type: 'number', value: tok.value};
    }

    // Identifier — possibly a function call
    if (tok && tok.type === 'IDENTIFIER') {
      this.advance();
      const next = this.peek();
      if (next && next.type === 'LPAREN') {
        // Function call
        this.advance(); // skip (
        const args: ASTNode[] = [];
        if (this.peek() && this.peek()!.type !== 'RPAREN') {
          args.push(this.parseTernary());
          while (this.peek() && this.peek()!.type === 'COMMA') {
            this.advance(); // skip ,
            args.push(this.parseTernary());
          }
        }
        this.expect('RPAREN');
        return {type: 'call', name: tok.value, args};
      }
      return {type: 'identifier', name: tok.value};
    }

    throw new Error(ERROR_MSG.unexpectedToken(tok));
  }
}
