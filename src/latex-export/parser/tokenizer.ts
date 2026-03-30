/* Expression tokenizer: converts a mathematical expression string into a flat array of tokens. */

import {Token} from '../types';

export type {Token} from '../types';

/** Tokenize a mathematical expression string.
 *  @param expr  expression string (RHS of a formula)
 *  @returns     array of tokens
 */
export function tokenize(expr: string): Token[] {
  const tokens: Token[] = [];
  let i = 0;
  const len = expr.length;

  while (i < len) {
    const ch = expr[i];

    // Skip whitespace
    if (ch === ' ' || ch === '\t' || ch === '\r' || ch === '\n') {
      i++;
      continue;
    }

    // Numbers: digits or dot followed by digit
    if (isDigit(ch) || (ch === '.' && i + 1 < len && isDigit(expr[i + 1]))) {
      const start = i;
      // Integer part
      while (i < len && isDigit(expr[i]))
        i++;
      // Decimal part
      if (i < len && expr[i] === '.' && i + 1 < len && isDigit(expr[i + 1])) {
        i++; // skip dot
        while (i < len && isDigit(expr[i]))
          i++;
      }
      // Scientific notation: e/E followed by optional +/- and digits
      if (i < len && (expr[i] === 'e' || expr[i] === 'E')) {
        const next = i + 1 < len ? expr[i + 1] : '';
        if (isDigit(next) || ((next === '-' || next === '+') && i + 2 < len && isDigit(expr[i + 2]))) {
          i++; // skip e/E
          if (expr[i] === '-' || expr[i] === '+')
            i++; // skip sign
          while (i < len && isDigit(expr[i]))
            i++;
        }
      }
      tokens.push({type: 'NUMBER', value: expr.slice(start, i), position: start});
      continue;
    }

    // Identifiers: letters, digits, underscores; also Math.ceil / Math.floor
    if (isAlpha(ch) || ch === '_') {
      const start = i;
      while (i < len && (isAlphaNum(expr[i]) || expr[i] === '_'))
        i++;
      // Handle Math.xxx
      if (expr.slice(start, i) === 'Math' && i < len && expr[i] === '.') {
        i++; // skip dot
        while (i < len && (isAlphaNum(expr[i]) || expr[i] === '_'))
          i++;
      }
      tokens.push({type: 'IDENTIFIER', value: expr.slice(start, i), position: start});
      continue;
    }

    // Two-character tokens
    if (i + 1 < len) {
      const two = expr[i] + expr[i + 1];
      if (two === '=>') {
        tokens.push({type: 'ARROW', value: '=>', position: i});
        i += 2;
        continue;
      }
      if (two === '**') {
        tokens.push({type: 'OPERATOR', value: '**', position: i});
        i += 2;
        continue;
      }
      if (two === '>=' || two === '<=' || two === '==' || two === '!=') {
        tokens.push({type: 'COMPARISON', value: two, position: i});
        i += 2;
        continue;
      }
      if (two === '+=') {
        tokens.push({type: 'ASSIGN', value: '+=', position: i});
        i += 2;
        continue;
      }
    }

    // Single-character tokens
    if (ch === '(') {
      tokens.push({type: 'LPAREN', value: '(', position: i});
      i++;
      continue;
    }
    if (ch === ')') {
      tokens.push({type: 'RPAREN', value: ')', position: i});
      i++;
      continue;
    }
    if (ch === ',') {
      tokens.push({type: 'COMMA', value: ',', position: i});
      i++;
      continue;
    }
    if (ch === '?') {
      tokens.push({type: 'QUESTION', value: '?', position: i});
      i++;
      continue;
    }
    if (ch === ':') {
      tokens.push({type: 'COLON', value: ':', position: i});
      i++;
      continue;
    }
    if (ch === ';') {
      tokens.push({type: 'SEMICOLON', value: ';', position: i});
      i++;
      continue;
    }
    if (ch === '=') {
      tokens.push({type: 'ASSIGN', value: '=', position: i});
      i++;
      continue;
    }

    // Comparison: single < or >
    if (ch === '>' || ch === '<') {
      tokens.push({type: 'COMPARISON', value: ch, position: i});
      i++;
      continue;
    }

    // Operators: + - * /
    if (ch === '+' || ch === '-' || ch === '*' || ch === '/') {
      tokens.push({type: 'OPERATOR', value: ch, position: i});
      i++;
      continue;
    }

    // Unknown character — skip
    i++;
  }

  return tokens;
}

function isDigit(ch: string): boolean {
  return ch >= '0' && ch <= '9';
}

function isAlpha(ch: string): boolean {
  return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z');
}

function isAlphaNum(ch: string): boolean {
  return isDigit(ch) || isAlpha(ch);
}
