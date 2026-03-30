import { tokenize, Token } from '../../../latex-export/parser/tokenizer';

/** Helper: extract just the type-value pairs for easier assertion */
function tv(tokens: Token[]): Array<[string, string]> {
  return tokens.map(t => [t.type, t.value]);
}

describe('tokenizer', () => {
  describe('numbers', () => {
    it('should tokenize integer', () => {
      const tokens = tokenize('42');
      expect(tv(tokens)).toEqual([['NUMBER', '42']]);
    });

    it('should tokenize decimal', () => {
      const tokens = tokenize('0.04');
      expect(tv(tokens)).toEqual([['NUMBER', '0.04']]);
    });

    it('should tokenize scientific notation: 1e4', () => {
      const tokens = tokenize('1e4');
      expect(tv(tokens)).toEqual([['NUMBER', '1e4']]);
    });

    it('should tokenize scientific notation: 1E-2', () => {
      const tokens = tokenize('1E-2');
      expect(tv(tokens)).toEqual([['NUMBER', '1E-2']]);
    });

    it('should tokenize scientific notation: 9.2E-2', () => {
      const tokens = tokenize('9.2E-2');
      expect(tv(tokens)).toEqual([['NUMBER', '9.2E-2']]);
    });

    it('should tokenize scientific notation: 1.23e4', () => {
      const tokens = tokenize('1.23e4');
      expect(tv(tokens)).toEqual([['NUMBER', '1.23e4']]);
    });

    it('should tokenize scientific notation: 3e7', () => {
      const tokens = tokenize('3e7');
      expect(tv(tokens)).toEqual([['NUMBER', '3e7']]);
    });

    it('should tokenize negative exponent after operator: x + 1E-2', () => {
      const tokens = tokenize('x + 1E-2');
      expect(tokens).toHaveLength(3);
      expect(tokens[2].type).toBe('NUMBER');
      expect(tokens[2].value).toBe('1E-2');
    });
  });

  describe('identifiers', () => {
    it('should tokenize simple identifier', () => {
      expect(tv(tokenize('x'))).toEqual([['IDENTIFIER', 'x']]);
    });

    it('should tokenize identifier with digits', () => {
      expect(tv(tokenize('x1'))).toEqual([['IDENTIFIER', 'x1']]);
    });

    it('should tokenize multi-letter identifier', () => {
      expect(tv(tokenize('FFox'))).toEqual([['IDENTIFIER', 'FFox']]);
    });

    it('should tokenize identifier with dots: Math.ceil', () => {
      expect(tv(tokenize('Math.ceil'))).toEqual([['IDENTIFIER', 'Math.ceil']]);
    });

    it('should tokenize PI', () => {
      expect(tv(tokenize('PI'))).toEqual([['IDENTIFIER', 'PI']]);
    });

    it('should tokenize complex identifiers like pKa2MEA', () => {
      expect(tv(tokenize('pKa2MEA'))).toEqual([['IDENTIFIER', 'pKa2MEA']]);
    });
  });

  describe('operators', () => {
    it('should tokenize basic arithmetic', () => {
      const tokens = tokenize('a + b');
      expect(tv(tokens)).toEqual([
        ['IDENTIFIER', 'a'],
        ['OPERATOR', '+'],
        ['IDENTIFIER', 'b'],
      ]);
    });

    it('should tokenize ** as single operator', () => {
      const tokens = tokenize('x**2');
      expect(tv(tokens)).toEqual([
        ['IDENTIFIER', 'x'],
        ['OPERATOR', '**'],
        ['NUMBER', '2'],
      ]);
    });

    it('should tokenize * and ** correctly in sequence', () => {
      const tokens = tokenize('a * b**2');
      expect(tokens[1].value).toBe('*');
      expect(tokens[3].value).toBe('**');
    });

    it('should tokenize unary minus', () => {
      const tokens = tokenize('-x');
      expect(tv(tokens)).toEqual([
        ['OPERATOR', '-'],
        ['IDENTIFIER', 'x'],
      ]);
    });
  });

  describe('comparisons', () => {
    it('should tokenize >=', () => {
      const tokens = tokenize('x >= 0');
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'COMPARISON', value: '>=' }));
    });

    it('should tokenize <', () => {
      const tokens = tokenize('p < 8');
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'COMPARISON', value: '<' }));
    });

    it('should tokenize ==', () => {
      const tokens = tokenize('a == b');
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'COMPARISON', value: '==' }));
    });

    it('should tokenize !=', () => {
      const tokens = tokenize('a != b');
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'COMPARISON', value: '!=' }));
    });
  });

  describe('ternary and arrow', () => {
    it('should tokenize ternary operator', () => {
      const tokens = tokenize('(x < 0) ? a : b');
      const types = tokens.map(t => t.type);
      expect(types).toContain('QUESTION');
      expect(types).toContain('COLON');
    });

    it('should tokenize arrow function', () => {
      const tokens = tokenize('(p, t) => p + t');
      const types = tokens.map(t => t.type);
      expect(types).toContain('ARROW');
    });
  });

  describe('parentheses and commas', () => {
    it('should tokenize function call', () => {
      const tokens = tokenize('sin(t)');
      expect(tv(tokens)).toEqual([
        ['IDENTIFIER', 'sin'],
        ['LPAREN', '('],
        ['IDENTIFIER', 't'],
        ['RPAREN', ')'],
      ]);
    });

    it('should tokenize pow with two args', () => {
      const tokens = tokenize('pow(x, 2)');
      expect(tv(tokens)).toEqual([
        ['IDENTIFIER', 'pow'],
        ['LPAREN', '('],
        ['IDENTIFIER', 'x'],
        ['COMMA', ','],
        ['NUMBER', '2'],
        ['RPAREN', ')'],
      ]);
    });
  });

  describe('complex expressions from real models', () => {
    it('should tokenize Robertson equation RHS', () => {
      const tokens = tokenize('-0.04 * A + 1e4 * B * C');
      expect(tokens[0]).toEqual(expect.objectContaining({ type: 'OPERATOR', value: '-' }));
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'NUMBER', value: '0.04' }));
      expect(tokens[5]).toEqual(expect.objectContaining({ type: 'NUMBER', value: '1e4' }));
    });

    it('should tokenize expression with pow and **', () => {
      const tokens = tokenize('k2Fa * (Ffree * E0)**2 * E1');
      const ops = tokens.filter(t => t.type === 'OPERATOR').map(t => t.value);
      expect(ops).toContain('**');
      expect(ops).toContain('*');
    });

    it('should tokenize ternary with function calls', () => {
      const expr = '(E70 >= 0) ? sqrt(E70) : 0';
      const tokens = tokenize(expr);
      expect(tokens.some(t => t.type === 'COMPARISON' && t.value === '>=')).toBe(true);
      expect(tokens.some(t => t.type === 'QUESTION')).toBe(true);
      expect(tokens.some(t => t.type === 'COLON')).toBe(true);
    });

    it('should tokenize Fin/Fper conditional', () => {
      const expr = 't < switchTime ? 0 : 0.025';
      const tokens = tokenize(expr);
      expect(tokens[0]).toEqual(expect.objectContaining({ type: 'IDENTIFIER', value: 't' }));
      expect(tokens[1]).toEqual(expect.objectContaining({ type: 'COMPARISON', value: '<' }));
    });
  });
});
