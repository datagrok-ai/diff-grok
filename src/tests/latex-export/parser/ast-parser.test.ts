import {parseExpression, ASTNode} from '../../../latex-export/parser/ast-parser';

/** Helper to build expected nodes concisely */
const num = (v: string): ASTNode => ({type: 'number', value: v});
const id = (n: string): ASTNode => ({type: 'identifier', name: n});
const bin = (op: string, left: ASTNode, right: ASTNode): ASTNode => ({type: 'binary', op, left, right});
const un = (op: string, operand: ASTNode): ASTNode => ({type: 'unary', op, operand});
const call = (name: string, args: ASTNode[]): ASTNode => ({type: 'call', name, args});

describe('parseExpression', () => {
  describe('atoms', () => {
    it('should parse a number', () => {
      expect(parseExpression('42')).toEqual(num('42'));
    });

    it('should parse scientific notation', () => {
      expect(parseExpression('1e4')).toEqual(num('1e4'));
    });

    it('should parse an identifier', () => {
      expect(parseExpression('x')).toEqual(id('x'));
    });
  });

  describe('binary operations', () => {
    it('should parse addition', () => {
      const ast = parseExpression('a + b');
      expect(ast).toEqual(bin('+', id('a'), id('b')));
    });

    it('should parse subtraction', () => {
      const ast = parseExpression('a - b');
      expect(ast).toEqual(bin('-', id('a'), id('b')));
    });

    it('should parse multiplication', () => {
      const ast = parseExpression('a * b');
      expect(ast).toEqual(bin('*', id('a'), id('b')));
    });

    it('should parse division', () => {
      const ast = parseExpression('a / b');
      expect(ast).toEqual(bin('/', id('a'), id('b')));
    });

    it('should respect * over + precedence', () => {
      const ast = parseExpression('a + b * c');
      expect(ast).toEqual(bin('+', id('a'), bin('*', id('b'), id('c'))));
    });

    it('should respect left-to-right for same precedence', () => {
      // a - b + c → (a - b) + c
      const ast = parseExpression('a - b + c');
      expect(ast).toEqual(bin('+', bin('-', id('a'), id('b')), id('c')));
    });
  });

  describe('exponentiation', () => {
    it('should parse **', () => {
      const ast = parseExpression('x**2');
      expect(ast).toEqual(bin('**', id('x'), num('2')));
    });

    it('should be right-associative: a**b**c → a**(b**c)', () => {
      const ast = parseExpression('a**b**c');
      expect(ast).toEqual(bin('**', id('a'), bin('**', id('b'), id('c'))));
    });

    it('should have higher precedence than *', () => {
      // a * b**2 → a * (b**2)
      const ast = parseExpression('a * b**2');
      expect(ast).toEqual(bin('*', id('a'), bin('**', id('b'), num('2'))));
    });
  });

  describe('unary', () => {
    it('should parse unary minus', () => {
      const ast = parseExpression('-x');
      expect(ast).toEqual(un('-', id('x')));
    });

    it('should parse unary minus before expression', () => {
      // -a * b → (-a) * b
      const ast = parseExpression('-a * b');
      expect(ast).toEqual(bin('*', un('-', id('a')), id('b')));
    });

    it('should parse double unary minus', () => {
      const ast = parseExpression('--x');
      expect(ast).toEqual(un('-', un('-', id('x'))));
    });
  });

  describe('parentheses', () => {
    it('should handle grouped expression', () => {
      const ast = parseExpression('(a + b) * c');
      expect(ast).toEqual(bin('*', bin('+', id('a'), id('b')), id('c')));
    });

    it('should handle nested parentheses', () => {
      const ast = parseExpression('((a + b))');
      expect(ast).toEqual(bin('+', id('a'), id('b')));
    });

    it('should handle power of grouped expression', () => {
      // (x + y)**2
      const ast = parseExpression('(x + y)**2');
      expect(ast).toEqual(bin('**', bin('+', id('x'), id('y')), num('2')));
    });
  });

  describe('function calls', () => {
    it('should parse sin(t)', () => {
      const ast = parseExpression('sin(t)');
      expect(ast).toEqual(call('sin', [id('t')]));
    });

    it('should parse pow(x, 2)', () => {
      const ast = parseExpression('pow(x, 2)');
      expect(ast).toEqual(call('pow', [id('x'), num('2')]));
    });

    it('should parse nested functions: exp(-t)', () => {
      const ast = parseExpression('exp(-t)');
      expect(ast).toEqual(call('exp', [un('-', id('t'))]));
    });

    it('should parse sqrt with complex arg', () => {
      const ast = parseExpression('sqrt(E70)');
      expect(ast).toEqual(call('sqrt', [id('E70')]));
    });

    it('should parse Math.ceil(x)', () => {
      const ast = parseExpression('Math.ceil(x)');
      expect(ast).toEqual(call('Math.ceil', [id('x')]));
    });
  });

  describe('ternary', () => {
    it('should parse simple ternary', () => {
      const ast = parseExpression('x < 0 ? a : b');
      expect(ast.type).toBe('ternary');
      if (ast.type === 'ternary') {
        expect(ast.condition).toEqual(bin('<', id('x'), num('0')));
        expect(ast.consequent).toEqual(id('a'));
        expect(ast.alternate).toEqual(id('b'));
      }
    });

    it('should parse ternary with function calls', () => {
      const ast = parseExpression('(E70 >= 0) ? sqrt(E70) : 0');
      expect(ast.type).toBe('ternary');
    });

    it('should parse nested ternary', () => {
      const ast = parseExpression('a < 1 ? x : b < 2 ? y : z');
      expect(ast.type).toBe('ternary');
      if (ast.type === 'ternary')
        expect(ast.alternate.type).toBe('ternary');
    });
  });

  describe('complex expressions from real models', () => {
    it('should parse: -0.04 * A + 1e4 * B * C', () => {
      const ast = parseExpression('-0.04 * A + 1e4 * B * C');
      expect(ast.type).toBe('binary');
      expect((ast as any).op).toBe('+');
    });

    it('should parse: 0.04 * A - 1e4 * B * C - 3e7 * B**2', () => {
      const ast = parseExpression('0.04 * A - 1e4 * B * C - 3e7 * B**2');
      expect(ast.type).toBe('binary');
    });

    it('should parse: V * S / (K + S) * X', () => {
      const ast = parseExpression('V * S / (K + S) * X');
      expect(ast.type).toBe('binary');
    });

    it('should parse: k2Fa * (Ffree * E0)**2 * E1', () => {
      const ast = parseExpression('k2Fa * (Ffree * E0)**2 * E1');
      expect(ast.type).toBe('binary');
    });

    it('should parse: 2 * (-E11 + E12 - E21 + E22)', () => {
      const ast = parseExpression('2 * (-E11 + E12 - E21 + E22)');
      expect(ast.type).toBe('binary');
      expect((ast as any).op).toBe('*');
    });

    it('should parse bioreactor multi-line expression', () => {
      // This tests a joined multi-line formula
      const expr = '(-(CL * A3 / V1 + Q / V1) * A1 + Q / V2 * A2 - kint * Rtot * A1 / (Kss + A1 / V1)) / (1 + Rtot * Kss / (Kss + A1 / V1)**2)';
      const ast = parseExpression(expr);
      expect(ast.type).toBe('binary');
      expect((ast as any).op).toBe('/');
    });
  });
});
