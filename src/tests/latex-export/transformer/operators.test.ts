import {operatorToLatex} from '../../../latex-export/transformer/operators';

/**
 * These tests verify the LaTeX rendering of binary and unary operators.
 * `operatorToLatex` takes the operator, rendered left/right operands,
 * and context about whether brackets are needed.
 */

describe('operatorToLatex', () => {
  describe('addition and subtraction', () => {
    it('a + b → a + b', () => {
      expect(operatorToLatex('+', 'a', 'b')).toBe('a + b');
    });

    it('a - b → a - b', () => {
      expect(operatorToLatex('-', 'a', 'b')).toBe('a - b');
    });
  });

  describe('multiplication', () => {
    it('a * b → a \\cdot b (with cdot enabled)', () => {
      expect(operatorToLatex('*', 'a', 'b', {useCdot: true})).toBe('a \\cdot b');
    });

    it('2 * x → 2 \\cdot x (number × identifier)', () => {
      expect(operatorToLatex('*', '2', 'x', {useCdot: true})).toBe('2 \\cdot x');
    });

    it('a * b → a \\, b (with cdot disabled — juxtaposition)', () => {
      expect(operatorToLatex('*', 'a', 'b', {useCdot: false})).toBe('a \\, b');
    });
  });

  describe('division → fraction', () => {
    it('a / b → \\frac{a}{b}', () => {
      expect(operatorToLatex('/', 'a', 'b')).toBe('\\frac{a}{b}');
    });

    it('x + y / z in frac form', () => {
      // When / is top-level, it should produce a fraction
      expect(operatorToLatex('/', 'x + y', 'z')).toBe('\\frac{x + y}{z}');
    });
  });

  describe('exponentiation', () => {
    it('x ** 2 → x^{2}', () => {
      expect(operatorToLatex('**', 'x', '2')).toBe('x^{2}');
    });

    it('compound base needs parens: (x + y) ** 2', () => {
      expect(operatorToLatex('**', 'x + y', '2', {baseIsCompound: true}))
        .toBe('\\left(x + y\\right)^{2}');
    });

    it('simple base no parens: x ** n', () => {
      expect(operatorToLatex('**', 'x', 'n', {baseIsCompound: false}))
        .toBe('x^{n}');
    });
  });

  describe('comparison operators', () => {
    it('x < 0 → x < 0', () => {
      expect(operatorToLatex('<', 'x', '0')).toBe('x < 0');
    });

    it('x >= 0 → x \\geq 0', () => {
      expect(operatorToLatex('>=', 'x', '0')).toBe('x \\geq 0');
    });

    it('x <= 0 → x \\leq 0', () => {
      expect(operatorToLatex('<=', 'x', '0')).toBe('x \\leq 0');
    });

    it('a != b → a \\neq b', () => {
      expect(operatorToLatex('!=', 'a', 'b')).toBe('a \\neq b');
    });

    it('a == b → a = b (in math context)', () => {
      expect(operatorToLatex('==', 'a', 'b')).toBe('a = b');
    });
  });
});
