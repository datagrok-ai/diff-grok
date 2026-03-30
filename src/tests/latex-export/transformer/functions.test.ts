import { functionToLatex } from '../../../latex-export/transformer/functions';

/**
 * These tests verify the LaTeX rendering of function calls.
 * The `functionToLatex` function takes a function name and 
 * an array of already-rendered LaTeX argument strings.
 */

describe('functionToLatex', () => {
  describe('trigonometric functions', () => {
    it('sin(t) → \\sin\\!\\left(t\\right)', () => {
      const result = functionToLatex('sin', ['t']);
      expect(result).toBe('\\sin\\!\\left(t\\right)');
    });

    it('cos(t) → \\cos\\!\\left(t\\right)', () => {
      const result = functionToLatex('cos', ['t']);
      expect(result).toBe('\\cos\\!\\left(t\\right)');
    });

    it('tan(x) → \\tan\\!\\left(x\\right)', () => {
      const result = functionToLatex('tan', ['x']);
      expect(result).toBe('\\tan\\!\\left(x\\right)');
    });

    it('sin with complex argument', () => {
      const result = functionToLatex('sin', ['\\pi \\cdot P_{1}']);
      expect(result).toBe('\\sin\\!\\left(\\pi \\cdot P_{1}\\right)');
    });
  });

  describe('exp', () => {
    it('exp(-t) → e^{-t} (simple argument)', () => {
      const result = functionToLatex('exp', ['-t']);
      expect(result).toBe('e^{-t}');
    });

    it('exp(x) → e^{x}', () => {
      const result = functionToLatex('exp', ['x']);
      expect(result).toBe('e^{x}');
    });

    it('exp with complex argument should still use e^{...}', () => {
      const result = functionToLatex('exp', ['-t \\cdot k']);
      // Should use e^{...} notation
      expect(result).toMatch(/^e\^{/);
    });
  });

  describe('sqrt', () => {
    it('sqrt(x) → \\sqrt{x}', () => {
      const result = functionToLatex('sqrt', ['x']);
      expect(result).toBe('\\sqrt{x}');
    });

    it('sqrt with complex argument', () => {
      const result = functionToLatex('sqrt', ['E_{70}']);
      expect(result).toBe('\\sqrt{E_{70}}');
    });
  });

  describe('pow', () => {
    it('pow(x, 2) → x^{2}', () => {
      const result = functionToLatex('pow', ['x', '2']);
      expect(result).toBe('x^{2}');
    });

    it('pow with complex base: pow(a + b, 2) → \\left(a + b\\right)^{2}', () => {
      const result = functionToLatex('pow', ['a + b', '2']);
      expect(result).toBe('\\left(a + b\\right)^{2}');
    });

    it('pow with complex exponent: pow(x, n + 1) → x^{n + 1}', () => {
      const result = functionToLatex('pow', ['x', 'n + 1']);
      expect(result).toBe('x^{n + 1}');
    });

    it('pow(VL, -0.65) → \\mathrm{VL}^{-0.65} (base needs wrapping check)', () => {
      // Note: the base is already a LaTeX-rendered identifier
      const result = functionToLatex('pow', ['\\mathrm{VL}', '-0.65']);
      expect(result).toBe('\\mathrm{VL}^{-0.65}');
    });
  });

  describe('log', () => {
    it('log(x) → \\ln(x)', () => {
      const result = functionToLatex('log', ['x']);
      expect(result).toContain('\\ln');
    });
  });

  describe('abs', () => {
    it('abs(x) → \\left\\lvert x \\right\\rvert', () => {
      const result = functionToLatex('abs', ['x']);
      expect(result).toBe('\\left\\lvert x \\right\\rvert');
    });
  });

  describe('ceil / floor', () => {
    it('Math.ceil(x) → \\left\\lceil x \\right\\rceil', () => {
      const result = functionToLatex('Math.ceil', ['x']);
      expect(result).toBe('\\left\\lceil x \\right\\rceil');
    });

    it('ceil(t) → \\left\\lceil t \\right\\rceil (aliased)', () => {
      const result = functionToLatex('ceil', ['t']);
      expect(result).toBe('\\left\\lceil t \\right\\rceil');
    });

    it('Math.floor(x) → \\left\\lfloor x \\right\\rfloor', () => {
      const result = functionToLatex('Math.floor', ['x']);
      expect(result).toBe('\\left\\lfloor x \\right\\rfloor');
    });
  });

  describe('unknown functions', () => {
    it('should render unknown function as \\mathrm{name}(...)', () => {
      const result = functionToLatex('customFunc', ['x', 'y']);
      expect(result).toBe('\\mathrm{customFunc}\\!\\left(x, y\\right)');
    });
  });
});
