import { describe, it, expect } from 'vitest';
import { expressionToLatex, derivativeToLatex } from '../../../latex-export/generator/latex-generator';

describe('expressionToLatex', () => {
  describe('simple expressions', () => {
    it('-y + sin(t) / t', () => {
      const result = expressionToLatex('-y + sin(t) / t');
      // Should contain fraction and sin
      expect(result).toContain('\\frac');
      expect(result).toContain('\\sin');
    });

    it('3e7 * B**2', () => {
      const result = expressionToLatex('3e7 * B**2');
      expect(result).toContain('3 \\times 10^{7}');
      expect(result).toContain('B^{2}');
    });

    it('x**2 + y**2', () => {
      const result = expressionToLatex('x**2 + y**2');
      expect(result).toContain('x^{2}');
      expect(result).toContain('y^{2}');
    });
  });

  describe('scientific notation rendering', () => {
    it('1e4 → 1 \\times 10^{4} or 10^{4}', () => {
      const result = expressionToLatex('1e4');
      expect(result).toMatch(/10\^{4}/);
    });

    it('9.2E-2 → 9.2 \\times 10^{-2}', () => {
      const result = expressionToLatex('9.2E-2');
      expect(result).toContain('9.2 \\times 10^{-2}');
    });

    it('1.23e4 → 1.23 \\times 10^{4}', () => {
      const result = expressionToLatex('1.23e4');
      expect(result).toContain('1.23 \\times 10^{4}');
    });

    it('should render small numbers without scientific notation: 0.04', () => {
      const result = expressionToLatex('0.04');
      expect(result).toBe('0.04');
    });
  });

  describe('Robertson model equations', () => {
    it('-0.04 * A + 1e4 * B * C', () => {
      const result = expressionToLatex('-0.04 * A + 1e4 * B * C');
      expect(result).toContain('A');
      expect(result).toContain('B');
      expect(result).toContain('C');
    });

    it('0.04 * A - 1e4 * B * C - 3e7 * B**2', () => {
      const result = expressionToLatex('0.04 * A - 1e4 * B * C - 3e7 * B**2');
      expect(result).toContain('B^{2}');
    });
  });

  describe('chem-react model expressions', () => {
    it('-k1 * x1 + k2 * (x2)**2 + k3 * x1 * x3 - k4 * (x1)**2', () => {
      const result = expressionToLatex('-k1 * x1 + k2 * (x2)**2 + k3 * x1 * x3 - k4 * (x1)**2');
      expect(result).toContain('k_{1}');
      expect(result).toContain('x_{1}');
      expect(result).toContain('x_{2}');
      expect(result).toContain('^{2}');
    });
  });

  describe('expressions with functions', () => {
    it('C1 * exp(-t) + P1', () => {
      const result = expressionToLatex('C1 * exp(-t) + P1');
      expect(result).toContain('e^{-t}');
      expect(result).toContain('C_{1}');
      expect(result).toContain('P_{1}');
    });

    it('C2 * cos(2 * t) + P2', () => {
      const result = expressionToLatex('C2 * cos(2 * t) + P2');
      expect(result).toContain('\\cos');
      expect(result).toContain('C_{2}');
    });

    it('sin(PI * P1)', () => {
      const result = expressionToLatex('sin(PI * P1)');
      expect(result).toContain('\\sin');
      expect(result).toContain('\\pi');
    });

    it('pow(VL, -0.65) * 0.065', () => {
      const result = expressionToLatex('pow(VL, -0.65) * 0.065');
      expect(result).toContain('^{-0.65}');
    });

    it('sqrt(E70)', () => {
      const result = expressionToLatex('sqrt(E70)');
      expect(result).toContain('\\sqrt');
    });
  });

  describe('ternary → cases', () => {
    it('(E70 >= 0) ? sqrt(E70) : 0', () => {
      const result = expressionToLatex('(E70 >= 0) ? sqrt(E70) : 0');
      expect(result).toContain('\\begin{cases}');
      expect(result).toContain('\\end{cases}');
      expect(result).toContain('\\sqrt');
      expect(result).toContain('\\geq');
    });

    it('t < switchTime ? 0 : 0.025', () => {
      const result = expressionToLatex('t < switchTime ? 0 : 0.025');
      expect(result).toContain('\\begin{cases}');
      expect(result).toContain('\\mathrm{switchTime}');
    });
  });

  describe('complex bioreactor expressions', () => {
    it('k1red * FFox * E0 * E1', () => {
      const result = expressionToLatex('k1red * FFox * E0 * E1');
      expect(result).toContain('\\mathrm{k1red}');
      expect(result).toContain('\\mathrm{FFox}');
    });

    it('k2Fa * (Ffree * E0)**2 * E1', () => {
      const result = expressionToLatex('k2Fa * (Ffree * E0)**2 * E1');
      expect(result).toContain('^{2}');
    });

    it('(MA * E0)**2', () => {
      const result = expressionToLatex('(MA * E0)**2');
      expect(result).toContain('^{2}');
    });
  });

  describe('fermentation expressions', () => {
    it('V * S / (K + S) * X', () => {
      const result = expressionToLatex('V * S / (K + S) * X');
      expect(result).toContain('\\frac');
    });
  });
});

describe('derivativeToLatex', () => {
  it('dy/dt → \\frac{dy}{dt}', () => {
    const result = derivativeToLatex('dy/dt');
    expect(result).toBe('\\frac{dy}{dt}');
  });

  it('dx1/dt → \\frac{dx_{1}}{dt}', () => {
    const result = derivativeToLatex('dx1/dt');
    expect(result).toBe('\\frac{dx_{1}}{dt}');
  });

  it('d(FFox)/dt → \\frac{d\\,\\mathrm{FFox}}{dt}', () => {
    const result = derivativeToLatex('d(FFox)/dt');
    expect(result).toBe('\\frac{d\\,\\mathrm{FFox}}{dt}');
  });

  it('d(depot)/dt → \\frac{d\\,\\mathrm{depot}}{dt}', () => {
    const result = derivativeToLatex('d(depot)/dt');
    expect(result).toBe('\\frac{d\\,\\mathrm{depot}}{dt}');
  });

  it('dA/dt → \\frac{dA}{dt}', () => {
    const result = derivativeToLatex('dA/dt');
    expect(result).toBe('\\frac{dA}{dt}');
  });

  it('d(MEAthiol)/dt → \\frac{d\\,\\mathrm{MEAthiol}}{dt}', () => {
    const result = derivativeToLatex('d(MEAthiol)/dt');
    expect(result).toBe('\\frac{d\\,\\mathrm{MEAthiol}}{dt}');
  });

  it('dy1/dt → \\frac{dy_{1}}{dt}', () => {
    const result = derivativeToLatex('dy1/dt');
    expect(result).toBe('\\frac{dy_{1}}{dt}');
  });

  it('dy20/dt → \\frac{dy_{20}}{dt}', () => {
    const result = derivativeToLatex('dy20/dt');
    expect(result).toBe('\\frac{dy_{20}}{dt}');
  });

  it('dP/dt → \\frac{dP}{dt}', () => {
    const result = derivativeToLatex('dP/dt');
    expect(result).toBe('\\frac{dP}{dt}');
  });
});
