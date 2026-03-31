import {readFileSync} from 'fs';
import {join} from 'path';
import {convertIvpToLatex} from '../../../latex-export/index';

/** Load an example IVP file from the latex-export examples/ directory */
function loadExample(name: string): string {
  return readFileSync(join(__dirname, '../../../latex-export/examples', name), 'utf-8');
}

describe('integration: convertIvpToLatex', () => {
  describe('basic.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain the model name', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('Template');
    });

    it('should contain a derivative', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\frac{dy}{dt}');
    });

    it('should contain sin function', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\sin');
    });

    it('should contain a fraction (sin(t)/t)', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\frac');
    });
  });

  describe('robertson.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain 3 derivatives', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      const derivCount = (result.match(/\\frac{d[A-Z]}{dt}/g) || []).length;
      expect(derivCount).toBe(3);
    });

    it('should render scientific notation', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('10^{');
    });

    it('should render B^{2}', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('B^{2}');
    });
  });

  describe('chem-react.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('chem-react.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain subscripted variables x_{1}, x_{2}, x_{3}, x_{4}', () => {
      const input = loadExample('chem-react.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('x_{1}');
      expect(result).toContain('x_{2}');
      expect(result).toContain('x_{3}');
      expect(result).toContain('x_{4}');
    });

    it('should contain subscripted parameters k_{1} through k_{6}', () => {
      const input = loadExample('chem-react.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      for (let i = 1; i <= 6; i++)
        expect(result).toContain(`k_{${i}}`);
    });

    it('should contain squared terms', () => {
      const input = loadExample('chem-react.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('^{2}');
    });
  });

  describe('extended.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('extended.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain exp and cos functions', () => {
      const input = loadExample('extended.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('e^{');
      expect(result).toContain('\\cos');
    });

    it('should contain pow(t, 5) rendered as t^{5}', () => {
      const input = loadExample('extended.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('t^{5}');
    });
  });

  describe('energy-n-control.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('energy-n-control.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should handle ternary → cases environment', () => {
      const input = loadExample('energy-n-control.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\begin{cases}');
      expect(result).toContain('\\end{cases}');
    });

    it('should handle PI → \\pi', () => {
      const input = loadExample('energy-n-control.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\pi');
    });

    it('should handle x**2 + y**2 for energy expression', () => {
      const input = loadExample('energy-n-control.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('x^{2}');
      expect(result).toContain('y^{2}');
    });

    it('should strip // comments from output', () => {
      const input = loadExample('energy-n-control.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).not.toContain('// this makes further code shorter');
      expect(result).not.toContain('// simple if-then-else');
    });
  });

  describe('fermentation.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('fermentation.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain V * S / (K + S) as fraction', () => {
      const input = loadExample('fermentation.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\frac');
    });
  });

  describe('pk.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('pk.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should handle d(depot)/dt and d(centr)/dt', () => {
      const input = loadExample('pk.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\mathrm{depot}');
      expect(result).toContain('\\mathrm{centr}');
    });
  });

  describe('pk-pd.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('pk-pd.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain 4 equations', () => {
      const input = loadExample('pk-pd.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      const derivCount = (result.match(/\\frac{d/g) || []).length;
      expect(derivCount).toBeGreaterThanOrEqual(4);
    });
  });

  describe('pollution.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('pollution.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain y_{1} through y_{20}', () => {
      const input = loadExample('pollution.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('y_{1}');
      expect(result).toContain('y_{20}');
    });

    it('should contain r_{1} through r_{25} in expressions', () => {
      const input = loadExample('pollution.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('r_{1}');
      expect(result).toContain('r_{25}');
    });

    it('should contain k_{1} through k_{25} in expressions', () => {
      const input = loadExample('pollution.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('k_{1}');
      expect(result).toContain('k_{25}');
    });
  });

  describe('ga-production.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('ga-production.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain Monod-type expression with fraction', () => {
      const input = loadExample('ga-production.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\frac');
    });

    it('should handle Greek-looking identifiers: alpha, beta, gamma, etc.', () => {
      const input = loadExample('ga-production.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\alpha');
      expect(result).toContain('\\beta');
      expect(result).toContain('\\gamma');
      expect(result).toContain('\\delta');
      expect(result).toContain('\\lambda');
      expect(result).toContain('\\phi');
    });
  });

  describe('nimotuzumab.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('nimotuzumab.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should contain gamma as Greek letter', () => {
      const input = loadExample('nimotuzumab.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\gamma');
    });

    it('should handle complex fraction (nimotuzumab dA1/dt)', () => {
      const input = loadExample('nimotuzumab.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      const fracCount = (result.match(/\\frac/g) || []).length;
      expect(fracCount).toBeGreaterThanOrEqual(3);
    });
  });

  describe('bioreactor.ivp', () => {
    it('should convert without errors', () => {
      const input = loadExample('bioreactor.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toBeTruthy();
    });

    it('should handle multi-line equations', () => {
      const input = loadExample('bioreactor.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\mathrm{MEAthiol}');
    });

    it('should handle ternary in expressions (Fin, Fper)', () => {
      const input = loadExample('bioreactor.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('\\begin{cases}');
    });

    it('should handle pow() calls', () => {
      const input = loadExample('bioreactor.ivp');
      const result = convertIvpToLatex(input, {format: 'latex'});
      expect(result).toContain('^{');
    });
  });

  describe('markdown output', () => {
    it('should produce valid markdown with $$ delimiters', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'markdown'});
      expect(result).toContain('$$');
      expect(result).toContain('##');
    });

    it('should produce inline math for initial conditions table', () => {
      const input = loadExample('chem-react.ivp');
      const result = convertIvpToLatex(input, {format: 'markdown', includeInits: true});
      expect(result).toContain('$');
    });

    it('should contain aligned environment', () => {
      const input = loadExample('basic.ivp');
      const result = convertIvpToLatex(input, {format: 'markdown'});
      expect(result).toContain('\\begin{aligned}');
    });
  });

  describe('options', () => {
    it('should exclude metadata when includeMetadata: false', () => {
      const input = loadExample('robertson.ivp');
      const result = convertIvpToLatex(input, {format: 'latex', includeMetadata: false});
      expect(result).not.toContain('Robertson');
    });

    it('should exclude parameters when includeParameters: false', () => {
      const input = loadExample('chem-react.ivp');
      const full = convertIvpToLatex(input, {format: 'latex', includeParameters: true});
      const without = convertIvpToLatex(input, {format: 'latex', includeParameters: false});
      expect(full.length).toBeGreaterThan(without.length);
    });

    it('should exclude constants when includeConstants: false', () => {
      const input = loadExample('bioreactor.ivp');
      const full = convertIvpToLatex(input, {format: 'latex', includeConstants: true});
      const without = convertIvpToLatex(input, {format: 'latex', includeConstants: false});
      expect(full.length).toBeGreaterThan(without.length);
    });
  });
});
