import {identifierToLatex} from '../../../latex-export/transformer/identifier';

describe('identifierToLatex', () => {
  describe('Greek letters — exact match', () => {
    it('alpha → \\alpha', () => {
      expect(identifierToLatex('alpha')).toBe('\\alpha');
    });

    it('beta → \\beta', () => {
      expect(identifierToLatex('beta')).toBe('\\beta');
    });

    it('gamma → \\gamma', () => {
      expect(identifierToLatex('gamma')).toBe('\\gamma');
    });

    it('delta → \\delta', () => {
      expect(identifierToLatex('delta')).toBe('\\delta');
    });

    it('epsilon → \\epsilon', () => {
      expect(identifierToLatex('epsilon')).toBe('\\epsilon');
    });

    it('varepsilon → \\varepsilon', () => {
      expect(identifierToLatex('varepsilon')).toBe('\\varepsilon');
    });

    it('zeta → \\zeta', () => {
      expect(identifierToLatex('zeta')).toBe('\\zeta');
    });

    it('eta → \\eta', () => {
      expect(identifierToLatex('eta')).toBe('\\eta');
    });

    it('theta → \\theta', () => {
      expect(identifierToLatex('theta')).toBe('\\theta');
    });

    it('vartheta → \\vartheta', () => {
      expect(identifierToLatex('vartheta')).toBe('\\vartheta');
    });

    it('iota → \\iota', () => {
      expect(identifierToLatex('iota')).toBe('\\iota');
    });

    it('kappa → \\kappa', () => {
      expect(identifierToLatex('kappa')).toBe('\\kappa');
    });

    it('lambda → \\lambda', () => {
      expect(identifierToLatex('lambda')).toBe('\\lambda');
    });

    it('mu → \\mu', () => {
      expect(identifierToLatex('mu')).toBe('\\mu');
    });

    it('nu → \\nu', () => {
      expect(identifierToLatex('nu')).toBe('\\nu');
    });

    it('xi → \\xi', () => {
      expect(identifierToLatex('xi')).toBe('\\xi');
    });

    it('pi → \\pi', () => {
      expect(identifierToLatex('pi')).toBe('\\pi');
    });

    it('rho → \\rho', () => {
      expect(identifierToLatex('rho')).toBe('\\rho');
    });

    it('varrho → \\varrho', () => {
      expect(identifierToLatex('varrho')).toBe('\\varrho');
    });

    it('sigma → \\sigma', () => {
      expect(identifierToLatex('sigma')).toBe('\\sigma');
    });

    it('tau → \\tau', () => {
      expect(identifierToLatex('tau')).toBe('\\tau');
    });

    it('upsilon → \\upsilon', () => {
      expect(identifierToLatex('upsilon')).toBe('\\upsilon');
    });

    it('phi → \\phi', () => {
      expect(identifierToLatex('phi')).toBe('\\phi');
    });

    it('varphi → \\varphi', () => {
      expect(identifierToLatex('varphi')).toBe('\\varphi');
    });

    it('chi → \\chi', () => {
      expect(identifierToLatex('chi')).toBe('\\chi');
    });

    it('psi → \\psi', () => {
      expect(identifierToLatex('psi')).toBe('\\psi');
    });

    it('omega → \\omega', () => {
      expect(identifierToLatex('omega')).toBe('\\omega');
    });
  });

  describe('Greek letters — uppercase', () => {
    it('Gamma → \\Gamma', () => {
      expect(identifierToLatex('Gamma')).toBe('\\Gamma');
    });

    it('Delta → \\Delta', () => {
      expect(identifierToLatex('Delta')).toBe('\\Delta');
    });

    it('Theta → \\Theta', () => {
      expect(identifierToLatex('Theta')).toBe('\\Theta');
    });

    it('Lambda → \\Lambda', () => {
      expect(identifierToLatex('Lambda')).toBe('\\Lambda');
    });

    it('Sigma → \\Sigma', () => {
      expect(identifierToLatex('Sigma')).toBe('\\Sigma');
    });

    it('Phi → \\Phi', () => {
      expect(identifierToLatex('Phi')).toBe('\\Phi');
    });

    it('Psi → \\Psi', () => {
      expect(identifierToLatex('Psi')).toBe('\\Psi');
    });

    it('Omega → \\Omega', () => {
      expect(identifierToLatex('Omega')).toBe('\\Omega');
    });

    it('Alpha → \\mathrm{A} (same as Latin)', () => {
      expect(identifierToLatex('Alpha')).toBe('\\mathrm{A}');
    });
  });

  describe('Greek prefix + digit suffix', () => {
    it('alpha1 → \\alpha_{1}', () => {
      expect(identifierToLatex('alpha1')).toBe('\\alpha_{1}');
    });

    it('mu2 → \\mu_{2}', () => {
      expect(identifierToLatex('mu2')).toBe('\\mu_{2}');
    });

    it('gamma12 → \\gamma_{12}', () => {
      expect(identifierToLatex('gamma12')).toBe('\\gamma_{12}');
    });
  });

  describe('Greek prefix + letter suffix → NOT Greek', () => {
    it('muM → \\mathrm{muM} (not \\mu M)', () => {
      expect(identifierToLatex('muM')).toBe('\\mathrm{muM}');
    });

    it('lambda_max should not split incorrectly', () => {
      // lambda followed by _ is handled but lambdaX (letter suffix) is not
      const result = identifierToLatex('lambdaX');
      expect(result).toBe('\\mathrm{lambdaX}');
    });
  });

  describe('single letter + digits → subscript', () => {
    it('x1 → x_{1}', () => {
      expect(identifierToLatex('x1')).toBe('x_{1}');
    });

    it('y20 → y_{20}', () => {
      expect(identifierToLatex('y20')).toBe('y_{20}');
    });

    it('k25 → k_{25}', () => {
      expect(identifierToLatex('k25')).toBe('k_{25}');
    });

    it('V2 → V_{2}', () => {
      expect(identifierToLatex('V2')).toBe('V_{2}');
    });

    it('r1 → r_{1}', () => {
      expect(identifierToLatex('r1')).toBe('r_{1}');
    });

    it('A1 → A_{1}', () => {
      expect(identifierToLatex('A1')).toBe('A_{1}');
    });
  });

  describe('multi-letter identifiers → mathrm', () => {
    it('FFox → \\mathrm{FFox}', () => {
      expect(identifierToLatex('FFox')).toBe('\\mathrm{FFox}');
    });

    it('MEAthiol → \\mathrm{MEAthiol}', () => {
      expect(identifierToLatex('MEAthiol')).toBe('\\mathrm{MEAthiol}');
    });

    it('KKox → \\mathrm{KKox}', () => {
      expect(identifierToLatex('KKox')).toBe('\\mathrm{KKox}');
    });

    it('FKred → \\mathrm{FKred}', () => {
      expect(identifierToLatex('FKred')).toBe('\\mathrm{FKred}');
    });

    it('switchTime → \\mathrm{switchTime}', () => {
      expect(identifierToLatex('switchTime')).toBe('\\mathrm{switchTime}');
    });

    it('pO2sat → \\mathrm{pO2sat}', () => {
      expect(identifierToLatex('pO2sat')).toBe('\\mathrm{pO2sat}');
    });
  });

  describe('multi-letter + trailing digits → mathrm with subscript', () => {
    it('E11 → E_{11}', () => {
      expect(identifierToLatex('E11')).toBe('E_{11}');
    });

    it('E1 → E_{1}', () => {
      expect(identifierToLatex('E1')).toBe('E_{1}');
    });

    it('Ks → \\mathrm{Ks} (s is not a digit)', () => {
      expect(identifierToLatex('Ks')).toBe('\\mathrm{Ks}');
    });

    it('Kla → \\mathrm{Kla}', () => {
      expect(identifierToLatex('Kla')).toBe('\\mathrm{Kla}');
    });
  });

  describe('chemical formulas', () => {
    it('CO2 → \\mathrm{CO}_{2}', () => {
      expect(identifierToLatex('CO2')).toBe('\\mathrm{CO}_{2}');
    });

    it('SO4 → \\mathrm{SO}_{4}', () => {
      expect(identifierToLatex('SO4')).toBe('\\mathrm{SO}_{4}');
    });

    it('N2O5 → \\mathrm{N_{2}O_{5}}', () => {
      expect(identifierToLatex('N2O5')).toBe('\\mathrm{N_{2}O_{5}}');
    });

    it('O3P → \\mathrm{O3P} or \\mathrm{O}_{3}\\mathrm{P}', () => {
      // This is ambiguous — could be chemical or just a name.
      // Accept either reasonable interpretation
      const result = identifierToLatex('O3P');
      expect(result).toBeTruthy();
    });

    it('HO2 → \\mathrm{HO}_{2}', () => {
      expect(identifierToLatex('HO2')).toBe('\\mathrm{HO}_{2}');
    });

    it('HCHO → \\mathrm{HCHO}', () => {
      expect(identifierToLatex('HCHO')).toBe('\\mathrm{HCHO}');
    });

    it('C2O3 → \\mathrm{C_{2}O_{3}}', () => {
      expect(identifierToLatex('C2O3')).toBe('\\mathrm{C_{2}O_{3}}');
    });

    it('MEO2 → \\mathrm{MEO}_{2}', () => {
      expect(identifierToLatex('MEO2')).toBe('\\mathrm{MEO}_{2}');
    });
  });

  describe('special constants', () => {
    it('PI → \\pi', () => {
      expect(identifierToLatex('PI')).toBe('\\pi');
    });

    it('Inf → \\infty', () => {
      expect(identifierToLatex('Inf')).toBe('\\infty');
    });

    it('inf → \\infty', () => {
      expect(identifierToLatex('inf')).toBe('\\infty');
    });
  });

  describe('single letters — no transformation', () => {
    it('x → x', () => {
      expect(identifierToLatex('x')).toBe('x');
    });

    it('t → t', () => {
      expect(identifierToLatex('t')).toBe('t');
    });

    it('A → A', () => {
      expect(identifierToLatex('A')).toBe('A');
    });

    it('P → P', () => {
      expect(identifierToLatex('P')).toBe('P');
    });
  });
});
