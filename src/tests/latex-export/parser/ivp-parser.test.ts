import {parseIvp} from '../../../latex-export/parser/ivp-parser';

describe('parseIvp', () => {
  describe('basic structure', () => {
    it('should parse #name', () => {
      const result = parseIvp('#name: Robertson');
      expect(result.name).toBe('Robertson');
    });

    it('should parse #description', () => {
      const result = parseIvp('#name: Test\n#description: A test model');
      expect(result.description).toBe('A test model');
    });

    it('should parse multi-line #comment', () => {
      const input = `#name: Test
#comment:
  Source: https://example.com
  This is a multi-line comment.`;
      const result = parseIvp(input);
      expect(result.comment).toContain('Source: https://example.com');
      expect(result.comment).toContain('multi-line comment');
    });
  });

  describe('equations', () => {
    it('should parse simple derivative', () => {
      const input = `#equations:
  dy/dt = -y + sin(t) / t`;
      const result = parseIvp(input);
      expect(result.equations).toHaveLength(1);
      expect(result.equations[0].lhs).toBe('dy/dt');
      expect(result.equations[0].rhs).toBe('-y + sin(t) / t');
      expect(result.equations[0].isDerivative).toBe(true);
    });

    it('should parse derivative with parenthesized variable', () => {
      const input = `#equations:
  d(FFox)/dt = -E11 + E12`;
      const result = parseIvp(input);
      expect(result.equations[0].lhs).toBe('d(FFox)/dt');
      expect(result.equations[0].isDerivative).toBe(true);
    });

    it('should parse multiple equations', () => {
      const input = `#equations:
  dx/dt = -k1 * x + k2 * y**2
  dy/dt = 2 * k1 * x - 2 * k2 * y**2`;
      const result = parseIvp(input);
      expect(result.equations).toHaveLength(2);
    });

    it('should handle multi-line equation', () => {
      const input = `#equations:
  d(MEAthiol)/dt = 2 * (-E11 + E12 - E21 + E22)
                   - (MEAthiol + MA) * (Fin + Fper) / VL`;
      const result = parseIvp(input);
      expect(result.equations).toHaveLength(1);
      expect(result.equations[0].rhs).toContain('- (MEAthiol + MA)');
    });
  });

  describe('expressions', () => {
    it('should parse simple expression', () => {
      const input = `#expressions:
  E1 = C1 * exp(-t) + P1`;
      const result = parseIvp(input);
      expect(result.expressions).toHaveLength(1);
      expect(result.expressions[0].lhs).toBe('E1');
      expect(result.expressions[0].rhs).toBe('C1 * exp(-t) + P1');
    });

    it('should detect arrow function', () => {
      const input = `#expressions:
  func1 = (p, t) => (p < 8) ? ceil(t) : ceil(t**2 / 20)`;
      const result = parseIvp(input);
      expect(result.expressions[0].isArrowFunction).toBe(true);
      expect(result.expressions[0].arrowParams).toEqual(['p', 't']);
    });

    it('should strip comments from expressions', () => {
      const input = `#expressions:
  // this is a comment
  E1 = x + y
  E2 = x - y  // inline comment`;
      const result = parseIvp(input);
      expect(result.expressions).toHaveLength(2);
      expect(result.expressions[0].lhs).toBe('E1');
      expect(result.expressions[1].rhs).toBe('x - y');
    });

    it('should handle alias expressions', () => {
      const input = `#expressions:
  ceil = Math.ceil`;
      const result = parseIvp(input);
      expect(result.expressions).toHaveLength(1);
      expect(result.expressions[0].rhs).toBe('Math.ceil');
    });
  });

  describe('inits', () => {
    it('should parse initial values with annotations', () => {
      const input = `#inits:
  x = 2 {units: mol/L; category: Initial values; min: 0; max: 5} [Initial value of x]
  y = 0 {units: mol/L; category: Initial values; min: -2; max: 2} [Initial value of y]`;
      const result = parseIvp(input);
      expect(result.inits).toHaveLength(2);
      expect(result.inits[0].name).toBe('x');
      expect(result.inits[0].value).toBe('2');
      expect(result.inits[0].units).toBe('mol/L');
      expect(result.inits[0].category).toBe('Initial values');
      expect(result.inits[0].tooltip).toBe('Initial value of x');
    });

    it('should parse init without annotations', () => {
      const input = `#inits:
  y = 0`;
      const result = parseIvp(input);
      expect(result.inits[0].name).toBe('y');
      expect(result.inits[0].value).toBe('0');
      expect(result.inits[0].units).toBeUndefined();
    });

    it('should parse scientific notation values', () => {
      const input = `#inits:
  X = 9.2E-2 {units: mg/ml}`;
      const result = parseIvp(input);
      expect(result.inits[0].value).toBe('9.2E-2');
    });
  });

  describe('parameters', () => {
    it('should parse parameters with full annotations', () => {
      const input = `#parameters:
  k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}
  k2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}`;
      const result = parseIvp(input);
      expect(result.parameters).toHaveLength(2);
      expect(result.parameters[0].name).toBe('k1');
      expect(result.parameters[0].value).toBe('0.7');
    });

    it('should parse parameter with caption and tooltip', () => {
      const input = `#parameters:
  P1 = 1.1 {caption: frequency; category: Parameters; min: 1; max: 10} [The control parameter]`;
      const result = parseIvp(input);
      expect(result.parameters[0].caption).toBe('frequency');
      expect(result.parameters[0].tooltip).toBe('The control parameter');
    });
  });

  describe('constants', () => {
    it('should parse constants', () => {
      const input = `#constants:
  C1 = 1
  C2 = 3`;
      const result = parseIvp(input);
      expect(result.constants).toHaveLength(2);
      expect(result.constants[0].name).toBe('C1');
      expect(result.constants[0].value).toBe('1');
    });

    it('should parse constants with various formats', () => {
      const input = `#constants:
   VLinit = 7.2
      VTV = 10
    speed = 400`;
      const result = parseIvp(input);
      expect(result.constants).toHaveLength(3);
      expect(result.constants[0].name).toBe('VLinit');
      expect(result.constants[0].value).toBe('7.2');
    });
  });

  describe('argument', () => {
    it('should parse argument block', () => {
      const input = `#argument: t
  start = 0 {caption: Initial; category: Time}
  finish = 10 {caption: Final; category: Time}
  step = 0.01 {caption: Step; category: Time}`;
      const result = parseIvp(input);
      expect(result.argument.name).toBe('t');
      expect(result.argument.entries).toHaveLength(3);
    });

    it('should parse argument with stage name', () => {
      const input = `#argument: t, 1-st stage
  start = 0
  stage1 = 60
  step = 0.1`;
      const result = parseIvp(input);
      expect(result.argument.name).toBe('t');
    });
  });

  describe('loop', () => {
    it('should parse loop block', () => {
      const input = `#loop:
  count = 10 {caption: count; category: Dosing}
  depot += dose`;
      const result = parseIvp(input);
      expect(result.loops).toHaveLength(1);
    });
  });

  describe('update', () => {
    it('should parse update block with stage name', () => {
      const input = `#update: 2-nd stage
  duration = overall - _t1
  S += 70`;
      const result = parseIvp(input);
      expect(result.updates).toHaveLength(1);
    });
  });

  describe('tolerance', () => {
    it('should parse tolerance', () => {
      const input = '#tolerance: 1e-7';
      const result = parseIvp(input);
      expect(result.tolerance).toBe('1e-7');
    });
  });

  describe('full model', () => {
    it('should parse basic.ivp correctly', () => {
      const input = `#name: Template
#equations:
  dy/dt = -y + sin(t) / t

#argument: t
  initial = 0.01 {min: 0.01; max: 10}
    final = 15   {min: 15; max: 150}
     step = 0.01 {min: 0.001; max: 0.1}

#inits:
  y = 0 {min: 0; max: 9}`;

      const result = parseIvp(input);
      expect(result.name).toBe('Template');
      expect(result.equations).toHaveLength(1);
      expect(result.argument.name).toBe('t');
      expect(result.inits).toHaveLength(1);
    });

    it('should parse robertson.ivp correctly', () => {
      const input = `#name: Robertson
#description: Robertson's chemical reaction model
#comment: This is classic example of stiff ODEs.
#equations:
  dA/dt = -0.04 * A + 1e4 * B * C
  dB/dt = 0.04 * A - 1e4 * B * C - 3e7 * B**2
  dC/dt = 3e7 * B**2

#inits:
  A = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  B = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  C = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}

#argument: t
  start = 0 {units: sec; caption: Initial; category: Time}
  finish = 40 {units: sec; caption: Final; category: Time}
  step = 0.01 {units: sec; caption: Step; category: Time}

#tolerance: 1e-7`;

      const result = parseIvp(input);
      expect(result.name).toBe('Robertson');
      expect(result.equations).toHaveLength(3);
      expect(result.equations[1].rhs).toContain('3e7 * B**2');
      expect(result.inits).toHaveLength(3);
      expect(result.tolerance).toBe('1e-7');
    });
  });
});
