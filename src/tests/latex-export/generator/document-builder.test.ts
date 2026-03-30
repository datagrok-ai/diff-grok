import { describe, it, expect } from 'vitest';
import { buildLatexDocument, buildMarkdownDocument } from '../../../latex-export/generator/document-builder';
import { ParsedModel } from '../../../latex-export/types';

/** Minimal model fixture */
function minimalModel(): ParsedModel {
  return {
    name: 'Test Model',
    description: 'A test',
    comment: 'Some comment',
    equations: [
      { lhs: 'dy/dt', rhs: '-y + sin(t) / t', isDerivative: true, isArrowFunction: false },
    ],
    expressions: [],
    inits: [
      { name: 'y', value: '0', units: 'mol/L', tooltip: 'Initial y' },
    ],
    parameters: [],
    constants: [],
    argument: { name: 't', entries: [
      { name: 'start', value: '0' },
      { name: 'finish', value: '10' },
      { name: 'step', value: '0.01' },
    ]},
    loops: [],
    updates: [],
    output: [],
  };
}

/** Model with expressions and parameters */
function richModel(): ParsedModel {
  return {
    name: 'Extended Model',
    description: '2D system',
    equations: [
      { lhs: 'dx/dt', rhs: 'E1 * y + sin(t)', isDerivative: true, isArrowFunction: false },
      { lhs: 'dy/dt', rhs: 'E2 * x - pow(t, 5)', isDerivative: true, isArrowFunction: false },
    ],
    expressions: [
      { lhs: 'E1', rhs: 'C1 * exp(-t) + P1', isDerivative: false, isArrowFunction: false },
      { lhs: 'E2', rhs: 'C2 * cos(2 * t) + P2', isDerivative: false, isArrowFunction: false },
    ],
    inits: [
      { name: 'x', value: '2', units: 'C', category: 'Initial values', tooltip: 'Initial x' },
      { name: 'y', value: '0', units: 'C', category: 'Initial values', tooltip: 'Initial y' },
    ],
    parameters: [
      { name: 'P1', value: '1', category: 'Parameters', tooltip: 'P1 parameter' },
      { name: 'P2', value: '-1', category: 'Parameters', tooltip: 'P2 parameter' },
    ],
    constants: [
      { name: 'C1', value: '1' },
      { name: 'C2', value: '3' },
    ],
    argument: { name: 't', entries: [
      { name: 'start', value: '0', caption: 'Initial' },
      { name: 'finish', value: '10', caption: 'Final' },
      { name: 'step', value: '0.01', caption: 'Step' },
    ]},
    loops: [],
    updates: [],
    output: [],
    tolerance: '5e-5',
  };
}

describe('buildLatexDocument', () => {
  it('should include section with model name', () => {
    const result = buildLatexDocument(minimalModel());
    expect(result).toContain('\\section{Test Model}');
  });

  it('should include description', () => {
    const result = buildLatexDocument(minimalModel());
    expect(result).toContain('A test');
  });

  it('should include equations subsection with align environment', () => {
    const result = buildLatexDocument(minimalModel());
    expect(result).toContain('\\subsection{Equations}');
    expect(result).toContain('\\begin{align}');
    expect(result).toContain('\\end{align}');
  });

  it('should include derivative in equations', () => {
    const result = buildLatexDocument(minimalModel());
    expect(result).toContain('\\frac{dy}{dt}');
  });

  it('should include initial conditions table', () => {
    const result = buildLatexDocument(minimalModel());
    expect(result).toContain('Initial Conditions');
    expect(result).toContain('mol/L');
  });

  it('should include expressions subsection', () => {
    const result = buildLatexDocument(richModel());
    expect(result).toContain('\\subsection{Expressions}');
    expect(result).toContain('E_{1}');
    expect(result).toContain('E_{2}');
  });

  it('should include parameters table', () => {
    const result = buildLatexDocument(richModel());
    expect(result).toContain('Parameters');
    expect(result).toContain('P_{1}');
    expect(result).toContain('P_{2}');
  });

  it('should include constants table', () => {
    const result = buildLatexDocument(richModel());
    expect(result).toContain('Constants');
    expect(result).toContain('C_{1}');
    expect(result).toContain('C_{2}');
  });

  it('should skip empty sections', () => {
    const result = buildLatexDocument(minimalModel());
    // Minimal model has no expressions — no expressions section
    expect(result).not.toContain('\\subsection{Expressions}');
  });
});

describe('buildMarkdownDocument', () => {
  it('should use ## for model name', () => {
    const result = buildMarkdownDocument(minimalModel());
    expect(result).toContain('## Test Model');
  });

  it('should use ### for subsections', () => {
    const result = buildMarkdownDocument(minimalModel());
    expect(result).toContain('### Equations');
  });

  it('should wrap equations in $$ delimiters', () => {
    const result = buildMarkdownDocument(minimalModel());
    expect(result).toContain('$$');
  });

  it('should use aligned environment', () => {
    const result = buildMarkdownDocument(richModel());
    expect(result).toContain('\\begin{aligned}');
    expect(result).toContain('\\end{aligned}');
  });

  it('should use markdown table for parameters', () => {
    const result = buildMarkdownDocument(richModel());
    expect(result).toContain('|');
    expect(result).toContain('---');
  });

  it('should use inline $...$ for variables in tables', () => {
    const result = buildMarkdownDocument(richModel());
    // Table should have $P_{1}$ etc.
    expect(result).toMatch(/\$[^$]+\$/);
  });

  it('should render description as italic', () => {
    const result = buildMarkdownDocument(minimalModel());
    // Comment should be italic
    expect(result).toContain('*');
  });
});
