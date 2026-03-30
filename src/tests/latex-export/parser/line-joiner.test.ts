import { describe, it, expect } from 'vitest';
import { stripComments, joinMultiLineFormulas, stripAnnotations } from '../../../latex-export/parser/line-joiner';

describe('stripComments', () => {
  it('should remove trailing comment', () => {
    expect(stripComments('dx/dt = x + y  // 1st equation')).toBe('dx/dt = x + y');
  });

  it('should remove full-line comment', () => {
    expect(stripComments('  // this makes further code shorter')).toBe('');
  });

  it('should not touch lines without comments', () => {
    expect(stripComments('  dx/dt = x + y')).toBe('  dx/dt = x + y');
  });

  it('should preserve single slash (division)', () => {
    expect(stripComments('  a / b / c')).toBe('  a / b / c');
  });

  it('should handle comment immediately after expression (no space)', () => {
    expect(stripComments('x = 1//comment')).toBe('x = 1');
  });

  it('should handle line with only //', () => {
    expect(stripComments('//')).toBe('');
  });

  it('should handle expression with multiple divisions and a comment', () => {
    expect(stripComments('  a / b + c / d // note')).toBe('  a / b + c / d');
  });

  it('should trim trailing whitespace after stripping comment', () => {
    const result = stripComments('E1 = x   // decay');
    expect(result).toBe('E1 = x');
  });
});

describe('stripAnnotations', () => {
  it('should remove curly brace annotations', () => {
    expect(stripAnnotations('x = 2 {units: mol/L; min: 0; max: 5}')).toBe('x = 2');
  });

  it('should remove square bracket tooltips', () => {
    expect(stripAnnotations('P1 = 1 [The control parameter]')).toBe('P1 = 1');
  });

  it('should remove both annotations and tooltips', () => {
    expect(stripAnnotations('x = 2 {units: mol/L; category: Initial values; min: 0; max: 5} [Initial value of x]'))
      .toBe('x = 2');
  });

  it('should handle lines without annotations', () => {
    expect(stripAnnotations('C1 = 1')).toBe('C1 = 1');
  });

  it('should handle complex annotation with nested content', () => {
    expect(stripAnnotations('k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}')).toBe('k1 = 0.7');
  });
});

describe('joinMultiLineFormulas', () => {
  it('should join continuation line (no = sign)', () => {
    const lines = [
      'd(MEAthiol)/dt = 2 * (-E11 + E12)',
      '                 - (MEAthiol + MA) * Fin / VL',
    ];
    const result = joinMultiLineFormulas(lines);
    expect(result).toHaveLength(1);
    expect(result[0]).toContain('2 * (-E11 + E12)');
    expect(result[0]).toContain('- (MEAthiol + MA) * Fin / VL');
  });

  it('should not join lines that contain =', () => {
    const lines = [
      'E1 = C1 * exp(-t) + P1',
      'E2 = C2 * cos(2 * t) + P2',
    ];
    const result = joinMultiLineFormulas(lines);
    expect(result).toHaveLength(2);
  });

  it('should handle three continuation lines', () => {
    const lines = [
      'dx/dt = a + b',
      '        + c + d',
      '        + e + f',
    ];
    const result = joinMultiLineFormulas(lines);
    expect(result).toHaveLength(1);
    expect(result[0]).toContain('+ e + f');
  });

  it('should skip empty lines', () => {
    const lines = [
      'E1 = x',
      '',
      'E2 = y',
    ];
    const result = joinMultiLineFormulas(lines);
    expect(result).toHaveLength(2);
  });
});
