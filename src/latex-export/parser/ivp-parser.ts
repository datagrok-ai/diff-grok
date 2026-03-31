/* IVP file parser: parses the full text of an .ivp file into a structured ParsedModel. */

import {ParsedModel, FormulaLine, AnnotatedLine, OutputLine} from '../types';
import {stripComments, stripAnnotations, parseAnnotations, parseTooltip, joinMultiLineFormulas} from './line-joiner';

interface RawBlock {
  keyword: string;
  headerValue: string;
  lines: string[];
}

/** Split IVP text into raw blocks. */
function splitBlocks(text: string): RawBlock[] {
  const blocks: RawBlock[] = [];
  let current: RawBlock | null = null;

  for (const raw of text.split('\n')) {
    const headerMatch = raw.match(/^#([\w.]+):\s*(.*)/);
    if (headerMatch) {
      current = {keyword: headerMatch[1], headerValue: headerMatch[2].trim(), lines: []};
      blocks.push(current);
    } else if (current)
      current.lines.push(raw);
  }
  return blocks;
}

/** Parse a formula line (equation or expression). */
function parseFormulaLine(joined: string): FormulaLine {
  const trimmed = joined.trim();

  // Derivative pattern: d(var)/dt or dvar/dt
  const derivMatch = trimmed.match(/^(d\([^)]+\)\/d\w+|d\w+\/d\w+)\s*=\s*(.*)/);
  if (derivMatch) {
    return {
      lhs: derivMatch[1].trim(),
      rhs: derivMatch[2].trim(),
      isDerivative: true,
      isArrowFunction: false,
    };
  }

  // Regular assignment: name = expr
  const eqIdx = trimmed.indexOf('=');
  if (eqIdx === -1)
    return {lhs: trimmed, rhs: '', isDerivative: false, isArrowFunction: false};

  const lhs = trimmed.substring(0, eqIdx).trim();
  const rhs = trimmed.substring(eqIdx + 1).trim();

  // Arrow function: (params) => expr
  const arrowMatch = rhs.match(/^\(([^)]*)\)\s*=>\s*(.*)/);
  if (arrowMatch) {
    const params = arrowMatch[1].split(',').map((p) => p.trim()).filter(Boolean);
    return {
      lhs,
      rhs: arrowMatch[2].trim(),
      isDerivative: false,
      isArrowFunction: true,
      arrowParams: params,
    };
  }

  return {lhs, rhs, isDerivative: false, isArrowFunction: false};
}

/** Parse an annotated line (inits, parameters, constants, argument entries). */
function parseAnnotatedLine(raw: string): AnnotatedLine {
  const trimmed = raw.trim();

  // Extract tooltip
  const tooltip = parseTooltip(trimmed);

  // Extract annotations
  const annMatch = trimmed.match(/\{([^}]*)\}/);
  const annotations = annMatch ? parseAnnotations(annMatch[0]) : {};

  // Strip annotations and tooltip to get name = value
  const clean = stripAnnotations(trimmed);
  const eqIdx = clean.indexOf('=');

  let name: string;
  let value: string;
  if (eqIdx !== -1) {
    name = clean.substring(0, eqIdx).trim();
    value = clean.substring(eqIdx + 1).trim();
  } else {
    name = clean.trim();
    value = '';
  }

  const result: AnnotatedLine = {name, value};
  if (annotations.units) result.units = annotations.units;
  if (annotations.caption) result.caption = annotations.caption;
  if (annotations.category) result.category = annotations.category;
  if (tooltip) result.tooltip = tooltip;
  if (annotations.min) result.min = annotations.min;
  if (annotations.max) result.max = annotations.max;
  if (annotations.step) result.step = annotations.step;

  return result;
}

/** Parse an output line. */
function parseOutputLine(raw: string): OutputLine {
  const trimmed = raw.trim();
  const annMatch = trimmed.match(/\{([^}]*)\}/);
  const annotations = annMatch ? parseAnnotations(annMatch[0]) : {};
  const name = trimmed.replace(/\s*\{[^}]*\}/, '').trim();
  const result: OutputLine = {name};
  if (annotations.caption) result.caption = annotations.caption;
  return result;
}

/** Preprocess block content lines: strip comments, filter blanks, join multi-line. */
function preprocessLines(lines: string[]): string[] {
  const stripped = lines.map(stripComments);
  return joinMultiLineFormulas(stripped);
}

/** Parse the full text of an IVP file into a structured model.
 *  @param ivpText  raw content of an .ivp file
 *  @returns        parsed model structure
 */
export function parseIvp(ivpText: string): ParsedModel {
  const blocks = splitBlocks(ivpText);

  const model: ParsedModel = {
    equations: [],
    expressions: [],
    inits: [],
    parameters: [],
    constants: [],
    argument: {name: 't', entries: []},
    loops: [],
    updates: [],
    output: [],
  };

  for (const block of blocks) {
    switch (block.keyword) {
    case 'name':
      model.name = block.headerValue;
      break;

    case 'description':
      model.description = block.headerValue;
      break;

    case 'comment': {
      const commentLines = block.lines
        .map((l) => l.trimStart())
        .filter((l) => l.length > 0);
      if (block.headerValue)
        commentLines.unshift(block.headerValue);
      model.comment = commentLines.join('\n');
      break;
    }

    case 'equations': {
      const joined = preprocessLines(block.lines);
      for (const line of joined) {
        const trimmed = line.trim();
        if (trimmed) model.equations.push(parseFormulaLine(trimmed));
      }
      break;
    }

    case 'expressions': {
      const joined = preprocessLines(block.lines);
      for (const line of joined) {
        const trimmed = line.trim();
        if (trimmed) model.expressions.push(parseFormulaLine(trimmed));
      }
      break;
    }

    case 'inits': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      for (const line of lines)
        model.inits.push(parseAnnotatedLine(line));
      break;
    }

    case 'parameters': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      for (const line of lines)
        model.parameters.push(parseAnnotatedLine(line));
      break;
    }

    case 'constants': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      for (const line of lines)
        model.constants.push(parseAnnotatedLine(line));
      break;
    }

    case 'argument': {
      const parts = block.headerValue.split(',');
      model.argument.name = parts[0].trim();
      if (parts.length > 1) model.argument.stageName = parts.slice(1).join(',').trim();
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      for (const line of lines)
        model.argument.entries.push(parseAnnotatedLine(line));
      break;
    }

    case 'loop': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      const entries: AnnotatedLine[] = [];
      for (const line of lines)
        entries.push(parseAnnotatedLine(line));
      model.loops.push({entries});
      break;
    }

    case 'update': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      const entries: AnnotatedLine[] = [];
      for (const line of lines)
        entries.push(parseAnnotatedLine(line));
      model.updates.push({stageName: block.headerValue, entries});
      break;
    }

    case 'output': {
      const lines = block.lines.map(stripComments).filter((l) => l.trim());
      for (const line of lines)
        model.output.push(parseOutputLine(line));
      break;
    }

    case 'tolerance':
      model.tolerance = block.headerValue;
      break;

    case 'meta.solver':
      model.solverMeta = block.headerValue;
      break;

    default:
      break;
    }
  }

  return model;
}
