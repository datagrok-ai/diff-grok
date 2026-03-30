/* Document builder: assembles converted blocks into a complete LaTeX or Markdown document. */

import {ParsedModel, ConvertOptions, FormulaLine, AnnotatedLine} from '../types';
import {expressionToLatex, derivativeToLatex} from './latex-generator';
import {identifierToLatex} from '../transformer/identifier';

const DEFAULT_OPTS: ConvertOptions = {
  format: 'latex',
  includeMetadata: true,
  includeInits: true,
  includeParameters: true,
  includeConstants: true,
  useCdot: true,
};

/** Build a complete LaTeX document from a parsed model.
 *  @param model    parsed IVP model
 *  @param options  conversion options
 *  @returns        LaTeX document string
 */
export function buildLatexDocument(
  model: ParsedModel,
  options?: Partial<ConvertOptions>,
): string {
  const opts = {...DEFAULT_OPTS, ...options};
  const parts: string[] = [];

  if (opts.includeMetadata) {
    if (model.name)
      parts.push(`\\section{${model.name}}`);
    if (model.description)
      parts.push(model.description);
    if (model.comment)
      parts.push(model.comment);
  }

  if (model.equations.length > 0) {
    parts.push('\\subsection{Equations}');
    parts.push(buildLatexAlign(model.equations));
  }

  if (model.expressions.length > 0) {
    parts.push('\\subsection{Expressions}');
    parts.push(buildLatexAlign(model.expressions));
  }

  if (opts.includeInits && model.inits.length > 0) {
    parts.push('\\subsection{Initial Conditions}');
    parts.push(buildLatexTable(model.inits));
  }

  if (opts.includeParameters && model.parameters.length > 0) {
    parts.push('\\subsection{Parameters}');
    parts.push(buildLatexTable(model.parameters));
  }

  if (opts.includeConstants && model.constants.length > 0) {
    parts.push('\\subsection{Constants}');
    parts.push(buildLatexTable(model.constants));
  }

  return parts.join('\n\n');
}

/** Build a complete Markdown document from a parsed model.
 *  @param model    parsed IVP model
 *  @param options  conversion options
 *  @returns        Markdown document string
 */
export function buildMarkdownDocument(
  model: ParsedModel,
  options?: Partial<ConvertOptions>,
): string {
  const opts = {...DEFAULT_OPTS, ...options};
  const parts: string[] = [];

  if (opts.includeMetadata) {
    if (model.name)
      parts.push(`## ${model.name}`);
    if (model.description)
      parts.push(`*${model.description}*`);
    if (model.comment)
      parts.push(`*${model.comment}*`);
  }

  if (model.equations.length > 0) {
    parts.push('### Equations');
    parts.push(buildMarkdownAlign(model.equations));
  }

  if (model.expressions.length > 0) {
    parts.push('### Expressions');
    parts.push(buildMarkdownAlign(model.expressions));
  }

  if (opts.includeInits && model.inits.length > 0) {
    parts.push('### Initial Conditions');
    parts.push(buildMarkdownTable(model.inits));
  }

  if (opts.includeParameters && model.parameters.length > 0) {
    parts.push('### Parameters');
    parts.push(buildMarkdownTable(model.parameters));
  }

  if (opts.includeConstants && model.constants.length > 0) {
    parts.push('### Constants');
    parts.push(buildMarkdownTable(model.constants));
  }

  return parts.join('\n\n');
}

function renderFormulaLhs(f: FormulaLine): string {
  if (f.isDerivative)
    return derivativeToLatex(f.lhs);
  const lhs = identifierToLatex(f.lhs);
  if (f.isArrowFunction && f.arrowParams)
    return `${lhs}(${f.arrowParams.join(', ')})`;
  return lhs;
}

function buildLatexAlign(formulas: FormulaLine[]): string {
  const lines = formulas.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs);
    return `  ${lhs} &= ${rhs}`;
  });
  return `\\begin{align}\n${lines.join(' \\\\\n')}\n\\end{align}`;
}

function buildMarkdownAlign(formulas: FormulaLine[]): string {
  const lines = formulas.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs);
    return `  ${lhs} &= ${rhs}`;
  });
  return `$$\n\\begin{aligned}\n${lines.join(' \\\\\n')}\n\\end{aligned}\n$$`;
}

function buildLatexTable(entries: AnnotatedLine[]): string {
  const hasUnits = entries.some((e) => e.units);
  const cols = hasUnits ? 'lll' : 'll';
  const header = hasUnits ?
    '\\toprule\nName & Value & Units \\\\\n\\midrule' :
    '\\toprule\nName & Value \\\\\n\\midrule';

  const rows = entries.map((e) => {
    const name = `$${identifierToLatex(e.name)}$`;
    const value = `$${e.value}$`;
    if (hasUnits)
      return `${name} & ${value} & ${e.units || ''}`;
    return `${name} & ${value}`;
  });

  return `\\begin{tabular}{${cols}}\n${header}\n${rows.join(' \\\\\n')}\n\\bottomrule\n\\end{tabular}`;
}

function buildMarkdownTable(entries: AnnotatedLine[]): string {
  const hasUnits = entries.some((e) => e.units);

  let header: string;
  let separator: string;
  if (hasUnits) {
    header = '| Name | Value | Units |';
    separator = '| --- | --- | --- |';
  } else {
    header = '| Name | Value |';
    separator = '| --- | --- |';
  }

  const rows = entries.map((e) => {
    const name = `$${identifierToLatex(e.name)}$`;
    const value = `$${e.value}$`;
    if (hasUnits)
      return `| ${name} | ${value} | ${e.units || ''} |`;
    return `| ${name} | ${value} |`;
  });

  return `${header}\n${separator}\n${rows.join('\n')}`;
}
