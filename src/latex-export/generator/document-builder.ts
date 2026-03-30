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
  compact: false,
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
  if (opts.compact)
    return buildLatexCompact(model, opts);

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
    parts.push(`\\subsection{${plural('Equation', model.equations.length)}}`);
    parts.push(buildLatexAlign(model.equations));
    if (model.argument.entries.length >= 3)
      parts.push(buildLatexArgRange(model));
  }

  if (model.expressions.length > 0) {
    parts.push(`\\subsection{${plural('Expression', model.expressions.length)}}`);
    parts.push(buildLatexAlign(model.expressions));
  }

  if (opts.includeInits && model.inits.length > 0) {
    parts.push(`\\subsection{${initsTitle(model)}}`);
    parts.push(buildLatexTable(model.inits));
  }

  if (opts.includeParameters && model.parameters.length > 0) {
    parts.push(`\\subsection{${plural('Parameter', model.parameters.length)}}`);
    parts.push(buildLatexTable(model.parameters));
  }

  if (opts.includeConstants && model.constants.length > 0) {
    parts.push(`\\subsection{${plural('Constant', model.constants.length)}}`);
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
  if (opts.compact)
    return buildMarkdownCompact(model, opts);

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
    parts.push(`### ${plural('Equation', model.equations.length)}`);
    parts.push(buildMarkdownAlign(model.equations));
    if (model.argument.entries.length >= 3)
      parts.push(buildMarkdownArgRange(model));
  }

  if (model.expressions.length > 0) {
    parts.push(`### ${plural('Expression', model.expressions.length)}`);
    parts.push(buildMarkdownAlign(model.expressions));
  }

  if (opts.includeInits && model.inits.length > 0) {
    parts.push(`### ${initsTitle(model)}`);
    parts.push(buildMarkdownTable(model.inits));
  }

  if (opts.includeParameters && model.parameters.length > 0) {
    parts.push(`### ${plural('Parameter', model.parameters.length)}`);
    parts.push(buildMarkdownTable(model.parameters));
  }

  if (opts.includeConstants && model.constants.length > 0) {
    parts.push(`### ${plural('Constant', model.constants.length)}`);
    parts.push(buildMarkdownTable(model.constants));
  }

  return parts.join('\n\n');
}

function plural(singular: string, count: number): string {
  if (count === 1) return singular;
  if (singular === 'Constant') return 'Constants';
  return singular + 's';
}

function buildArgRangeLatex(model: ParsedModel): string {
  const arg = model.argument.name;
  const start = model.argument.entries[0].value;
  const finish = model.argument.entries[1].value;
  const step = model.argument.entries[2].value;
  return `${arg} \\in \\left[${start},\\, ${finish}\\right], \\quad \\Delta ${arg} = ${step}`;
}

function buildLatexArgRange(model: ParsedModel): string {
  return `\\[\n  ${buildArgRangeLatex(model)}\n\\]`;
}

function buildMarkdownArgRange(model: ParsedModel): string {
  return `$$${buildArgRangeLatex(model)}$$`;
}

function initsTitle(model: ParsedModel): string {
  const base = plural('Initial Condition', model.inits.length);
  const argName = model.argument.name;
  const initEntry = model.argument.entries[0];
  if (argName && initEntry)
    return `${base} (${argName} = ${initEntry.value})`;
  return base;
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

function buildInlineInits(model: ParsedModel): string {
  const argName = model.argument.name;
  const t0 = model.argument.entries.length > 0 ? model.argument.entries[0].value : '0';
  const items = model.inits.map((e) =>
    `$${identifierToLatex(e.name)}(${t0}) = ${e.value}$`);
  return items.join(', ');
}

function compactInitsLabel(model: ParsedModel): string {
  return plural('initial condition', model.inits.length);
}

function compactInlineItems(model: ParsedModel): string[] {
  const t0 = model.argument.entries.length > 0 ? model.argument.entries[0].value : '0';
  return model.inits.map((e) => {
    const name = identifierToLatex(e.name);
    return `${name}(${t0}) = ${e.value}`;
  });
}

function compactExprItems(model: ParsedModel): string[] {
  return model.expressions.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs);
    return `${lhs} = ${rhs}`;
  });
}

function buildLatexCompact(model: ParsedModel, opts: ConvertOptions): string {
  const parts: string[] = [];

  if (model.name)
    parts.push(`\\textbf{${model.name}}`);

  if (model.equations.length > 0) {
    parts.push(`${plural('Equation', model.equations.length)}:`);
    parts.push(buildLatexAlign(model.equations));
    if (model.argument.entries.length >= 3)
      parts.push(buildLatexArgRange(model));
  }

  if (opts.includeInits && model.inits.length > 0) {
    const lines = compactInlineItems(model).map((item) => `\\[ ${item} \\]`);
    parts.push(`${compactInitsLabel(model)}:\n${lines.join('\n')}`);
  }

  if (model.expressions.length > 0) {
    const lines = compactExprItems(model).map((item) => `\\[ ${item} \\]`);
    parts.push(`where\n${lines.join('\n')}`);
  }

  if (opts.includeParameters && model.parameters.length > 0)
    parts.push(`${plural('Parameter', model.parameters.length)}:\n\n${buildLatexTable(model.parameters)}`);

  if (opts.includeConstants && model.constants.length > 0)
    parts.push(`${plural('Constant', model.constants.length)}:\n\n${buildLatexTable(model.constants)}`);

  return parts.join('\n\n');
}

function buildMarkdownCompact(model: ParsedModel, opts: ConvertOptions): string {
  const parts: string[] = [];

  if (model.name)
    parts.push(`**${model.name}**`);

  if (model.equations.length > 0) {
    parts.push(`${plural('Equation', model.equations.length)}:`);
    parts.push(buildMarkdownAlign(model.equations));
    if (model.argument.entries.length >= 3)
      parts.push(buildMarkdownArgRange(model));
  }

  if (opts.includeInits && model.inits.length > 0) {
    const lines = compactInlineItems(model).map((item) => `$$${item}$$`);
    parts.push(`${compactInitsLabel(model)}:\n${lines.join('\n')}`);
  }

  if (model.expressions.length > 0) {
    const lines = compactExprItems(model).map((item) => `$$${item}$$`);
    parts.push(`where\n${lines.join('\n')}`);
  }

  if (opts.includeParameters && model.parameters.length > 0)
    parts.push(`${plural('Parameter', model.parameters.length)}:\n\n${buildMarkdownTable(model.parameters)}`);

  if (opts.includeConstants && model.constants.length > 0)
    parts.push(`${plural('Constant', model.constants.length)}:\n\n${buildMarkdownTable(model.constants)}`);

  return parts.join('\n\n');
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

  return `\\begin{tabular}{${cols}}\n${header}\n${rows.join(' \\\\\n')} \\\\\n\\bottomrule\n\\end{tabular}`;
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
