/* Document builder: assembles converted blocks into a complete LaTeX or Markdown document. */

import {ParsedModel, ConvertOptions, FormulaLine, AnnotatedLine} from '../types';
import {expressionToLatex, derivativeToLatex, RenderOptions} from './latex-generator';
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
    parts.push(buildLatexAlign(model.equations, opts));
    if (model.argument.entries.length >= 3)
      parts.push(buildLatexArgRange(model));
  }

  if (model.expressions.length > 0) {
    parts.push(`\\subsection{${plural('Expression', model.expressions.length)}}`);
    parts.push(buildLatexAlign(model.expressions, opts));
  }

  const loopInfo = extractLoopInfo(model);
  if (loopInfo && loopInfo.updates.length > 0) {
    parts.push('\\subsection{Cyclic Update}');
    const lines = loopInfo.updates.map((u) => `  ${renderLoopUpdateLine(u, opts)}`);
    parts.push(`Before each cycle:\n\\begin{align}\n${lines.join(' \\\\\n')}\n\\end{align}`);
  }

  const stages = extractStageInfos(model);
  if (stages.length > 0) {
    parts.push('\\subsection{Stage Transitions}');
    for (const stage of stages) {
      const dur = stage.duration ? ` (duration: $${expressionToLatex(stage.duration, opts)}$)` : '';
      parts.push(`\\textbf{${stage.stageName}}${dur}`);
      if (stage.transitions.length > 0) {
        const lines = stage.transitions.map((t) => `  ${renderTransitionLine(t, opts)}`);
        parts.push(`Before this stage:\n\\begin{align}\n${lines.join(' \\\\\n')}\n\\end{align}`);
      }
    }
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
    parts.push(buildMarkdownAlign(model.equations, opts));
    if (model.argument.entries.length >= 3)
      parts.push(buildMarkdownArgRange(model));
  }

  if (model.expressions.length > 0) {
    parts.push(`### ${plural('Expression', model.expressions.length)}`);
    parts.push(buildMarkdownAlign(model.expressions, opts));
  }

  const loopInfo = extractLoopInfo(model);
  if (loopInfo && loopInfo.updates.length > 0) {
    parts.push('### Cyclic Update');
    const lines = loopInfo.updates.map((u) => `  ${renderLoopUpdateLine(u, opts)}`);
    parts.push(`Before each cycle:\n\n$$\n\\begin{aligned}\n${lines.join(' \\\\\n')}\n\\end{aligned}\n$$`);
  }

  const stages = extractStageInfos(model);
  if (stages.length > 0) {
    parts.push('### Stage Transitions');
    for (const stage of stages) {
      const dur = stage.duration ? ` (duration: $${expressionToLatex(stage.duration, opts)}$)` : '';
      parts.push(`**${stage.stageName}**${dur}`);
      if (stage.transitions.length > 0) {
        const lines = stage.transitions.map((t) => `  ${renderTransitionLine(t, opts)}`);
        parts.push(`Before this stage:\n\n$$\n\\begin{aligned}\n${lines.join(' \\\\\n')}\n\\end{aligned}\n$$`);
      }
    }
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

  const loopInfo = extractLoopInfo(model);
  const stages = extractStageInfos(model);
  let effectiveFinish = finish;
  let suffix = '';

  if (loopInfo) {
    // Cyclic model: total = start + count * (finish - start)
    const s = parseFloat(start);
    const f = parseFloat(finish);
    const n = parseFloat(loopInfo.count);
    if (!isNaN(s) && !isNaN(f) && !isNaN(n))
      effectiveFinish = String(s + n * (f - s));
    else
      effectiveFinish = `${start} + ${loopInfo.count} \\cdot (${finish} - ${start})`;
    const duration = (!isNaN(parseFloat(start)) && !isNaN(parseFloat(finish))) ?
      String(parseFloat(finish) - parseFloat(start)) : `${finish} - ${start}`;
    suffix = `, \\quad ${arg}_{\\text{cycle}} = ${duration}`;
  } else if (stages.length > 0) {
    // Multistage model: total = start + stage1_duration + stage2_duration + ...
    const durations = [finish, ...stages.map((s) => s.duration).filter(Boolean) as string[]];
    const totalExpr = `${start} + ${durations.join(' + ')}`;
    const total = tryEvalNumeric(totalExpr, model);
    if (total !== undefined)
      effectiveFinish = String(total);
    else
      effectiveFinish = totalExpr;
  }

  return `${arg} \\in \\left[${start},\\, ${effectiveFinish}\\right], \\quad \\Delta ${arg} = ${step}${suffix}`;
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

function buildLatexAlign(formulas: FormulaLine[], opts: RenderOptions): string {
  const lines = formulas.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs, opts);
    return `  ${lhs} &= ${rhs}`;
  });
  return `\\begin{align}\n${lines.join(' \\\\\n')}\n\\end{align}`;
}

function buildMarkdownAlign(formulas: FormulaLine[], opts: RenderOptions): string {
  const lines = formulas.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs, opts);
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

function compactExprItems(model: ParsedModel, opts: RenderOptions): string[] {
  return model.expressions.map((f) => {
    const lhs = renderFormulaLhs(f);
    const rhs = expressionToLatex(f.rhs, opts);
    return `${lhs} = ${rhs}`;
  });
}

interface LoopUpdate { variable: string; value: string }

interface LoopInfo { count: string; updates: LoopUpdate[] }

function extractLoopInfo(model: ParsedModel): LoopInfo | undefined {
  if (model.loops.length === 0) return undefined;
  const loop = model.loops[0];
  let count = '';
  const updates: LoopUpdate[] = [];
  for (const entry of loop.entries) {
    if (entry.name === 'count')
      count = entry.value;
    else {
      const variable = entry.name.replace(/\s*\+$/, '');
      updates.push({variable, value: entry.value});
    }
  }
  return count ? {count, updates} : undefined;
}

function renderLoopUpdateLine(u: LoopUpdate, opts: RenderOptions): string {
  const varLatex = identifierToLatex(u.variable);
  const valLatex = expressionToLatex(u.value, opts);
  return `${varLatex} \\leftarrow ${varLatex} + ${valLatex}`;
}

interface StageTransition {
  variable: string;
  value: string;
  isIncrement: boolean;
}

interface StageInfo {
  stageName: string;
  duration?: string;
  transitions: StageTransition[];
}

/** Replace internal argument refs (_t0, _t1, _h) with actual values from #argument. */
function resolveArgRefs(expr: string, model: ParsedModel): string {
  const entries = model.argument.entries;
  let result = expr;
  if (entries.length > 0) result = result.replace(/\b_t0\b/g, entries[0].value);
  if (entries.length > 1) result = result.replace(/\b_t1\b/g, entries[1].value);
  if (entries.length > 2) result = result.replace(/\b_h\b/g, entries[2].value);
  return result;
}

/** Try to evaluate a simple arithmetic expression numerically, resolving parameter/constant names. */
function tryEvalNumeric(expr: string, model: ParsedModel): number | undefined {
  let resolved = expr;
  for (const p of [...model.parameters, ...model.constants])
    resolved = resolved.replace(new RegExp(`\\b${p.name}\\b`, 'g'), p.value);
  const val = parseFloat(resolved);
  if (!isNaN(val) && /^[\d.eE+\-\s*/()]+$/.test(resolved.trim())) {
    try { return Function(`"use strict"; return (${resolved})`)() as number; } catch { /* skip */ }
  }
  return undefined;
}

function extractStageInfos(model: ParsedModel): StageInfo[] {
  return model.updates.map((block) => {
    let duration: string | undefined;
    const transitions: StageTransition[] = [];
    for (const entry of block.entries) {
      if (entry.name === 'duration') {
        duration = resolveArgRefs(entry.value, model);
      } else {
        const isIncrement = entry.name.endsWith('+');
        const variable = isIncrement ? entry.name.replace(/\s*\+$/, '') : entry.name;
        transitions.push({variable, value: entry.value, isIncrement});
      }
    }
    return {stageName: block.stageName, duration, transitions};
  });
}

function renderTransitionLine(t: StageTransition, opts: RenderOptions): string {
  const varLatex = identifierToLatex(t.variable);
  const valLatex = expressionToLatex(t.value, opts);
  if (t.isIncrement)
    return `${varLatex} \\leftarrow ${varLatex} + ${valLatex}`;
  return `${varLatex} \\leftarrow ${valLatex}`;
}

function buildLatexCompact(model: ParsedModel, opts: ConvertOptions): string {
  const parts: string[] = [];

  if (model.name)
    parts.push(`\\textbf{${model.name}}`);

  if (model.equations.length > 0) {
    parts.push(`${plural('Equation', model.equations.length)}:`);
    parts.push(buildLatexAlign(model.equations, opts));
    if (model.argument.entries.length >= 3)
      parts.push(buildLatexArgRange(model));
  }

  if (opts.includeInits && model.inits.length > 0) {
    const lines = compactInlineItems(model).map((item) => `\\[ ${item} \\]`);
    parts.push(`${compactInitsLabel(model)}:\n${lines.join('\n')}`);
  }

  if (model.expressions.length > 0) {
    const lines = compactExprItems(model, opts).map((item) => `\\[ ${item} \\]`);
    parts.push(`where\n${lines.join('\n')}`);
  }

  const loopInfo = extractLoopInfo(model);
  if (loopInfo && loopInfo.updates.length > 0) {
    const lines = loopInfo.updates.map((u) => `\\[ ${renderLoopUpdateLine(u, opts)} \\]`);
    parts.push(`before each cycle:\n${lines.join('\n')}`);
  }

  const stages = extractStageInfos(model);
  for (const stage of stages) {
    if (stage.transitions.length > 0) {
      const dur = stage.duration ? ` ($${expressionToLatex(stage.duration, opts)}$)` : '';
      const lines = stage.transitions.map((t) => `\\[ ${renderTransitionLine(t, opts)} \\]`);
      parts.push(`before ${stage.stageName}${dur}:\n${lines.join('\n')}`);
    }
  }

  if (opts.includeParameters && model.parameters.length > 0)
    parts.push(`${plural('Parameter', model.parameters.length)}:\n${buildLatexList(model.parameters)}`);

  if (opts.includeConstants && model.constants.length > 0)
    parts.push(`${plural('Constant', model.constants.length)}:\n${buildLatexList(model.constants)}`);

  return parts.join('\n\n');
}

function buildMarkdownCompact(model: ParsedModel, opts: ConvertOptions): string {
  const parts: string[] = [];

  if (model.name)
    parts.push(`**${model.name}**`);

  if (model.equations.length > 0) {
    parts.push(`${plural('Equation', model.equations.length)}:`);
    parts.push(buildMarkdownAlign(model.equations, opts));
    if (model.argument.entries.length >= 3)
      parts.push(buildMarkdownArgRange(model));
  }

  if (opts.includeInits && model.inits.length > 0) {
    const lines = compactInlineItems(model).map((item) => `$$${item}$$`);
    parts.push(`${compactInitsLabel(model)}:\n${lines.join('\n')}`);
  }

  if (model.expressions.length > 0) {
    const lines = compactExprItems(model, opts).map((item) => `$$${item}$$`);
    parts.push(`where\n${lines.join('\n')}`);
  }

  const loopInfo = extractLoopInfo(model);
  if (loopInfo && loopInfo.updates.length > 0) {
    const lines = loopInfo.updates.map((u) => `$$${renderLoopUpdateLine(u, opts)}$$`);
    parts.push(`before each cycle:\n${lines.join('\n')}`);
  }

  const stages = extractStageInfos(model);
  for (const stage of stages) {
    if (stage.transitions.length > 0) {
      const dur = stage.duration ? ` ($${expressionToLatex(stage.duration, opts)}$)` : '';
      const lines = stage.transitions.map((t) => `$$${renderTransitionLine(t, opts)}$$`);
      parts.push(`before ${stage.stageName}${dur}:\n${lines.join('\n')}`);
    }
  }

  if (opts.includeParameters && model.parameters.length > 0)
    parts.push(`${plural('Parameter', model.parameters.length)}:\n${buildMarkdownList(model.parameters)}`);

  if (opts.includeConstants && model.constants.length > 0)
    parts.push(`${plural('Constant', model.constants.length)}:\n${buildMarkdownList(model.constants)}`);

  return parts.join('\n\n');
}

function buildLatexList(entries: AnnotatedLine[]): string {
  const items = entries.map((e) => {
    const name = `$${identifierToLatex(e.name)}$`;
    const value = `$${e.value}$`;
    const units = e.units ? ` (${e.units})` : '';
    return `  \\item ${name} = ${value}${units}`;
  }).join('\n');
  return `\\begin{itemize}\n${items}\n\\end{itemize}`;
}

function buildMarkdownList(entries: AnnotatedLine[]): string {
  return entries.map((e) => {
    const name = `$${identifierToLatex(e.name)}$`;
    const value = `$${e.value}$`;
    const units = e.units ? ` (${e.units})` : '';
    return `- ${name} = ${value}${units}`;
  }).join('\n');
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
