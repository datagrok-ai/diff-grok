// ============================================================
// latex-export — Public API
// ============================================================
// Part of the diff-grok library.
// Converts IVP (Diff Studio) ODE model files to LaTeX / Markdown.

import { ConvertOptions } from './types';
// import { parseIvp } from './parser/ivp-parser';
// import { buildLatexDocument } from './generator/document-builder';
// import { buildMarkdownDocument } from './generator/document-builder';

export type { ConvertOptions, ParsedModel, ASTNode, Token } from './types';

const DEFAULT_OPTIONS: ConvertOptions = {
  format: 'latex',
  includeMetadata: true,
  includeInits: true,
  includeParameters: true,
  includeConstants: true,
  useCdot: true,
};

/**
 * Convert an IVP model file content to LaTeX markup.
 *
 * @param ivpText - Raw text content of an .ivp file
 * @param options - Conversion options
 * @returns LaTeX or Markdown string
 */
export function convertIvpToLatex(
  ivpText: string,
  options?: Partial<ConvertOptions>,
): string {
  const opts = { ...DEFAULT_OPTIONS, ...options };

  // TODO: Implement — follow this pipeline:
  //
  // 1. Parse IVP text → ParsedModel
  //      const model = parseIvp(ivpText);
  //
  // 2. Build document (internally tokenizes, parses AST, generates LaTeX for each formula)
  //      if (opts.format === 'markdown')
  //        return buildMarkdownDocument(model, opts);
  //      else
  //        return buildLatexDocument(model, opts);

  throw new Error('Not implemented yet');
}
