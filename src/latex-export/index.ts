/* Converts IVP (Diff Studio) ODE model files to LaTeX / Markdown. */

import {ConvertOptions} from './types';
import {parseIvp} from './parser/ivp-parser';
import {buildLatexDocument, buildMarkdownDocument} from './generator/document-builder';

export type {ConvertOptions, ParsedModel, ASTNode, Token} from './types';

const DEFAULT_OPTIONS: ConvertOptions = {
  format: 'latex',
  includeMetadata: true,
  includeInits: true,
  includeParameters: true,
  includeConstants: true,
  useCdot: true,
  compact: false,
};

/** Convert an IVP model file content to LaTeX markup.
 *  @param ivpText  raw text content of an .ivp file
 *  @param options  conversion options
 *  @returns        LaTeX or Markdown string
 */
export function convertIvpToLatex(
  ivpText: string,
  options?: Partial<ConvertOptions>,
): string {
  const opts = {...DEFAULT_OPTIONS, ...options};
  const model = parseIvp(ivpText);

  if (opts.format === 'markdown')
    return buildMarkdownDocument(model, opts);
  return buildLatexDocument(model, opts);
}
