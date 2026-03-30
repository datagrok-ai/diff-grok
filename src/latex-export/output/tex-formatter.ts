/* TeX output formatter: wraps content in a document preamble for standalone .tex files. */

/** Wrap LaTeX content in a minimal document preamble for a self-contained .tex file.
 *  @param content  LaTeX body content
 *  @returns        complete .tex file content
 */
export function wrapTexDocument(content: string): string {
  return [
    '\\documentclass{article}',
    '\\usepackage{amsmath}',
    '\\usepackage{amssymb}',
    '\\usepackage{booktabs}',
    '',
    '\\begin{document}',
    '',
    content,
    '',
    '\\end{document}',
  ].join('\n');
}
