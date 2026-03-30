// ============================================================
// TeX Output Formatter
// ============================================================
// Wraps the document content for standalone .tex file usage.
// Uses amsmath, booktabs packages.

/**
 * Wrap LaTeX content in a minimal document preamble.
 * Useful if the user wants a self-contained .tex file.
 *
 * @param content - LaTeX body content
 * @returns Complete .tex file content
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
