// ============================================================
// Markdown Output Formatter
// ============================================================
// Ensures LaTeX math is properly wrapped for Markdown rendering.
// Uses $...$ for inline math and $$...$$ for display math.

/**
 * Wrap a LaTeX expression as inline math for Markdown.
 */
export function inlineMath(latex: string): string {
  return `$${latex}$`;
}

/**
 * Wrap a LaTeX block as display math for Markdown.
 */
export function displayMath(latex: string): string {
  return `$$\n${latex}\n$$`;
}

/**
 * Build a Markdown table from rows.
 *
 * @param headers - Column headers
 * @param rows - Array of row arrays
 * @returns Markdown table string
 */
export function markdownTable(headers: string[], rows: string[][]): string {
  const headerLine = '| ' + headers.join(' | ') + ' |';
  const separatorLine = '| ' + headers.map(() => '---').join(' | ') + ' |';
  const dataLines = rows.map(row => '| ' + row.join(' | ') + ' |');
  return [headerLine, separatorLine, ...dataLines].join('\n');
}
