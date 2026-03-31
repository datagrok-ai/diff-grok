/* Markdown output formatter: wraps LaTeX math for Markdown rendering. */

/** Wrap a LaTeX expression as inline math for Markdown.
 *  @param latex  LaTeX expression
 *  @returns      string wrapped in $...$
 */
export function inlineMath(latex: string): string {
  return `$${latex}$`;
}

/** Wrap a LaTeX block as display math for Markdown.
 *  @param latex  LaTeX block
 *  @returns      string wrapped in $$...$$
 */
export function displayMath(latex: string): string {
  return `$$\n${latex}\n$$`;
}

/** Build a Markdown table from rows.
 *  @param headers  column headers
 *  @param rows     array of row arrays
 *  @returns        Markdown table string
 */
export function markdownTable(headers: string[], rows: string[][]): string {
  const headerLine = '| ' + headers.join(' | ') + ' |';
  const separatorLine = '| ' + headers.map(() => '---').join(' | ') + ' |';
  const dataLines = rows.map((row) => '| ' + row.join(' | ') + ' |');
  return [headerLine, separatorLine, ...dataLines].join('\n');
}
