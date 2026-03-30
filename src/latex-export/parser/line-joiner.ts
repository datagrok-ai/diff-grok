/* Line preprocessor: strips comments, annotations, and joins multi-line formulas. */

/** Remove everything from `//` to end of line.
 *  @param line  raw line of text
 *  @returns     line with comment removed
 */
export function stripComments(line: string): string {
  const idx = line.indexOf('//');
  if (idx === -1) return line;
  return line.substring(0, idx).trimEnd();
}

/** Remove `{...}` annotation blocks and `[...]` tooltip blocks from a line.
 *  @param line  line containing annotations
 *  @returns     line with only the name = value portion
 */
export function stripAnnotations(line: string): string {
  return line.replace(/\s*\{[^}]*\}/, '').replace(/\s*\[[^\]]*\]/, '').trimEnd();
}

/** Parse `{key1: val1; key2: val2}` annotation string into a Record.
 *  @param annotation  annotation string including braces
 *  @returns           key-value pairs (units, caption, category, min, max, step)
 */
export function parseAnnotations(annotation: string): Record<string, string> {
  const result: Record<string, string> = {};
  // Remove surrounding braces
  const inner = annotation.replace(/^\{|\}$/g, '').trim();
  if (!inner) return result;
  const pairs = inner.split(';');
  for (const pair of pairs) {
    const colonIdx = pair.indexOf(':');
    if (colonIdx === -1) continue;
    const key = pair.substring(0, colonIdx).trim();
    const val = pair.substring(colonIdx + 1).trim();
    if (key) result[key] = val;
  }
  return result;
}

/** Extract tooltip text from `[...]` at the end of a line.
 *  @param line  line that may contain a tooltip
 *  @returns     tooltip text or undefined
 */
export function parseTooltip(line: string): string | undefined {
  const match = line.match(/\[([^\]]*)\]\s*$/);
  return match ? match[1] : undefined;
}

/** Join multi-line formulas (continuation lines without top-level `=`).
 *  @param lines  pre-stripped lines (comments already removed)
 *  @returns      array of joined formula strings
 */
export function joinMultiLineFormulas(lines: string[]): string[] {
  const result: string[] = [];
  for (const line of lines) {
    if (line.trim() === '') continue;
    const isIndented = /^\s/.test(line);
    const hasTopLevelEquals = hasTopLevelEqual(line);
    if (isIndented && !hasTopLevelEquals && result.length > 0)
      result[result.length - 1] += ' ' + line.trim();
    else
      result.push(line);
  }
  return result;
}

/** Check if a line contains `=` at the top level (outside parentheses). */
function hasTopLevelEqual(line: string): boolean {
  let depth = 0;
  for (let i = 0; i < line.length; i++) {
    if (line[i] === '(') depth++;
    else if (line[i] === ')') depth--;
    else if (line[i] === '=' && depth === 0) return true;
  }
  return false;
}
