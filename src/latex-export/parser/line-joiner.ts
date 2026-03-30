// ============================================================
// Line Preprocessor
// ============================================================
// Responsibilities:
//   1. Strip // comments (everything from // to end of line)
//   2. Strip {annotations} and [tooltips] (but parse them first for metadata)
//   3. Join multi-line formulas (continuation lines without top-level =)

/**
 * Remove everything from `//` to end of line.
 * A single `/` is division — only `//` is a comment.
 * Trims trailing whitespace after removal.
 */
export function stripComments(line: string): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Remove `{...}` annotation blocks and `[...]` tooltip blocks from a line.
 * Returns the line with only the `name = value` portion.
 * Trims trailing whitespace.
 */
export function stripAnnotations(line: string): string {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Parse `{key1: val1; key2: val2}` annotation string into a Record.
 * Returns an object with keys: units, caption, category, min, max, step.
 */
export function parseAnnotations(annotation: string): Record<string, string> {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Extract tooltip text from `[...]` at the end of a line.
 */
export function parseTooltip(line: string): string | undefined {
  // TODO: Implement
  throw new Error('Not implemented');
}

/**
 * Join multi-line formulas. A continuation line is:
 *   - Indented (starts with whitespace)
 *   - Does NOT contain `=` at the top level (outside parentheses)
 *   - Does NOT start with a derivative pattern `d.../dt`
 *
 * @param lines - Pre-stripped lines (comments already removed)
 * @returns Array of joined formula strings
 */
export function joinMultiLineFormulas(lines: string[]): string[] {
  // TODO: Implement
  throw new Error('Not implemented');
}
