// ============================================================
// IVP File Parser
// ============================================================
// Parses the full text of an .ivp file into a structured ParsedModel.
//
// Block detection:
//   Lines starting with `#keyword` begin a new block.
//   Content lines are indented with spaces.
//
// See docs/SYNTAX.md for the complete format reference.

import { ParsedModel } from '../types';

/**
 * Parse the full text of an IVP file into a structured model.
 *
 * Steps:
 *   1. Split into lines
 *   2. Identify block boundaries (#name, #equations, etc.)
 *   3. For each block, preprocess lines (strip comments, join multi-line)
 *   4. Parse block-specific content (equations, parameters, etc.)
 *
 * @param ivpText - Raw content of an .ivp file
 * @returns Parsed model structure
 */
export function parseIvp(ivpText: string): ParsedModel {
  // TODO: Implement
  throw new Error('Not implemented');
}
