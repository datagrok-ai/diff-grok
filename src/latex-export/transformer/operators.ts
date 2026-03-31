/* Operator-to-LaTeX transformer: renders binary and comparison operators as LaTeX. */

import {OperatorContext} from '../types';

/** Render a binary operator with its operands as LaTeX.
 *  @param op       operator string (+, -, *, /, **, <, >, <=, >=, ==, !=)
 *  @param left     already-rendered LaTeX for the left operand
 *  @param right    already-rendered LaTeX for the right operand
 *  @param context  optional context flags
 *  @returns        LaTeX string
 */
export function operatorToLatex(
  op: string,
  left: string,
  right: string,
  context?: OperatorContext,
): string {
  switch (op) {
  case '+':
  case '-':
    return `${left} ${op} ${right}`;

  case '*':
    if (context?.useCdot === false)
      return `${left} \\, ${right}`;
    return `${left} \\cdot ${right}`;

  case '/':
    if (context?.insideExponent)
      return `${left} / ${right}`;
    return `\\frac{${left}}{${right}}`;

  case '**': {
    const base = context?.baseIsCompound ? `\\left(${left}\\right)` : left;
    return `${base}^{${right}}`;
  }

  case '>=':
    return `${left} \\geq ${right}`;
  case '<=':
    return `${left} \\leq ${right}`;
  case '!=':
    return `${left} \\neq ${right}`;
  case '==':
    return `${left} = ${right}`;

  case '<':
  case '>':
    return `${left} ${op} ${right}`;

  default:
    return `${left} ${op} ${right}`;
  }
}
