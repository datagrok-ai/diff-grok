/* LaTeX generator: recursively walks the AST and produces LaTeX strings. */

import {ASTNode} from '../types';
import {parseExpression} from '../parser/ast-parser';
import {identifierToLatex} from '../transformer/identifier';
import {functionToLatex} from '../transformer/functions';
import {operatorToLatex} from '../transformer/operators';
import {needsParentheses, wrapParens} from './bracket-manager';

/** Convert an expression string to LaTeX by parsing and rendering the AST.
 *  @param expr  expression string (RHS of a formula)
 *  @returns     LaTeX string
 */
export function expressionToLatex(expr: string): string {
  const ast = parseExpression(expr);
  return nodeToLatex(ast);
}

/** Render an AST node to LaTeX.
 *  @param node       AST node to render
 *  @param _parentOp  parent operator (unused, reserved for context)
 *  @param _position  child position (unused, reserved for context)
 *  @returns          LaTeX string
 */
export function nodeToLatex(
  node: ASTNode,
  _parentOp?: string,
  _position?: 'left' | 'right',
): string {
  switch (node.type) {
  case 'number':
    return numberToLatex(node.value);

  case 'identifier':
    return identifierToLatex(node.name);

  case 'unary': {
    const operand = nodeToLatex(node.operand);
    const needsWrap = node.operand.type === 'binary' || node.operand.type === 'ternary';
    return `${node.op}${needsWrap ? `\\left(${operand}\\right)` : operand}`;
  }

  case 'binary':
    return renderBinary(node);

  case 'call': {
    const args = node.args.map((a) => nodeToLatex(a));
    return functionToLatex(node.name, args);
  }

  case 'ternary': {
    const cond = nodeToLatex(node.condition);
    const cons = nodeToLatex(node.consequent);
    const alt = nodeToLatex(node.alternate);
    return `\\begin{cases} ${cons}, & \\text{if } ${cond} \\\\ ${alt}, & \\text{otherwise} \\end{cases}`;
  }
  }
}

function renderBinary(node: ASTNode & { type: 'binary' }): string {
  const {op, left, right} = node;

  // Division → \frac (no extra brackets needed)
  if (op === '/') {
    const l = nodeToLatex(left);
    const r = nodeToLatex(right);
    return operatorToLatex('/', l, r);
  }

  // Exponentiation → base^{exp}
  if (op === '**') {
    const isCompound = left.type === 'binary' || left.type === 'unary' || left.type === 'ternary';
    const l = nodeToLatex(left);
    const r = nodeToLatex(right);
    return operatorToLatex('**', l, r, {baseIsCompound: isCompound});
  }

  // For *, +, -, comparisons: use precedence-based bracketing
  const leftOp = getNodeOp(left);
  const rightOp = getNodeOp(right);

  const leftNeedsParen = needsParentheses(leftOp, op, 'left');
  const rightNeedsParen = needsParentheses(rightOp, op, 'right');

  const l = wrapParens(nodeToLatex(left), leftNeedsParen);
  const r = wrapParens(nodeToLatex(right), rightNeedsParen);

  if (op === '*')
    return operatorToLatex('*', l, r, {useCdot: true});

  return operatorToLatex(op, l, r);
}

function getNodeOp(node: ASTNode): string | undefined {
  if (node.type === 'binary') return node.op;
  if (node.type === 'unary') return 'unary';
  if (node.type === 'ternary') return '?';
  return undefined;
}

/** Convert a derivative LHS to LaTeX (e.g. "dy/dt" → "\frac{dy}{dt}").
 *  @param lhs  derivative string, e.g. "d(X)/dt" or "dX/dt"
 *  @returns    LaTeX string
 */
export function derivativeToLatex(lhs: string): string {
  // Pattern: d(varName)/darg or dvarName/darg
  const parenMatch = lhs.match(/^d\(([^)]+)\)\/d(.+)$/);
  if (parenMatch) {
    const varLatex = identifierToLatex(parenMatch[1]);
    const argLatex = parenMatch[2];
    return `\\frac{d\\,${varLatex}}{d${argLatex}}`;
  }

  const simpleMatch = lhs.match(/^d([A-Za-z]\w*)\/d(.+)$/);
  if (simpleMatch) {
    const varLatex = identifierToLatex(simpleMatch[1]);
    const argLatex = simpleMatch[2];
    return `\\frac{d${varLatex}}{d${argLatex}}`;
  }

  return lhs;
}

/** Render a number, converting scientific notation (e.g. "1e4" → "10^{4}").
 *  @param value  number string from the tokenizer
 *  @returns      LaTeX string
 */
export function numberToLatex(value: string): string {
  const match = value.match(/^([^eE]+)[eE]([+-]?\d+)$/);
  if (!match) return value;

  const coeff = match[1];
  let exp = match[2];
  if (exp.startsWith('+')) exp = exp.slice(1);

  if (coeff === '1') return `10^{${exp}}`;
  return `${coeff} \\times 10^{${exp}}`;
}
