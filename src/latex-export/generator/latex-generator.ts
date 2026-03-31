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

  // Division → \frac; flatten into product of fractions when chain has multiple divisions
  if (op === '/') {
    const factors = flattenMulDivChain(node);
    const divisorCount = factors.filter((f) => f.isDivisor).length;
    if (divisorCount > 1)
      return renderProductOfFractions(factors);
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

/** Flatten a left-associative * / chain into factors with divisor flags. */
function flattenMulDivChain(node: ASTNode & {type: 'binary'}):
  Array<{node: ASTNode; isDivisor: boolean}> {
  const factors: Array<{node: ASTNode; isDivisor: boolean}> = [];
  function collect(n: ASTNode): void {
    if (n.type === 'binary' && (n.op === '*' || n.op === '/')) {
      collect(n.left);
      factors.push({node: n.right, isDivisor: n.op === '/'});
    } else
      factors.push({node: n, isDivisor: false});
  }
  collect(node);
  return factors;
}

/** Render a list of factors joined with \cdot, adding parens where needed. */
function renderProductTerm(nodes: ASTNode[]): string {
  if (nodes.length === 0) return '1';
  if (nodes.length === 1) return nodeToLatex(nodes[0]);
  return nodes.map((n) => {
    const latex = nodeToLatex(n);
    if (n.type === 'binary' && (n.op === '+' || n.op === '-'))
      return `\\left(${latex}\\right)`;
    if (n.type === 'ternary')
      return `\\left(${latex}\\right)`;
    return latex;
  }).join(' \\cdot ');
}

/** Render factors as a product of fractions. */
function renderProductOfFractions(
  factors: Array<{node: ASTNode; isDivisor: boolean}>,
): string {
  const groups: Array<{nums: ASTNode[]; dens: ASTNode[]}> = [];
  let cur: {nums: ASTNode[]; dens: ASTNode[]} = {nums: [], dens: []};

  for (const f of factors) {
    if (f.isDivisor) {
      cur.dens.push(f.node);
    } else {
      if (cur.dens.length > 0) {
        groups.push(cur);
        cur = {nums: [], dens: []};
      }
      cur.nums.push(f.node);
    }
  }
  groups.push(cur);

  const parts = groups.map((g) => {
    const num = renderProductTerm(g.nums);
    if (g.dens.length === 0) return num;
    const den = renderProductTerm(g.dens);
    return `\\frac{${num}}{${den}}`;
  });

  return parts.join(' \\cdot ');
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
