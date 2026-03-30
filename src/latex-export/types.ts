/* IVP-to-LaTeX converter — type definitions. */

export type TokenType =
  | 'NUMBER'
  | 'IDENTIFIER'
  | 'OPERATOR'
  | 'LPAREN'
  | 'RPAREN'
  | 'COMPARISON'
  | 'QUESTION'
  | 'COLON'
  | 'COMMA'
  | 'ARROW'
  | 'ASSIGN'
  | 'SEMICOLON';

export interface Token {
  type: TokenType;
  value: string;
  position: number;
}

export interface NumberNode {
  type: 'number';
  value: string;
}

export interface IdentifierNode {
  type: 'identifier';
  name: string;
}

export interface BinaryNode {
  type: 'binary';
  op: string;
  left: ASTNode;
  right: ASTNode;
}

export interface UnaryNode {
  type: 'unary';
  op: string;
  operand: ASTNode;
}

export interface CallNode {
  type: 'call';
  name: string;
  args: ASTNode[];
}

export interface TernaryNode {
  type: 'ternary';
  condition: ASTNode;
  consequent: ASTNode;
  alternate: ASTNode;
}

export type ASTNode =
  | NumberNode
  | IdentifierNode
  | BinaryNode
  | UnaryNode
  | CallNode
  | TernaryNode;

export interface FormulaLine {
  /** Left-hand side: e.g. "dx/dt", "E1", "func1" */
  lhs: string;
  /** Right-hand side expression */
  rhs: string;
  /** True if LHS is a derivative: d(...)/dt or dX/dt */
  isDerivative: boolean;
  /** True if RHS contains => (arrow function) */
  isArrowFunction: boolean;
  /** Parameter names for arrow functions, e.g. ["p", "t"] */
  arrowParams?: string[];
}

export interface AnnotatedLine {
  name: string;
  value: string;
  units?: string;
  caption?: string;
  category?: string;
  tooltip?: string;
  min?: string;
  max?: string;
  step?: string;
}

export interface OutputLine {
  name: string;
  caption?: string;
}

export interface LoopBlock {
  entries: AnnotatedLine[];
}

export interface UpdateBlock {
  stageName: string;
  entries: AnnotatedLine[];
}

export interface ParsedModel {
  name?: string;
  description?: string;
  comment?: string;
  equations: FormulaLine[];
  expressions: FormulaLine[];
  inits: AnnotatedLine[];
  parameters: AnnotatedLine[];
  constants: AnnotatedLine[];
  argument: {
    name: string;
    stageName?: string;
    entries: AnnotatedLine[];
  };
  loops: LoopBlock[];
  updates: UpdateBlock[];
  output: OutputLine[];
  tolerance?: string;
  solverMeta?: string;
}

export interface ConvertOptions {
  /** Output format: 'latex' or 'markdown'. Default: 'latex' */
  format: 'latex' | 'markdown';

  /** Include metadata sections (name, description, comment). Default: true */
  includeMetadata: boolean;

  /** Include initial conditions table. Default: true */
  includeInits: boolean;

  /** Include parameters table. Default: true */
  includeParameters: boolean;

  /** Include constants table. Default: true */
  includeConstants: boolean;

  /** Use \cdot for multiplication. Default: true. If false, use juxtaposition. */
  useCdot: boolean;
}

export interface OperatorContext {
  useCdot?: boolean;
  baseIsCompound?: boolean;
  insideExponent?: boolean;
}

// Centralized error messages (cf. solver-defs.ts ERROR_MSG)
function tokenStr(tok?: {type: string; value: string}): string {
  return tok ? `${tok.type} '${tok.value}'` : 'EOF';
}

export const ERROR_MSG = {
  expectedToken: (type: string, value: string | undefined, tok?: {type: string; value: string}) =>
    `Expected ${type}${value ? ` '${value}'` : ''}, got ${tokenStr(tok)}`,
  unexpectedToken: (tok?: {type: string; value: string}) =>
    `Unexpected token: ${tokenStr(tok)}`,
};
