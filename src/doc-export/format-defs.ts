export type DerivativeFormatting = (argName: string, funcName: string) => string;

export const fracDerForm: DerivativeFormatting = (argName: string, funcName: string) => {
  return `\\frac{d${funcName}}{d${argName}}`;
};

export const dashDerForm: DerivativeFormatting = (argName: string, funcName: string) => {
  return `${funcName}'(${argName})`;
};

export const divDerForm: DerivativeFormatting = (argName: string, funcName: string) => {
  return `d${funcName}/d${argName}`;
};

export enum DER_FORM {
  FRAC,
  DASH,
  DIV,
};

export function getDerForm(type?: DER_FORM): DerivativeFormatting {
  switch (type) {
  case DER_FORM.FRAC:
    return fracDerForm;

  case DER_FORM.DASH:
    return dashDerForm;

  case DER_FORM.DIV:
    return divDerForm;

  default:
    return divDerForm;
  }
}

export enum HEADER_TAG {
  SECTION = '##',
  SUBSECTION = '####',
};

export enum TITLE {
  DEQS = 'Differential equations',
  INITS = 'Initial Conditions',
  EXPR = 'Auxiliary Computations',
  PARAMS = 'Parameters',
  CONSTS = 'Constants',
};

export enum FORMULA_TAG {
  DOLLAR,
  DOUBLE_DOLLAR,
};

export type MathMLOptions = {
  toTagGreeks: boolean,
  toReplaceAsterixWithCdot: boolean,
  toTagMathFuncs: boolean,
};

export type FormattingOptions = {
  eqnTag: FORMULA_TAG,
  valTag: FORMULA_TAG,
  derForm: DER_FORM,
  mathMlOpts: MathMLOptions,
};

export type EqnTags = {
  open: string,
  close: string,
};

export function getOpenCloseEqTags(type?: FORMULA_TAG): EqnTags {
  switch (type) {
  case FORMULA_TAG.DOLLAR:
    return {open: '$', close: '$'};

  case FORMULA_TAG.DOUBLE_DOLLAR:
    return {open: '$$', close: '$$'};

  default:
    return {open: '$$', close: '$$'};
  }
}

export const BASE_GREEK_LETTERS = [
  'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
  'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'omicron', 'pi', 'rho',
  'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega',
  'varepsilon', 'vartheta', 'varpi', 'varrho', 'varsigma', 'varphi',
];

export const GREEKS_TAGGING_DEFAULT = true;
export const ASTERIX_REPLACING_DEFAULT = true;

export const MATH_FUNCS = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh'];

export const MATH_FUNC_REPLACING_DEFAULT = true;
