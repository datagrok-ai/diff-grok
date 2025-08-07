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
  DEQS = 'Differential Equations',
  INITS = 'Initial Conditions',
  EXPR = 'Auxiliary Computations',
  PARAMS = 'Parameters',
  CONSTS = 'Constants',
  IVP = 'Initial Value Problem',
};

export enum EQN_TAG {
  DOLLAR,
  DOUBLE_DOLLAR,
};

export const SHORT_EQN_TAG = EQN_TAG.DOUBLE_DOLLAR;

export type FormattingOptions = {
  eqnTag: EQN_TAG,
  derForm: DER_FORM,
};

export type EqnTags = {
  open: string,
  close: string,
};

export function getOpenCloseEqTags(type?: EQN_TAG): EqnTags {
  switch (type) {
  case EQN_TAG.DOLLAR:
    return {open: '$', close: '$'};

  case EQN_TAG.DOUBLE_DOLLAR:
    return {open: '$$', close: '$$'};

  default:
    return {open: '$$', close: '$$'};
  }
}
