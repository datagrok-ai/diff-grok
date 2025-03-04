export const ARG_INP_COUNT = 3;
export const ARG_COL_IDX = 0;
export const CORRECTION_FACTOR = 0.5;
export const DEFAULT_CORRECTION = 0.0000001;

/** Default argument names*/
export enum ARG {
  START = '_t0',
  FINISH = '_t1',
  STEP = '_h',
};

export const argName2IdxMap = new Map<string, number>([
  [ARG.START, 0],
  [ARG.FINISH, 1],
  [ARG.STEP, 2],
]);
