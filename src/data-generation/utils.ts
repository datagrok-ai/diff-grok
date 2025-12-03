import {ANNOT_SEPAR, BRACE_OPEN, BRACE_CLOSE, CONTROL_SEP, IVP, Input} from '../scripting-tools';

export type InputOpts = {
  caption?: string,
  min?: number,
  max?: number,
  defaultVal: number,
};

export type ModelInputs = Record<string, number>;

/** Return options of IVP inputs */
export function getInputOptionsMap(ivp: IVP): Map<string, InputOpts> {
  const opts = new Map<string, InputOpts>([
    ['_t0', getInputOpts(ivp.arg.initial)],
    ['_t1', getInputOpts(ivp.arg.final)],
    ['_h', getInputOpts(ivp.arg.step)],
  ]);

  ivp.inits.forEach((input, name) => opts.set(name, getInputOpts(input)));
  ivp.params?.forEach((input, name) => opts.set(name, getInputOpts(input)));

  if (ivp.loop != null)
    opts.set('_count', getInputOpts(ivp.loop.count));

  return opts;
}

function getInputOpts(input: Input): InputOpts {
  const annot = input.annot;
  const defaultVal = input.value;
  const opts: InputOpts = {defaultVal: defaultVal};

  if (annot == null)
    return opts;

  const posOpen = annot.indexOf(BRACE_OPEN);
  const posClose = annot.indexOf(BRACE_CLOSE);

  if ((posOpen < 0) || (posClose < 0) || (posOpen >= posClose))
    return opts;

  let pos: number;
  let key: string;
  let val;

  annot.slice(posOpen + 1, posClose).split(ANNOT_SEPAR).forEach((str) => {
    pos = str.indexOf(CONTROL_SEP);

    if (pos === -1)
      return;

    key = str.slice(0, pos).trim();
    val = str.slice(pos + 1).trim();

    switch (key) {
    case 'min':
    case 'max':
      opts[key] = strToVal(val);
      break;

    case 'caption':
      opts[key] = val;
      break;
    }
  });

  return opts;
}; // getInputOpts

const strToVal = (s: string) => {
  const num = Number(s);
  return !isNaN(num) ? num : undefined;
};

/** Return default inputs of a model */
export function getDefaultInputs(optsMap: Map<string, InputOpts>): ModelInputs {
  const inputs: ModelInputs = {};

  optsMap.forEach((opts, name) => inputs[name] = opts.defaultVal);

  return inputs;
}

/** Return mapping captions to model input names */
export function getCaptionToInputNameMap(optsMap: Map<string, InputOpts>): Map<string, string> {
  const res = new Map<string, string>();

  optsMap.forEach((opts, name) => res.set(opts.caption ?? name, name));

  return res;
}
