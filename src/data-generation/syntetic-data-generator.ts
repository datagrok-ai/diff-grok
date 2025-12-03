import {getIVP, IVP} from '../scripting-tools';
import {getIvp2WebWorker, IVP2WebWorker} from '../worker-tools';
import {InputOpts, ModelInputs, getCaptionToInputNameMap, getDefaultInputs, getInputOptionsMap} from './utils';

export class SyntheticDataGenerator {
  private ivp: IVP;
  private ivpWW: IVP2WebWorker;
  private inputOpts: Map<string, InputOpts>;
  private defaultInputs: ModelInputs;
  private captionToNameMap: Map<string, string>;

  constructor(model: string) {
    this.ivp = getIVP(model);
    this.ivpWW = getIvp2WebWorker(this.ivp);
    this.inputOpts = getInputOptionsMap(this.ivp);
    this.defaultInputs = getDefaultInputs(this.inputOpts);
    this.captionToNameMap = getCaptionToInputNameMap(this.inputOpts);

    console.log(this.ivp);
    console.log(this.inputOpts);
    console.log(this.defaultInputs);
    console.log(this.captionToNameMap);
  } // constructor

  public getDefaultInputs(): ModelInputs {
    return structuredClone(this.defaultInputs);
  }
} // SyntheticDataGenerator
