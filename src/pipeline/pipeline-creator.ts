import {IVP} from '../scripting-tools';
import {getOutputNames} from './output';
import {Pipeline} from './pipeline';

/** Creator of a pipeline for in-webworker computations */
export abstract class PipelineCreator {
  protected ivp: IVP;
  protected outputNames: string[];

  constructor(ivp: IVP) {
    this.ivp = ivp;
    this.outputNames = getOutputNames(ivp);
  };

  public abstract getPipeline(inputs: Float64Array): Pipeline;
}
