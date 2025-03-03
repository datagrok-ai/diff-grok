import {IVP} from '../scripting-tools';
import {getOutputCode} from './output';
import {Pipeline} from './pipeline';
import {PipelineCreator} from './pipeline-creator';

/** Pipeline creator for model without loops & updates */
export class BasicModelPipelineCreator extends PipelineCreator {
  constructor(ivp: IVP) {
    super(ivp);
  }

  public getPipeline(inputs: Float64Array): Pipeline {
    return {
      wrappers: [
        {
          preproc: null,
          out: null,
          postproc: null,
        },
      ],
      out: getOutputCode(this.ivp),
    };
  }
} // BasicModelPipelineCreator
