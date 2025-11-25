import {Callback} from './callback-base';
import {CallbackAction} from '../solver-defs';

/** This callback terminates computations if the maximum iterations is exceeded */
export class IterCheckerCallback extends Callback {
  private maxIter: number;
  private currentIter: number;

  /**
   * @param maxIter - the maximum number of iterations
   */
  constructor(maxIter: number) {
    super();
    this.maxIter = maxIter;
    this.currentIter = 0;
  }

  /** Action performed on iteration start of numerical method */
  public onIterationStart(): void {
    ++this.currentIter;

    if (this.currentIter > this.maxIter)
      throw new CallbackAction(`Max iterations count exceeded (${this.maxIter})`);
  }
}
