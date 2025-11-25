// Manager of callbacks

import {SolverOptions} from '../solver-defs';
import {Callback} from './callback-base';
import {TimeCheckerCallback} from './time-checker-callback';
import {IterCheckerCallback} from './iter-checker-callback';

/** Return callback corresponding to the given options
 * @param options - options of the numerical method
*/
export function getCallback(options?: Partial<SolverOptions>): Callback | undefined {
  if (options === undefined)
    return undefined;

  if (options.maxIterations !== undefined)
    return new IterCheckerCallback(options.maxIterations);

  if (options.maxTimeMs !== undefined)
    return new TimeCheckerCallback(options.maxTimeMs);

  return undefined;
}
