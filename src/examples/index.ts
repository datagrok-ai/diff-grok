import {robertson, robertsonReferencePoint} from './robertson';
import {hires, hiresReferencePoint} from './hires';
import {vdpol, vdpolReferencePoint} from './vdpol';
import {orego, oregoReferencePoint} from './orego';
import {e5, e5ReferencePoint} from './e5';
import {pollution, pollutionReferencePoint} from './pollution';

/** Problems for performance check
 * @internal
 */
export const perfProbs = [
  robertson,
  hires,
  vdpol,
  orego,
  e5,
  pollution,
];

/** Reference points for the performance problems
 * @internal
 */
export const refPoints = [
  robertsonReferencePoint,
  hiresReferencePoint,
  vdpolReferencePoint,
  oregoReferencePoint,
  e5ReferencePoint,
  pollutionReferencePoint,
];

export {corrProbs, CorrProblem} from './corr-probs';

export {printRobertson, printHires, printOrego, printE5, printVdpol, printPollution} from './print-benchmark';
