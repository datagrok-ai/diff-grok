// A script for checking performance

import {ODEs, corrProbs, CorrProblem, mrt, ros3prw, ros34prw, perfProbs, refPoints} from '../../index';

/** Return numerical solution error: maximum absolute deviation between approximate & exact solutions */
function getError(method: (odes: ODEs) => Float64Array[], corProb: CorrProblem): number {
  let error = 0;

  // Get numerical solution
  const approxSolution = method(corProb.odes);

  const exact = corProb.exact;

  const arg = approxSolution[0];

  const pointsCount = arg.length;
  const funcsCount = approxSolution.length - 1;

  // Compute error
  for (let i = 0; i < pointsCount; ++i) {
    const exactSolution = exact(arg[i]);

    for (let j = 0; j < funcsCount; ++j)
      error = Math.max(error, Math.abs(exactSolution[j] - approxSolution[j + 1][i]));
  }

  return error;
} // getError

function getDeviationFromReferencePoint(solution: Float64Array[], refPoint: Float64Array) {
  let absolute = 0;
  let relative = 0;
  const lastIdx = solution[0].length - 1;
  let cur = 0;

  for (let k = 0; k < solution.length; ++k) {
    //console.log(solution[k][lastIdx], refPoint[k]);
    cur = Math.abs(solution[k][lastIdx] - refPoint[k]);
    absolute = Math.max(cur, absolute);
    relative = Math.max(cur / (1 + Math.abs(refPoint[k])), relative);
  }

  return {
    absolute: absolute,
    relative: relative,
  };
} // getDeviationFromReferencePoint

const methods = new Map([
  ['MRT', mrt],
  ['ROS3PRw', ros3prw],
  ['ROS34PRw', ros34prw],
]);

console.log('Performance\n');

methods.forEach((method, name) => {
  console.log(' ', name);

  console.log('           PROBLEM   TIME,MS   ABS.ERR*    REL.ERR*');

  perfProbs.forEach((odes, idx) => {
    const start = Date.now();
    const solution = method(odes);
    const errs = getDeviationFromReferencePoint(solution.slice(1), refPoints[idx]);
    const absErr = errs.absolute.toExponential(2);
    const relErr = errs.relative.toExponential(2);
    const finish = Date.now();
    console.log(`  ${odes.name}     ${finish - start}     ${absErr}     ${relErr}`.padStart(50));
  });

  console.log();
});

console.log('* approximate vs. reference (also, numerical)\n\n');

console.log('Correctness (maximum absolute deviation): approximate vs. exact\n');

methods.forEach((method, name) => {
  console.log('  ', name);
  console.log('          PROBLEM   ABS. ERR');

  corrProbs.forEach((problem) => {
    const error = getError(method, problem);
    console.log(`     ${problem.odes.name}    ${error.toExponential(2)}`.padStart(28));
  });

  console.log();
});
