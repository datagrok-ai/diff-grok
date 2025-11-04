import {ODEs, CorrProblem} from '../../index';

/** Return numerical solution error: maximum absolute deviation between approximate & exact solutions */
export function getError(method: (odes: ODEs) => Float64Array[], corProb: CorrProblem): number {
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
