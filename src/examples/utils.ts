// Utils

import {ODEs, SolverMethod} from '../solver-tools';

/** Solve the given problem and output to console the solution */
export function printSolution(task: ODEs, method: SolverMethod, separator: string): void {
  try {
  // Solve the problem
    const solution = method(task);

    // Output results

    console.log(task.name);

    let line: string = [task.arg.name].concat(task.solutionColNames).join(separator);

    console.log(line);

    const length = solution[0].length;

    for (let i = 0; i < length; ++i) {
      line = solution.map((col) => col[i].toString()).join(separator);
      console.log(line);
    }
  } catch (err) {
    console.log('Solver failed: ', err instanceof Error ? err.message : 'Unknown problem!');
  }
}
