// Pipeline tests

import {problems, LOOP_OUTPUT_EXPR_MODEL, LOOP_OUTPUT_EXPR_INPUTS,
  UPDATE_OUTPUT_EXPR_MODEL, UPDATE_OUTPUT_EXPR_INPUTS, OUTPUT_EXPR_K} from './test-defs';
import {evalModel} from './test-utils';

problems.forEach((task, name) => {
  test(`Pipeline: model - ${name}`, () => {
    const solution = evalModel(task.model, task.inputs);
    expect(solution.length).toBe(task.outputsCount);
  });
});

// Expressions in #output with #loop / #update: verify the energy column (recomputed per stage)
// matches energy = k * exp(-t) * (x^2 + y^2) built from the returned t, x, y columns.
const outputExprCases = new Map<string, {model: string, inputs: Record<string, number>}>([
  ['loop', {model: LOOP_OUTPUT_EXPR_MODEL, inputs: LOOP_OUTPUT_EXPR_INPUTS}],
  ['update', {model: UPDATE_OUTPUT_EXPR_MODEL, inputs: UPDATE_OUTPUT_EXPR_INPUTS}],
]);

outputExprCases.forEach((task, name) => {
  test(`Pipeline: output expressions in ${name} model`, () => {
    const [t, x, y, energy] = evalModel(task.model, task.inputs);
    expect(energy.length).toBe(t.length);

    // Tolerance above the boundary argument-correction (DEFAULT_CORRECTION = 1e-7), which nudges
    // the returned t column after each stage computes energy from its own (uncorrected) argument.
    for (let i = 0; i < energy.length; ++i) {
      const expected = OUTPUT_EXPR_K * Math.exp(-t[i]) * (x[i] * x[i] + y[i] * y[i]);
      expect(Math.abs(energy[i] - expected)).toBeLessThan(1e-6);
    }
  });
});
