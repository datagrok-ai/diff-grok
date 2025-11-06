// Pipeline tests

import {problems} from './test-defs';
import {evalModel} from './test-utils';

problems.forEach((task, name) => {
  test(`Pipeline: model - ${name}`, () => {
    const solution = evalModel(task.model, task.inputs);
    expect(solution.length).toBe(task.outputsCount);
  });
});
