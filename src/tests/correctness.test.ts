// Correctness tests

import {corrProbs} from '../../index';
import {getError} from './test-utils';
import {methods, MAX_MAD} from './test-defs';

methods.forEach((method, name) => {
  describe(`${name}`, () => {
    corrProbs.forEach((problem) => {
      test(`${problem.odes.name}`, () => {
        const error = getError(method, problem);
        expect(error).toBeLessThan(MAX_MAD);
      });
    });
  });
});
