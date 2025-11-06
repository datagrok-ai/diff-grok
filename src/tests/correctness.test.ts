// Correctness tests

import {corrProbs} from '../../index';
import {getError} from './test-utils';
import {methods, MAX_MAD} from './test-defs';

methods.forEach((method, name) => {
  corrProbs.forEach((problem) => {
    test(`Correctness: method - ${name}, problem - ${problem.odes.name}`, () => {
      const error = getError(method, problem);
      expect(error).toBeLessThan(MAX_MAD);
    });
  });
});
