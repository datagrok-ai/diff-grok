// Performance tests

import {perfProbs} from '../../index';
import {methods, TIMEOUT_MS} from './test-defs';

methods.forEach((method, name) => {
  perfProbs.forEach((odes) => {
    test(`Performance: method - ${name}, problem - ${odes.name}`, () => {
      const start = performance.now();
      method(odes);
      const elapsed = performance.now() - start;
      expect(elapsed).toBeLessThan(TIMEOUT_MS);
    });
  });
});
