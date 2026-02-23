// Performance tests

import {perfProbs} from '../../index';
import {implicitMethods, TIMEOUT_MS} from './test-defs';

implicitMethods.forEach((method, name) => {
  describe(`${name}`, () => {
    perfProbs.forEach((odes) => {
      test(`${odes.name}`, () => {
        const start = performance.now();
        method(odes);
        const elapsed = performance.now() - start;
        expect(elapsed).toBeLessThan(TIMEOUT_MS);
      });
    });
  });
});
