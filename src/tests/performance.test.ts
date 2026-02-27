// Performance tests

import {perfProbs} from '../../index';
import {implicitMethods, perfExclusions, TIMEOUT_MS} from './test-defs';

implicitMethods.forEach((method, name) => {
  describe(`${name}`, () => {
    const excluded = perfExclusions.get(name);
    perfProbs.forEach((odes) => {
      const testFn = excluded?.has(odes.name) ? test.skip : test;
      testFn(`${odes.name}`, () => {
        const start = performance.now();
        method(odes);
        const elapsed = performance.now() - start;
        expect(elapsed).toBeLessThan(TIMEOUT_MS);
      });
    });
  });
});
