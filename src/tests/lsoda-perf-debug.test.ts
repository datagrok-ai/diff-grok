// Quick debug test for LSODA Robertson
import {lsodaWeb, getCallback} from '../../index';
import {robertson} from '../examples/robertson';

test('LSODA Robertson', () => {
  // const cb = getCallback({maxIterations: 500000});
  // try {
  //   const start = performance.now();
  //   const result = lsodaWeb(robertson, cb);
  //   const elapsed = performance.now() - start;
  //   console.log(`Robertson: ${elapsed.toFixed(0)}ms`);
  //   expect(elapsed).toBeLessThan(10000);
  // } catch (e: any) {
  //   console.log(`Robertson FAILED: ${e.message}`);
  //   throw e;
  // }
});
