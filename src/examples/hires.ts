/** High Irradiance Responses of photomorphogenesis (see https://archimede.uniba.it/~testset/problems/hires.php) */
export const hires = {
  name: 'HIRES',
  arg: {name: 't', start: 0, finish: 321.8122, step: 0.01},
  initial: [1, 0, 0, 0, 0, 0, 0, 0.0057],
  func: (t: number, y: Float64Array, output: Float64Array) => {
    // extract function values
    const y1 = y[0];
    const y2 = y[1];
    const y3 = y[2];
    const y4 = y[3];
    const y5 = y[4];
    const y6 = y[5];
    const y7 = y[6];
    const y8 = y[7];

    // compute output
    output[0] = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007;
    output[1] = 1.71 * y1 - 8.75 * y2;
    output[2] = -10.03 * y3 + 0.43 * y4 + 0.035 * y5;
    output[3] = 8.32 * y2 + 1.71 * y3 - 1.12 * y4;
    output[4] = -1.745 * y5 + 0.43 * y6 + 0.43 * y7;
    output[5] = -280 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7;
    output[6] = 280 * y6 * y8 - 1.81 * y7;
    output[7] = -280 * y6 * y8 + 1.81 * y7;
  },
  tolerance: 1e-10,
  solutionColNames: ['y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8'],
}; // hires

export const hiresReferencePoint = new Float64Array([
  0.7371312573325668E-3,
  0.1442485726316185E-3,
  0.5888729740967575E-4,
  0.1175651343283149E-2,
  0.2386356198831331E-2,
  0.6238968252742796E-2,
  0.2849998395185769E-2,
  0.2850001604814231E-2,
]);
