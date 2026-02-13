/** Robertson chemical reaction, updated version (see https://archimede.uniba.it/~testset/problems/rober.php) */
export const robertson = {
    name: 'Robertson',
    arg: { name: 't', start: 0, finish: 10e11, step: 2.5e6 },
    initial: [1, 0, 0],
    func: (t, y, output) => {
        output[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
        output[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * Math.pow(y[1], 2);
        output[2] = 3e7 * Math.pow(y[1], 2);
    },
    tolerance: 1e-7,
    solutionColNames: ['A', 'B', 'C'],
};
export const robertsonReferencePoint = new Float64Array([
    0.2083340149701255E-7,
    0.8333360770334713E-13,
    0.9999999791665050,
]);
//# sourceMappingURL=robertson.js.map