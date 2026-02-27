/** Kintetic constants for the E5 model */
var E5;
(function (E5) {
    E5[E5["K1"] = 7.89e-10] = "K1";
    E5[E5["K2"] = 1130000000] = "K2";
    E5[E5["K3"] = 11000000] = "K3";
    E5[E5["K4"] = 1130] = "K4";
})(E5 || (E5 = {}));
;
/** The E5 model (chemical pyrolysis: https://archimede.uniba.it/~testset/report/e5.pdf) */
export const e5 = {
    name: 'E5',
    arg: { name: 't', start: 0, finish: 1e13, step: 2.5e8 },
    initial: [0.00176, 0, 0, 0],
    func: (t, y, output) => {
        // extract function values
        const y1 = y[0];
        const y2 = y[1];
        const y3 = y[2];
        const y4 = y[3];
        // compute output
        output[0] = -E5.K1 * y1 - E5.K3 * y1 * y3;
        output[1] = E5.K1 * y1 - E5.K2 * y2 * y3;
        output[2] = E5.K1 * y1 - E5.K2 * y2 * y3 - E5.K3 * y1 * y3 + E5.K4 * y4;
        output[3] = E5.K3 * y1 * y3 - E5.K4 * y4;
    },
    tolerance: 1e-6,
    solutionColNames: ['y1', 'y2', 'y3', 'y4'],
};
export const e5ReferencePoint = new Float64Array([
    0.1152903278711829E-290,
    0.8867655517642120E-22,
    0.8854814626268838E-22,
    0,
]);
//# sourceMappingURL=e5.js.map