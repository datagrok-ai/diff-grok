/** Kintetic constants for the Pollution model */
var POL;
(function (POL) {
    POL[POL["K1"] = 0.35] = "K1";
    POL[POL["K2"] = 26.6] = "K2";
    POL[POL["K3"] = 12300] = "K3";
    POL[POL["K4"] = 0.00086] = "K4";
    POL[POL["K5"] = 0.00082] = "K5";
    POL[POL["K6"] = 15000] = "K6";
    POL[POL["K7"] = 0.00013] = "K7";
    POL[POL["K8"] = 24000] = "K8";
    POL[POL["K9"] = 16500] = "K9";
    POL[POL["K10"] = 9000] = "K10";
    POL[POL["K11"] = 0.022] = "K11";
    POL[POL["K12"] = 12000] = "K12";
    POL[POL["K13"] = 1.88] = "K13";
    POL[POL["K14"] = 16300] = "K14";
    POL[POL["K15"] = 4800000] = "K15";
    POL[POL["K16"] = 0.00035] = "K16";
    POL[POL["K17"] = 0.0175] = "K17";
    POL[POL["K18"] = 100000000] = "K18";
    POL[POL["K19"] = 444000000000] = "K19";
    POL[POL["K20"] = 1240] = "K20";
    POL[POL["K21"] = 2.1] = "K21";
    POL[POL["K22"] = 5.78] = "K22";
    POL[POL["K23"] = 0.0474] = "K23";
    POL[POL["K24"] = 1780] = "K24";
    POL[POL["K25"] = 3.12] = "K25";
})(POL || (POL = {}));
; // POL
/** The chemical reaction part of the air pollution model (https://archimede.uniba.it/~testset/report/pollu.pdf) */
export const pollution = {
    name: 'Pollution',
    arg: { name: 't', start: 0, finish: 60, step: 0.002 },
    initial: [0, 0.2, 0, 0.04, 0, 0, 0.1, 0.3, 0.01, 0, 0, 0, 0, 0, 0, 0, 0.007, 0, 0, 0],
    func: (t, y, output) => {
        // extract function values
        const y1 = y[0];
        const y2 = y[1];
        const y3 = y[2];
        const y4 = y[3];
        const y5 = y[4];
        const y6 = y[5];
        const y7 = y[6];
        const y8 = y[7];
        const y9 = y[8];
        const y10 = y[9];
        const y11 = y[10];
        const y12 = y[11];
        const y13 = y[12];
        const y14 = y[13];
        const y15 = y[14];
        const y16 = y[15];
        const y17 = y[16];
        const y18 = y[17];
        const y19 = y[18];
        const y20 = y[19];
        // evaluate expressions
        const r1 = POL.K1 * y1;
        const r2 = POL.K2 * y2 * y4;
        const r3 = POL.K3 * y5 * y2;
        const r4 = POL.K4 * y7;
        const r5 = POL.K5 * y7;
        const r6 = POL.K6 * y7 * y6;
        const r7 = POL.K7 * y9;
        const r8 = POL.K8 * y9 * y6;
        const r9 = POL.K9 * y11 * y2;
        const r10 = POL.K10 * y11 * y1;
        const r11 = POL.K11 * y13;
        const r12 = POL.K12 * y10 * y2;
        const r13 = POL.K13 * y14;
        const r14 = POL.K14 * y1 * y6;
        const r15 = POL.K15 * y3;
        const r16 = POL.K16 * y4;
        const r17 = POL.K17 * y4;
        const r18 = POL.K18 * y16;
        const r19 = POL.K19 * y16;
        const r20 = POL.K20 * y17 * y6;
        const r21 = POL.K21 * y19;
        const r22 = POL.K22 * y19;
        const r23 = POL.K23 * y1 * y4;
        const r24 = POL.K24 * y19 * y1;
        const r25 = POL.K25 * y20;
        // compute output
        output[0] = -(r1 + r10 + r14 + r23 + r24) + (r2 + r3 + r9 + r11 + r12 + r22 + r25);
        output[1] = -r2 - r3 - r9 - r12 + r1 + r21;
        output[2] = -r15 + r1 + r17 + r19 + r22;
        output[3] = -r2 - r16 - r17 - r23 + r15;
        output[4] = -r3 + 2 * r4 + r6 + r7 + r13 + r20;
        output[5] = -r6 - r8 - r14 - r20 + r3 + 2 * r18;
        output[6] = -r4 - r5 - r6 + r13;
        output[7] = r4 + r5 + r6 + r7;
        output[8] = -r7 - r8;
        output[9] = -r12 + r7 + r9;
        output[10] = -r9 - r10 + r8 + r11;
        output[11] = r9;
        output[12] = -r11 + r10;
        output[13] = -r13 + r12;
        output[14] = r14;
        output[15] = -r18 - r19 + r16;
        output[16] = -r20;
        output[17] = r20;
        output[18] = -r21 - r22 - r24 + r23 + r25;
        output[19] = -r25 + r24;
    },
    tolerance: 1e-6,
    solutionColNames: ['y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10', 'y11',
        'y12', 'y13', 'y14', 'y15', 'y16', 'y17', 'y18', 'y19', 'y20'],
}; // pollution
export const pollutionReferencePoint = new Float64Array([
    0.05646255480022769,
    0.1342484130422339,
    0.4139734331099427E-8,
    0.5523140207484359E-2,
    0.2018977262302196E-6,
    0.1464541863493966E-6,
    0.7784249118997964E-1,
    0.3245075353396018,
    0.7494013383880406E-2,
    0.1622293157301561E-7,
    0.1135863833257075E-7,
    0.2230505975721359E-2,
    0.2087162882798630E-3,
    0.1396921016840158E-4,
    0.8964884856898295E-2,
    0.4352846369330103E-17,
    0.6899219696263405E-2,
    0.1007803037365946E-3,
    0.1772146513969984E-5,
    0.5682943292316392E-4,
]);
//# sourceMappingURL=pollution.js.map