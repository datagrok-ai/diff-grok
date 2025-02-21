// Features for creating output code for a pipeline

import {IVP} from '../scripting-tools';

export function getOutputCode(ivp: IVP): string | null {
    const outputs = ivp.outputs;
    if (outputs === null)
        return null;

    checkApplicability(ivp);

    const lines: string[] = ['const solution = [];'];

    if (exprsInOutputs(ivp)) {
        return null;
    }   

    outputs.forEach((_, name) => {
        if (name === ivp.arg.name)
            lines.push('solution.push(arguments[0][0]);');
        else {
            const idx = ivp.deqs.solutionNames.indexOf(name);

            if (idx > -1)
                lines.push(`solution.push(arguments[0][${idx + 1}]);`);
        }
    });

    lines.push('return solution;');

    return lines.join('\n');
} // getOutputCode

/** Check applicability of model outputs */
function checkApplicability(ivp: IVP): void {
    const outputs = ivp.outputs;

    if (outputs === null)
        throw new Error('Model has no outputs');

    const exprs = ivp.exprs;
    if (exprs === null)
        return;

    if (ivp.loop !== null) {
        if (exprsInOutputs(ivp))
            throw new Error('Non-supported model: expressions in output & loops');

        return;
    }

    if (ivp.updates !== null) {
        if (exprsInOutputs(ivp))
            throw new Error('Non-supported model: expressions in output and updates');
    }
} // checkApplicability

/** Check whether outputs have items from expressions */
function exprsInOutputs(ivp: IVP): boolean {
    const outputs = ivp.outputs;
    if (outputs === null)
        return false;

    const exprs = ivp.exprs;
    if (exprs === null)
        return false;

    let flag = false;

    exprs.forEach((_, name) => flag ||= outputs.has(name));

    return flag;
}
