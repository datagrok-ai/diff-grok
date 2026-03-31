/* Function call to LaTeX transformer: maps math function calls to their LaTeX equivalents. */

/** Convert a function call to LaTeX.
 *  @param name  function name (e.g. 'sin', 'pow', 'Math.ceil')
 *  @param args  array of already-rendered LaTeX argument strings
 *  @returns     LaTeX string for the function call
 */
export function functionToLatex(name: string, args: string[]): string {
  switch (name) {
  case 'sin':
  case 'cos':
  case 'tan':
    return `\\${name}\\!\\left(${args[0]}\\right)`;

  case 'exp':
    return `e^{${args[0]}}`;

  case 'sqrt':
    return `\\sqrt{${args[0]}}`;

  case 'pow': {
    const base = isSimpleArg(args[0]) ? args[0] : `\\left(${args[0]}\\right)`;
    return `${base}^{${args[1]}}`;
  }

  case 'log':
    return `\\ln\\!\\left(${args[0]}\\right)`;

  case 'log10':
    return `\\log_{10}\\!\\left(${args[0]}\\right)`;

  case 'abs':
    return `\\left\\lvert ${args[0]} \\right\\rvert`;

  case 'Math.ceil':
  case 'ceil':
    return `\\left\\lceil ${args[0]} \\right\\rceil`;

  case 'Math.floor':
    return `\\left\\lfloor ${args[0]} \\right\\rfloor`;

  default:
    return `\\mathrm{${name}}\\!\\left(${args.join(', ')}\\right)`;
  }
}

/** Determine if an argument string is "simple" (single token, no binary ops).
 *  Used by pow() to decide whether the base needs parentheses.
 *  @param arg  rendered LaTeX argument string
 *  @returns    true if the argument has no top-level binary operators
 */
export function isSimpleArg(arg: string): boolean {
  // Simple if it doesn't contain binary operators at the top level
  // A simple arg: single identifier, number, or LaTeX command without + - * /
  // We check for unescaped +, -, spaces around operators
  return !/(?:^|[^\\])[+]/.test(arg) && !/\s-\s/.test(arg) && !/\s\*\s/.test(arg) && !/\s\/\s/.test(arg);
}
