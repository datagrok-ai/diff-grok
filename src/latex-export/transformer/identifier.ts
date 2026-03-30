/* Identifier-to-LaTeX transformer: converts identifier names to their LaTeX representations. */

const SPECIAL: Record<string, string> = {
  'PI': '\\pi',
  'Inf': '\\infty',
  'inf': '\\infty',
};

const GREEK_LOWER: Array<[string, string]> = [
  ['varepsilon', '\\varepsilon'],
  ['varsigma', '\\varsigma'],
  ['vartheta', '\\vartheta'],
  ['epsilon', '\\epsilon'],
  ['omicron', 'o'],
  ['upsilon', '\\upsilon'],
  ['lambda', '\\lambda'],
  ['varphi', '\\varphi'],
  ['varrho', '\\varrho'],
  ['alpha', '\\alpha'],
  ['delta', '\\delta'],
  ['gamma', '\\gamma'],
  ['kappa', '\\kappa'],
  ['omega', '\\omega'],
  ['sigma', '\\sigma'],
  ['theta', '\\theta'],
  ['beta', '\\beta'],
  ['iota', '\\iota'],
  ['zeta', '\\zeta'],
  ['chi', '\\chi'],
  ['eta', '\\eta'],
  ['phi', '\\phi'],
  ['psi', '\\psi'],
  ['rho', '\\rho'],
  ['tau', '\\tau'],
  ['mu', '\\mu'],
  ['nu', '\\nu'],
  ['pi', '\\pi'],
  ['xi', '\\xi'],
];

const GREEK_UPPER: Array<[string, string]> = [
  ['Epsilon', '\\mathrm{E}'],
  ['Omicron', '\\mathrm{O}'],
  ['Upsilon', '\\Upsilon'],
  ['Lambda', '\\Lambda'],
  ['Kappa', '\\mathrm{K}'],
  ['Sigma', '\\Sigma'],
  ['Theta', '\\Theta'],
  ['Omega', '\\Omega'],
  ['Alpha', '\\mathrm{A}'],
  ['Delta', '\\Delta'],
  ['Gamma', '\\Gamma'],
  ['Beta', '\\mathrm{B}'],
  ['Zeta', '\\mathrm{Z}'],
  ['Iota', '\\mathrm{I}'],
  ['Eta', '\\mathrm{H}'],
  ['Chi', '\\mathrm{X}'],
  ['Mu', '\\mathrm{M}'],
  ['Nu', '\\mathrm{N}'],
  ['Phi', '\\Phi'],
  ['Psi', '\\Psi'],
  ['Rho', '\\mathrm{P}'],
  ['Tau', '\\mathrm{T}'],
  ['Pi', '\\Pi'],
  ['Xi', '\\Xi'],
];

// Combined, sorted by key length descending
const ALL_GREEK: Array<[string, string]> = [...GREEK_LOWER, ...GREEK_UPPER]
  .sort((a, b) => b[0].length - a[0].length);

// Quick lookup for exact matches
const GREEK_MAP = new Map<string, string>(ALL_GREEK);

/** Convert an identifier to its LaTeX representation.
 *  @param name  identifier name (e.g. "alpha1", "CO2", "x")
 *  @returns     LaTeX string
 */
export function identifierToLatex(name: string): string {
  // 1. Special constants
  if (name in SPECIAL)
    return SPECIAL[name];

  // 2. Greek exact match
  const greek = matchGreekExact(name);
  if (greek !== undefined)
    return greek;

  // 3. Greek prefix + digit suffix
  const greekPrefix = matchGreekPrefix(name);
  if (greekPrefix !== undefined)
    return `${greekPrefix.latex}_{${greekPrefix.suffix}}`;

  // 4. Single letter (no digits)
  if (name.length === 1)
    return name;

  // 5. Single letter + digits
  if (name.length > 1 && isLetter(name[0]) && isAllDigits(name, 1))
    return `${name[0]}_{${name.slice(1)}}`;

  // 6. Chemical formula (all uppercase letters + digits, at least one digit)
  const chem = matchChemicalFormula(name);
  if (chem !== undefined)
    return chem;

  // 7. Multi-letter → \mathrm{name}
  return `\\mathrm{${name}}`;
}

/** Check if an identifier is a known Greek letter (exact match).
 *  @param name  identifier name
 *  @returns     LaTeX command or undefined
 */
export function matchGreekExact(name: string): string | undefined {
  return GREEK_MAP.get(name);
}

/** Check if an identifier starts with a Greek letter name followed by digits only.
 *  @param name  identifier name (e.g. "alpha1")
 *  @returns     latex command and digit suffix, or undefined
 */
export function matchGreekPrefix(name: string): { latex: string; suffix: string } | undefined {
  for (const [key, latex] of ALL_GREEK) {
    if (name.length > key.length && name.startsWith(key)) {
      const suffix = name.slice(key.length);
      // Only match if suffix is all digits (not starting with a letter)
      if (isAllDigits(suffix, 0))
        return {latex, suffix};
      // Suffix starts with a letter → not Greek
      return undefined;
    }
  }
  return undefined;
}

/** Detect and format chemical formula patterns like CO2, N2O5, SO4.
 *  @param name  identifier name
 *  @returns     LaTeX string or undefined if not a chemical formula
 */
export function matchChemicalFormula(name: string): string | undefined {
  // Chemical formula: all letters are uppercase, contains at least one digit
  let hasDigit = false;
  let hasMultipleLetters = false;
  let letterCount = 0;

  for (let i = 0; i < name.length; i++) {
    const ch = name[i];
    if (isDigit(ch))
      hasDigit = true;
    else if (isUpperCase(ch))
      letterCount++;
    else {
      // lowercase letter or other char → not a chemical formula
      return undefined;
    }
  }

  hasMultipleLetters = letterCount >= 2;
  if (!hasDigit || !hasMultipleLetters)
    return undefined;

  // Parse into segments: alternating letter groups and digit groups
  const segments: Array<{ letters: string; digits: string }> = [];
  let i = 0;
  while (i < name.length) {
    let letters = '';
    while (i < name.length && isUpperCase(name[i]))
      letters += name[i++];
    let digits = '';
    while (i < name.length && isDigit(name[i]))
      digits += name[i++];
    if (letters || digits)
      segments.push({letters, digits});
  }

  // Count how many segments have digits
  const segmentsWithDigits = segments.filter((s) => s.digits).length;

  if (segmentsWithDigits <= 1) {
    // Simple case: trailing digits only (CO2, HO2, MEO2, SO4)
    const allLetters = segments.map((s) => s.letters).join('');
    const trailingDigits = segments[segments.length - 1].digits;
    if (trailingDigits)
      return `\\mathrm{${allLetters}}_{${trailingDigits}}`;
    return `\\mathrm{${allLetters}}`;
  }

  // Interleaved case: N2O5, C2O3
  let inner = '';
  for (const seg of segments) {
    inner += seg.letters;
    if (seg.digits)
      inner += `_{${seg.digits}}`;
  }
  return `\\mathrm{${inner}}`;
}

function isLetter(ch: string): boolean {
  return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z');
}

function isUpperCase(ch: string): boolean {
  return ch >= 'A' && ch <= 'Z';
}

function isDigit(ch: string): boolean {
  return ch >= '0' && ch <= '9';
}

function isAllDigits(str: string, from: number): boolean {
  if (from >= str.length) return false;
  for (let i = from; i < str.length; i++)
    if (!isDigit(str[i])) return false;

  return true;
}
