import {BASE_GREEK_LETTERS, MathMLOptions, GREEKS_TAGGING_DEFAULT, ASTERIX_REPLACING_DEFAULT,
  MATH_FUNCS, MATH_FUNC_REPLACING_DEFAULT} from './format-defs';


// Helper: Capitalize first letter
function capitalize(word: string): string {
  return word.charAt(0).toUpperCase() + word.slice(1);
}

const ALL_GREEK_LETTERS = BASE_GREEK_LETTERS.flatMap((letter) => [letter, capitalize(letter)]);

function addBackslashToGreekLetters(text: string): string {
  const pattern = new RegExp(`(?<!\\\\)\\b(${ALL_GREEK_LETTERS.join('|')})(?=\\b|\\d|_)`, 'g');

  return text.replace(pattern, '\\$1');
}

function replaceAsteriskWithCdot(text: string): string {
  // This regex matches a single '*' that is not part of a '**' pair
  const pattern = /(?<!\*)\*(?!\*)/g;
  return text.replace(pattern, '\\cdot ');
}

function addBackslashToMathFunctions(text: string): string {
  // Build a regex to match any function name not already preceded by a backslash
  const pattern = new RegExp(`(?<!\\\\)\\b(${MATH_FUNCS.join('|')})\\b`, 'g');

  // Add backslash before each matched function
  return text.replace(pattern, '\\$1');
}

export function getFormatted(text: string, opts?: Partial<MathMLOptions>): string {
  let formatted = text;

  // Tag greeks
  if (opts?.toTagGreeks ?? GREEKS_TAGGING_DEFAULT)
    formatted = addBackslashToGreekLetters(formatted);

  // Replace '*'
  if (opts?.toReplaceAsterixWithCdot ?? ASTERIX_REPLACING_DEFAULT)
    formatted = replaceAsteriskWithCdot(formatted);

  // Tag math funcs
  if (opts?.toTagMathFuncs ?? MATH_FUNC_REPLACING_DEFAULT)
    formatted = addBackslashToMathFunctions(formatted);

  return formatted;
}
