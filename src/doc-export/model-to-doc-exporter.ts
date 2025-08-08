/* eslint-disable max-len */
import {IVP} from '../scripting-tools';
import {FormattingOptions, getDerForm, getOpenCloseEqTags, HEADER_TAG, TITLE} from './format-defs';
import {getFormatted} from './formatting-utils';


/** */
export class ModelToDocExporter {
  protected ivp: IVP;
  protected opts?: Partial<FormattingOptions>;

  constructor(ivp: IVP, opts?: Partial<FormattingOptions>) {
    this.ivp = ivp;
    this.opts = opts;
  };

  public getDocLines(): string[] {
    // Title
    const textLines: string[] = this.getDocBegin();

    // Differential equations
    textLines.push(...this.getEquationsSection());

    // Init conditions
    textLines.push(...this.getInitsSection());

    // Other components
    textLines.push(...this.getComponentsLines());

    return textLines;
  } // getDocLines

  public getShortDocLines(): string[] {
    // Title
    const textLines: string[] = this.getDocBegin();

    // Differential equations & init conditions
    textLines.push(...this.getIvpSection());

    // Other components
    textLines.push(...this.getComponentsLines());

    return textLines;
  } // getShortDocLines

  /** */
  private getComponentsLines(): string[] {
    const textLines: string[] = [];

    // Expressions
    if (this.ivp.exprs !== null)
      textLines.push(...this.getExpressionsSection());

    // Parameters
    if (this.ivp.params !== null)
      textLines.push(...this.getParamsSection());

    // Constants
    if (this.ivp.consts !== null)
      textLines.push(...this.getConstsSection());

    return textLines;
  }

  /** Returns initial value problem section */
  protected getIvpSection(): string[] {
    const textLines: string[] = [];
    const eqnForm = getDerForm(this.opts?.derForm);
    const tags = getOpenCloseEqTags(this.opts?.eqnTag);

    textLines.push(`${tags.open}\\begin{cases}`);

    // Equations
    const argName = this.ivp.arg.name;

    this.ivp.deqs.equations.forEach((equation, funcName) => {
      textLines.push(getFormatted(`${eqnForm(argName, funcName)} = ${equation} \\\\`));
    });

    // Init conditions
    const initArg = this.ivp.arg.initial.value;

    this.ivp.inits.forEach((input, name) => {
      textLines.push(getFormatted(`${name}(${initArg}) = ${input.value} \\\\`));
    });

    textLines.push(`\\end{cases}${tags.close}`);

    return textLines;
  }

  /** Return the model's begin */
  protected getSection(): string {
    return `${HEADER_TAG.SECTION} ${this.ivp.name}`;
  }

  /** Return the model's header: title & description */
  protected getDocBegin(): string[] {
    const textLines: string[] = [this.getSection()];

    // Description
    if (this.ivp.descr)
      textLines.push(`${this.ivp.descr}`);

    return textLines;
  }

  /** Return subsection begin */
  protected getSubsection(header: string): string {
    return `${HEADER_TAG.SUBSECTION} ${header}`;
  }

  /** Return the expressions section */
  protected getExpressionsSection(): string[] {
    const textLines: string[] = [];

    if (this.ivp.exprs === null)
      return textLines;

    const tags = getOpenCloseEqTags(this.opts?.eqnTag);

    // Section title
    textLines.push(this.getSubsection(TITLE.EXPR));

    this.ivp.exprs.forEach((right, left) => {
      textLines.push(getFormatted(`${tags.open}${left} = ${right}${tags.close}`));
    });

    return textLines;
  }

  /** Return the differential equations section */
  protected getEquationsSection(): string[] {
    const textLines: string[] = [];
    const eqnForm = getDerForm(this.opts?.derForm);
    const tags = getOpenCloseEqTags(this.opts?.eqnTag);

    // Equations
    const argName = this.ivp.arg.name;

    this.ivp.deqs.equations.forEach((equation, funcName) => {
      textLines.push(getFormatted(`${tags.open}${eqnForm(argName, funcName)} = ${equation}${tags.close}`));
    });

    return textLines;
  }

  /** Return the initial conditions section */
  protected getInitsSection(): string[] {
    const textLines: string[] = [];
    const tags = getOpenCloseEqTags(this.opts?.valTag);
    const initArg = this.ivp.arg.initial.value;

    // Section title
    textLines.push(this.getSubsection(TITLE.INITS));

    // Init conditions
    this.ivp.inits.forEach((input, name) => {
      textLines.push(getFormatted(`${tags.open}${name}(${initArg}) = ${input.value}${tags.close}`));
    });

    return textLines;
  }

  /** Return the parameters section */
  protected getParamsSection(): string[] {
    const textLines: string[] = [];

    if (this.ivp.params === null)
      return textLines;

    const tags = getOpenCloseEqTags(this.opts?.valTag);

    // Section title
    textLines.push(this.getSubsection(TITLE.PARAMS));

    this.ivp.params.forEach((inp, name) => {
      textLines.push(getFormatted(`${tags.open}${name} = ${inp.value}${tags.close}`));
    });

    return textLines;
  }

  /** Return the constants section */
  protected getConstsSection(): string[] {
    const textLines: string[] = [];

    if (this.ivp.consts === null)
      return textLines;

    const tags = getOpenCloseEqTags(this.opts?.valTag);

    // Section title
    textLines.push(this.getSubsection(TITLE.CONSTS));

    this.ivp.consts.forEach((inp, name) => {
      textLines.push(getFormatted(`${tags.open}${name} = ${inp.value}${tags.close}`));
    });

    return textLines;
  }
}; // ModelToDocExporter

/** */
export class ModelToLaTeXExporter extends ModelToDocExporter {
  constructor(ivp: IVP, opts?: Partial<FormattingOptions>) {
    super(ivp, opts);
  }

  protected getSection(): string {
    const ins = (this.ivp.descr !== null) ? '' : `\n\n${TITLE.DEQS}`;

    return `\\subsection*{${this.ivp.name}}${ins}`;
  }

  protected getSubsection(header: string): string {
    return `\\textbf{${header}}`;
  }
} // ModelToLaTeXExporter
