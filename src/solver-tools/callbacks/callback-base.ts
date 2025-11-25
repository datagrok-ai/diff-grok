// Solver's callback base

/** Solver callback */
export class Callback {
  constructor() {};

  /** Action performed on the numerical method iteration start */
  public onIterationStart(): void {}

  /** Action performed on computations end */
  public onComputationsCompleted(): void {}
};
