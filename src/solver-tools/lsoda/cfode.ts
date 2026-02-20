import {LsodaContext} from './common';

/**
 * Set coefficients for the current method.
 * Computes elco and tesco arrays for all orders.
 * meth=1: Adams (orders 1-12), meth=2: BDF (orders 1-5).
 */
export function cfode(ctx: LsodaContext, meth: number): void {
  const c = ctx.common!;
  const pc = new Float64Array(13);

  if (meth === 1) {
    c.elco[1][1] = 1.;
    c.elco[1][2] = 1.;
    c.tesco[1][1] = 0.;
    c.tesco[1][2] = 2.;
    c.tesco[2][1] = 1.;
    c.tesco[12][3] = 0.;
    pc[1] = 1.;
    let rqfac = 1.;

    for (let nq = 2; nq <= 12; nq++) {
      const rq1fac = rqfac;
      rqfac = rqfac / nq;
      const nqm1 = nq - 1;
      const fnqm1 = nqm1;
      const nqp1 = nq + 1;

      // Form coefficients of p(x)*(x+nq-1)
      pc[nq] = 0.;
      for (let i = nq; i >= 2; i--)
        pc[i] = pc[i - 1] + fnqm1 * pc[i];
      pc[1] = fnqm1 * pc[1];

      // Compute integral, -1 to 0, of p(x) and x*p(x)
      let pint = pc[1];
      let xpin = pc[1] / 2.;
      let tsign = 1.;
      for (let i = 2; i <= nq; i++) {
        tsign = -tsign;
        pint += tsign * pc[i] / i;
        xpin += tsign * pc[i] / (i + 1);
      }

      // Store coefficients in elco and tesco
      c.elco[nq][1] = pint * rq1fac;
      c.elco[nq][2] = 1.;
      for (let i = 2; i <= nq; i++)
        c.elco[nq][i + 1] = rq1fac * pc[i] / i;
      const agamq = rqfac * xpin;
      const ragq = 1. / agamq;
      c.tesco[nq][2] = ragq;
      if (nq < 12)
        c.tesco[nqp1][1] = ragq * rqfac / nqp1;
      c.tesco[nqm1][3] = ragq;
    }
    return;
  }

  // meth = 2 (BDF)
  pc[1] = 1.;
  let rq1fac = 1.;

  for (let nq = 1; nq <= 5; nq++) {
    const fnq = nq;
    const nqp1 = nq + 1;

    // Form coefficients of p(x)*(x+nq)
    pc[nqp1] = 0.;
    for (let i = nq + 1; i >= 2; i--)
      pc[i] = pc[i - 1] + fnq * pc[i];
    pc[1] *= fnq;

    // Store coefficients in elco and tesco
    for (let i = 1; i <= nqp1; i++)
      c.elco[nq][i] = pc[i] / pc[2];
    c.elco[nq][2] = 1.;
    c.tesco[nq][1] = rq1fac;
    c.tesco[nq][2] = nqp1 / c.elco[nq][1];
    c.tesco[nq][3] = (nq + 2) / c.elco[nq][1];
    rq1fac /= fnq;
  }
}
