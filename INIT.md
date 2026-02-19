# Initial Step Size Selection for the MRT Method

**Context:** The MRT solver uses a single scalar `tolerance` parameter with the error control scheme:

$$\frac{|e_i|}{\mathrm{sc}_i} \leq \mathrm{tolerance}, \quad \mathrm{sc}_i = |y_i| + h|f_i| + \mathrm{TINY}$$

where `TINY` is a machine-small constant preventing division by zero.

This document describes the standard algorithm for selecting the initial step size $h_0$, adapted to this error control convention.

**Reference:** E. Hairer, S.P. Nørsett, G. Wanner, *Solving Ordinary Differential Equations I*, Springer, 1993, Section II.4.

---

## Algorithm

Given:

- Initial point $(t_0, y_0)$, system dimension $m$
- Right-hand side $f(t, y)$
- Scalar tolerance `tol` (≈ `rtol`)
- Constant `TINY` (≈ `atol`)
- Method order $p = 2$

---

### Step 1. Compute the scale vector

$$\mathrm{sc}_i = |y_{0,i}| + \mathrm{TINY}, \quad i = 1, \ldots, m$$

> **Note:** At the initial point, we do not yet have $h$, so the term $h|f_i|$ from the runtime
> error control is omitted. This is consistent with Hairer's approach, where the scale
> is based on the solution magnitude alone.

---

### Step 2. Evaluate the right-hand side

$$f_0 = f(t_0, y_0)$$

---

### Step 3. Compute norms of the initial data

$$d_0 = \sqrt{\frac{1}{m} \sum_{i=1}^{m} \left(\frac{y_{0,i}}{\mathrm{sc}_i}\right)^2}$$

$$d_1 = \sqrt{\frac{1}{m} \sum_{i=1}^{m} \left(\frac{f_{0,i}}{\mathrm{sc}_i}\right)^2}$$

> These are RMS norms of the initial state and the initial derivative, scaled by `sc`.
> The $L_\infty$ norm (consistent with the `max`-based error control in the code) can also be used:
>
> $$d_0 = \max_i \frac{|y_{0,i}|}{\mathrm{sc}_i}, \qquad d_1 = \max_i \frac{|f_{0,i}|}{\mathrm{sc}_i}$$

---

### Step 4. First guess for $h_0$

$$h_0 = \begin{cases} \dfrac{0.01 \cdot d_0}{d_1}, & \text{if } d_0 \geq 10^{-5} \text{ and } d_1 \geq 10^{-5}, \\[8pt] 10^{-6}, & \text{otherwise}. \end{cases}$$

The ratio $d_0 / d_1$ estimates the time scale on which $y$ changes by a fraction of itself. The factor $0.01$ ensures a conservative start.

---

### Step 5. Explicit Euler probe

Perform one explicit Euler step of size $h_0$ to estimate the second derivative:

$$y_1 = y_0 + h_0 \, f_0$$

$$f_1 = f(t_0 + h_0, \; y_1)$$

$$d_2 = \frac{1}{h_0}\sqrt{\frac{1}{m} \sum_{i=1}^{m} \left(\frac{f_{1,i} - f_{0,i}}{\mathrm{sc}_i}\right)^2}$$

> $d_2$ approximates $\|y''\| / \mathrm{sc}$ — the scaled norm of the second derivative.

---

### Step 6. Refine the guess using the method order

$$h_1 = \left(\frac{\mathrm{tol}}{\max(d_1, \, d_2)}\right)^{1/(p+1)}$$

For the MRT method with $p = 2$:

$$h_1 = \left(\frac{\mathrm{tol}}{\max(d_1, \, d_2)}\right)^{1/3}$$

> **Rationale:** The local truncation error of an order-$p$ method scales as $O(h^{p+1})$.
> Setting $h^{p+1} \cdot \max(d_1, d_2) \approx \mathrm{tol}$ and solving for $h$ gives this formula.

---

### Step 7. Final initial step size

$$h_0 \leftarrow \min\!\Big(100 \, h_0, \;\; h_1\Big)$$

The $\min$ ensures that the refinement from Step 6 does not overshoot the first
guess by more than a factor of 100.

---

## Summary

| Step | Computation | Cost |
|------|------------|------|
| 1 | Scale vector $\mathrm{sc}$ | $O(m)$ |
| 2 | $f_0 = f(t_0, y_0)$ | 1 evaluation of $f$ |
| 3 | Norms $d_0$, $d_1$ | $O(m)$ |
| 4 | First guess $h_0$ | $O(1)$ |
| 5 | Euler probe: $y_1$, $f_1$, $d_2$ | 1 evaluation of $f$ + $O(m)$ |
| 6 | Refined guess $h_1$ | $O(1)$ |
| 7 | Final $h_0$ | $O(1)$ |

**Total overhead:** 2 evaluations of $f$ and $O(m)$ work — negligible compared to the integration itself.

---

## Pseudocode

```
function initialStepSize(t0, y0, f, tol, TINY, p):
    m = length(y0)

    // Step 1: scale vector
    for i = 1..m:
        sc[i] = abs(y0[i]) + TINY

    // Step 2: initial right-hand side
    f0 = f(t0, y0)

    // Step 3: scaled norms
    d0 = rmsNorm(y0 / sc)
    d1 = rmsNorm(f0 / sc)

    // Step 4: first guess
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    // Step 5: Euler probe
    y1 = y0 + h0 * f0
    f1 = f(t0 + h0, y1)
    d2 = rmsNorm((f1 - f0) / sc) / h0

    // Step 6: refine using method order
    h1 = (tol / max(d1, d2)) ^ (1 / (p + 1))

    // Step 7: final value
    h0 = min(100 * h0, h1)

    return h0
```

---

## Adaptation Notes

1. **Norm choice.** The pseudocode uses the RMS norm. To be fully consistent with the
   `max`-based error control in the MRT code, replace `rmsNorm` with `maxNorm`
   ($L_\infty$). Both choices are valid; the RMS norm is slightly less conservative.

2. **The $h|f_i|$ term.** The runtime error control includes $h|f_i|$ in the scale. This term
   is intentionally omitted during initial step selection because $h$ is unknown at that
   point. Once the integration starts, the adaptive controller naturally adjusts $h$
   to be consistent with the full scale definition.

3. **Guard against zero.** If both $d_1$ and $d_2$ are extremely small (e.g., constant
   initial data with $f_0 \approx 0$), set $h_0$ to a small default like $10^{-6}$ or
   $h_0 = |t_1 - t_0| \cdot \mathrm{tol}^{1/(p+1)}$ to avoid division by zero.

4. **Step size bounds.** The result should be clipped to $[h_{\min}, h_{\max}]$ where
   $h_{\max}$ is typically $|t_1 - t_0|$ (the full integration interval).
