# Modified Rosenbrock Triple (MRT)

**Source:** L.F. Shampine, M.W. Reichelt, "The MATLAB ODE Suite", SIAM J. Sci. Comput., Vol. 18, pp. 1–22, 1997, Section 3.1.

This method is the basis of the **`ode23s`** solver in MATLAB. It is a linearly implicit single-step method of order 2 for stiff ODEs, equipped with an embedded order-3 error estimate and a continuous extension (dense output).

---

## 1. Problem

$$y' = F(t, y), \quad y(t_0) = y_0$$

---

## 2. Method Constants

$$d = \frac{1}{2 + \sqrt{2}} \approx 0.2928932$$

$$e_{32} = 6 + \sqrt{2} \approx 7.4142136$$

---

## 3. Setup at Each Step

When stepping from $(t_n, y_n)$ to $t_{n+1} = t_n + h$:

- **Jacobian:** $J \approx \dfrac{\partial F}{\partial y}(t_n, y_n)$

- **Time derivative:** $T \approx \dfrac{\partial F}{\partial t}(t_n, y_n)$

- **Method matrix:** $W = I - h \, d \, J$

A **single LU factorization** of $W$ is performed.

---

## 4. Step Formulas

### Stage 1

$$F_0 = F(t_n, \; y_n)$$

$$W \, k_1 = F_0 + h \, d \, T$$

### Stage 2

$$F_1 = F\!\left(t_n + \frac{h}{2}, \;\; y_n + \frac{h}{2} \, k_1\right)$$

$$W \, k_2' = F_1 - k_1$$

$$k_2 = k_2' + k_1$$

### Solution Advance (order 2)

$$y_{n+1} = y_n + h \, k_2$$

### Stage 3 (for error estimation)

$$F_2 = F(t_{n+1}, \; y_{n+1})$$

$$W \, k_3 = F_2 - e_{32}(k_2 - F_1) - 2(k_1 - F_0) + h \, d \, T$$

### Local Error Estimate

$$\mathrm{err} = \frac{h}{6} \left( k_1 - 2 k_2 + k_3 \right)$$

---

## 5. FSAL (First Same As Last)

On a successful step, $F_2$ from the current step is reused as $F_0$ for the next step:

$$F_0^{\,\text{next}} = F_2$$

This saves one evaluation of $F$ per successful step.

---

## 6. Adaptive Step Size Control

Weighted error norm:

$$e = \sqrt{ \frac{1}{m} \sum_{i=1}^{m} \left( \frac{\mathrm{err}_i}{\mathrm{atol} + \mathrm{rtol} \cdot |y_{n+1,i}|} \right)^2 }$$

- **Step accepted** ($e \leq 1$): $\; h_{\text{new}} = h \cdot \min\!\left(\alpha_{\max}, \; \mathrm{safety} \cdot e^{-1/3}\right)$
- **Step rejected** ($e > 1$): $\; h_{\text{new}} = h \cdot \max\!\left(\alpha_{\min}, \; \mathrm{safety} \cdot e^{-1/3}\right)$

Typical parameters: $\mathrm{safety} \approx 0.9$, $\alpha_{\max} \approx 5$, $\alpha_{\min} \approx 0.2$.

---

## 7. Continuous Extension (Dense Output)

Approximation of $y(t_n + s \, h)$ for any $s \in [0, 1]$:

$$\boxed{ y_{n+s} = y_n + h \left[ \frac{s(1 - s)}{1 - 2d} \; k_1 \;+\; \frac{s(s - 2d)}{1 - 2d} \; k_2 \right] }$$

This is a quadratic polynomial in $s$, using only the already computed $k_1$ and $k_2$.

### Shorthand Notation

$$b_1(s) = \frac{s(1 - s)}{1 - 2d}, \qquad b_2(s) = \frac{s(s - 2d)}{1 - 2d}$$

$$y_{n+s} = y_n + h \big[ b_1(s) \, k_1 + b_2(s) \, k_2 \big]$$

### Numerical Values

With $d = \frac{1}{2+\sqrt{2}}$:

$$1 - 2d = \frac{\sqrt{2}}{2 + \sqrt{2}} = \sqrt{2} \, d \approx 0.4142136$$

$$2d \approx 0.5857864$$

### Boundary Conditions Check

| $s$ | $b_1(s)$ | $b_2(s)$ | $y_{n+s}$ |
|-----|----------|----------|-----------|
| $0$ | $0$      | $0$      | $y_n$ ✓  |
| $1$ | $0$      | $1$      | $y_n + h \, k_2$ ✓ |

---

## 8. Cost per Successful Step

| Operation | Count |
|-----------|-------|
| Evaluations of $F$ | 2 (or 1 with FSAL) |
| Jacobian evaluations | 1 (can be frozen) |
| LU factorization of $W$ | 1 |
| Back-solves | 3 |

### Cost per Interpolation Query

| Operation | Count |
|-----------|-------|
| Evaluations of $F$ | 0 |
| Linear system solves | 0 |
| Vector operations $O(m)$ | 2 multiplications + 2 additions |

---

## 9. Method Properties

- **Order:** 2 (solution advance), 3 (error estimate)
- **Stability:** L-stable when $J = \partial F / \partial y$
- **Type:** linearly implicit (Modified Rosenbrock / W-method)
- **Dense output:** quadratic interpolant, consistent with the method order
- **Best suited for:** stiff systems at moderate tolerances ($10^{-3}$–$10^{-6}$)
