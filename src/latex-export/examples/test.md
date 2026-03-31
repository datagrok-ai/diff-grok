## Acid Production

*Gluconic acid (GA) production by Aspergillus niger modeling*

### Equations

$$
\begin{aligned}
  \frac{dX}{dt} &= \mathrm{rX} \\
  \frac{dS}{dt} &= -\gamma \cdot \mathrm{rX} - \lambda \cdot X \\
  \frac{dO}{dt} &= \mathrm{Kla} \cdot \left(\mathrm{Cod} - O\right) - \delta \cdot \mathrm{rX} - \phi \cdot X \\
  \frac{dP}{dt} &= \alpha \cdot \mathrm{rX} + \beta \cdot X
\end{aligned}
$$

$$t \in \left[0,\, 100\right], \quad \Delta t = 0.1$$

### Expressions

$$
\begin{aligned}
  \mu &= \frac{\mathrm{muM} \cdot S}{\mathrm{Ks} + S} \cdot \frac{O}{\mathrm{Ko} + O} \\
  \mathrm{rX} &= \mu \cdot X
\end{aligned}
$$

### Stage Transitions

**2-nd stage** (duration: $\mathrm{overall} - 60$)

Before this stage:

$$
\begin{aligned}
  S \leftarrow S + 70
\end{aligned}
$$

### Initial Conditions (t = 0)

| Name | Value | Units |
| --- | --- | --- |
| $X$ | $5$ | kg/m³ |
| $S$ | $150$ | kg/m³ |
| $O$ | $7$ | kg/m³ |
| $P$ | $0$ | kg/m³ |

### Parameters

| Name | Value | Units |
| --- | --- | --- |
| $\mathrm{overall}$ | $100$ | h |
| $\mathrm{muM}$ | $0.668$ | 1/h |
| $\alpha$ | $2.92$ |  |
| $\beta$ | $0.131$ | 1/h |
| $\gamma$ | $2.12$ |  |
| $\lambda$ | $0.232$ | 1/h |
| $\delta$ | $0.278$ |  |
| $\phi$ | $4.87e-3$ | 1/h |
| $\mathrm{Ks}$ | $1.309e2$ | g/L |
| $\mathrm{Ko}$ | $3.63e-4$ | g/L |
| $\mathrm{Kla}$ | $1.7e-2$ | 1/s |
| $\mathrm{Cod}$ | $15$ | kg/m³ |
