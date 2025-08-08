## GA-production

Gluconic acid (GA) production by Aspergillus niger modeling

$$\begin{cases}

dX/dt = rX \\

dS/dt = -\gamma \cdot  rX - \lambda_1 \cdot  X \\

dO/dt = Kla \cdot  (Cod - O) - \delta2 \cdot  rX - \phi \cdot  X \\

dP/dt = \alpha \cdot  rX + \beta \cdot  X \\

X(0) = 5 \\

S(0) = 150 \\

O(0) = 7 \\

P(0) = 0 \\

\end{cases}$$

#### Auxiliary Computations

$$\mu = \mu_M \cdot  S / (Ks + S) \cdot  O / (Ko + O)$$

$$rX = \mu \cdot  X$$

#### Parameters

$$overall = 100$$

$$\mu_M = 0.668$$

$$\alpha = 2.92$$

$$\beta = 0.131$$

$$\gamma = 2.12$$

$$\lambda = 0.232$$

$$\delta = 0.278$$

$$\phi = 0.00487$$

$$Ks = 130.9$$

$$Ko = 0.000363$$

$$Kla = 0.017$$

$$Cod = 15$$

$$\frac{{x}^{2}}{y} - x
\int {x + x}\, dx + \mathrm{cos}\left(y+y\right)$$
