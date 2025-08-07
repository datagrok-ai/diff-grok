## Advanced

#### Initial Value Problem

$$\begin{cases}

dx/dt = E1 * y + sin(t) \\

dy/dt = E2 * x - pow(t, 5) \\

x(0) = 2 \\

y(0) = 0 \\

\end{cases}$$

#### Auxiliary Computations

* $E1 = C1 * exp(-t) + P1$

* $E2 = C2 * cos(2 * t) + P2$

#### Parameters

* $P1 = 1$

* $P2 = -1$

#### Constants

$C1 = 1$

$C2 = 3$
