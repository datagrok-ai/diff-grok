import {printSolution} from './utils';
import {robertson} from './robertson';
import {hires} from './hires';
import {vdpol} from './vdpol';
import {orego} from './orego';
import {e5} from './e5';
import {pollution} from './pollution';
import {mrt, ros34prw} from '../solver-tools';

const SEPARATOR = '  '; // use ',' to get csv

/** Print solution of the ROBER benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printRobertson();
 * ```
 *
 * ```python
 * # Python reference calculation for comparison
 * import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

METHOD = "Radau" # "BDF" or "LSODA"

# -------------------------
# 1. Define Robertson ODE
# -------------------------
def robertson(t, y):
    A, B, C = y
    dA = -0.04 * A + 1e4 * B * C
    dB = 0.04 * A - 1e4 * B * C - 3e7 * B**2
    dC = 3e7 * B**2
    return [dA, dB, dC]

# -------------------------
# 2. Initial conditions
# -------------------------
y0 = [1.0, 0.0, 0.0]
t_start = 0.0
t_end = 1e11
t_step = 2.5e6
t_span = (t_start, t_end)

# Output points
t_eval = np.arange(t_start, t_end + t_step, t_step)

# -------------------------
# 3. Solve ODE
# -------------------------
sol = solve_ivp(
    robertson,
    t_span,
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-7
)

print("Success:", sol.success)
print("Message:", sol.message)

# sol.t → time array, sol.y.T → concentrations (n_points, 3)
data = np.column_stack((sol.t, sol.y.T))

# Create DataFrame with headers
df = pd.DataFrame(data, columns=['time', 'A', 'B', 'C'])

# Save to CSV
df.to_csv(f'ROBER-using-{METHOD}.csv', index=False)
```
 */
export function printRobertson(): void {
  printSolution(robertson, mrt, SEPARATOR);
}

/** Print solution of the HIRES benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printHires();
 * ```
 *
 * ```python
 * # Python reference calculation for comparison
 * import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# -------------------------
# 0. Solver method
# -------------------------
METHOD = "Radau"  # stiff solver suitable for HIRES

# -------------------------
# 1. Define HIRES ODE
# -------------------------
def hires(t, y):
    y1, y2, y3, y4, y5, y6, y7, y8 = y
    dy = np.zeros(8)
    dy[0] = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy[1] = 1.71 * y1 - 8.75 * y2
    dy[2] = -10.03 * y3 + 0.43 * y4 + 0.035 * y5
    dy[3] = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    dy[4] = -1.745 * y5 + 0.43 * y6 + 0.43 * y7
    dy[5] = -280 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7
    dy[6] = 280 * y6 * y8 - 1.81 * y7
    dy[7] = -280 * y6 * y8 + 1.81 * y7
    return dy

# -------------------------
# 2. Initial conditions and time
# -------------------------
y0 = [1, 0, 0, 0, 0, 0, 0, 0.0057]

t_start = 0.0
t_end = 321.8122
t_step = 0.01
n_points = int(np.ceil((t_end - t_start) / t_step)) + 1
t_eval = np.linspace(t_start, t_end, n_points)

# -------------------------
# 3. Solve ODE
# -------------------------
sol = solve_ivp(
    hires,
    t_span,
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=1e-10,
    atol=1e-10
)

print("Success:", sol.success)
print("Message:", sol.message)

# -------------------------
# 4. Save solution to CSV
# -------------------------
# Combine time + concentrations
data = np.column_stack((sol.t, sol.y.T))

# Column names
col_names = ['t'] + ['y1','y2','y3','y4','y5','y6','y7','y8']

# Create DataFrame and save
df = pd.DataFrame(data, columns=col_names)
filename = f'HIRES-using-{METHOD}.csv'
df.to_csv(filename, index=False)

print(f"Solution saved to '{filename}'")

```
 */
export function printHires(): void {
  printSolution(hires, ros34prw, SEPARATOR);
}

/** Print solution of the VDPOL benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printVdpol();
 * ```
 *
 * ```python
 * # Python reference calculation for comparison
 * import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# -------------------------
# 0. Solver method
# -------------------------
METHOD = "Radau"  # stiff solver recommended for van der Pol with large mu

# -------------------------
# 1. Define van der Pol ODE (stiff)
# -------------------------
def vdpol(t, y):
    x1, x2 = y
    dx1 = x2
    dx2 = -x1 + 1000 * (1 - x1**2) * x2  # mu = 1000
    return [dx1, dx2]

# -------------------------
# 2. Initial conditions and time
# -------------------------
y0 = [-1, 1]
t_start = 0.0
t_end = 2000.0
t_step = 0.1

# Use linspace to avoid "t_eval not in t_span" error
n_points = int(np.ceil((t_end - t_start) / t_step)) + 1
t_eval = np.linspace(t_start, t_end, n_points)

# -------------------------
# 3. Solve ODE
# -------------------------
sol = solve_ivp(
    vdpol,
    (t_start, t_end),
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=1e-12,
    atol=1e-12
)

print("Success:", sol.success)
print("Message:", sol.message)

# -------------------------
# 4. Save solution to CSV
# -------------------------
data = np.column_stack((sol.t, sol.y.T))
col_names = ['t', 'x1', 'x2']

df = pd.DataFrame(data, columns=col_names)
filename = f"VANDERPOL-using-{METHOD}.csv"
df.to_csv(filename, index=False)

print(f"Solution saved to '{filename}'")

```
 */
export function printVdpol(): void {
  printSolution(vdpol, ros34prw, SEPARATOR);
}

/** Print solution of the OREGO benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printOrego();
 * ```
 *
 * ```python
 * # Python reference calculation for comparison
 * import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# -------------------------
# 0. Solver method
# -------------------------
METHOD = "Radau"  # stiff solver recommended

# -------------------------
# 1. Define OREGO ODE
# -------------------------
def orego(t, y):
    y1, y2, y3 = y
    dy1 = 77.27 * (y2 - y1 * y2 + y1 - 0.000008375 * y1**2)
    dy2 = ( -y2 - y1 * y2 + y3 ) / 77.27
    dy3 = 0.161 * (y1 - y3)
    return [dy1, dy2, dy3]

# -------------------------
# 2. Initial conditions and time
# -------------------------
y0 = [1, 2, 3]
t_start = 0.0
t_end = 360.0
t_step = 0.01

# Use linspace to ensure t_eval within t_span
n_points = int(np.ceil((t_end - t_start) / t_step)) + 1
t_eval = np.linspace(t_start, t_end, n_points)

# -------------------------
# 3. Solve ODE
# -------------------------
sol = solve_ivp(
    orego,
    (t_start, t_end),
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=1e-8,
    atol=1e-8
)

print("Success:", sol.success)
print("Message:", sol.message)

# -------------------------
# 4. Save solution to CSV
# -------------------------
data = np.column_stack((sol.t, sol.y.T))
col_names = ['t', 'y1', 'y2', 'y3']

df = pd.DataFrame(data, columns=col_names)
filename = f"OREGO-using-{METHOD}.csv"
df.to_csv(filename, index=False)

print(f"Solution saved to '{filename}'")
```
 */
export function printOrego(): void {
  printSolution(orego, ros34prw, SEPARATOR);
}

/** Print solution of the E5 benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printE5();
 * ```
 *
 * ```python
 * # Python reference calculation for comparison
 * import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# -------------------------
# 0. Solver method
# -------------------------
METHOD = "Radau"  # stiff solver recommended

# -------------------------
# 1. Define E5 constants
# -------------------------
K1 = 7.89e-10
K2 = 1.13e9
K3 = 1.1e7
K4 = 1.13e3

# -------------------------
# 2. Define E5 ODE
# -------------------------
def e5(t, y):
    y1, y2, y3, y4 = y
    dy = np.zeros(4)
    dy[0] = -K1 * y1 - K3 * y1 * y3
    dy[1] =  K1 * y1 - K2 * y2 * y3
    dy[2] =  K1 * y1 - K2 * y2 * y3 - K3 * y1 * y3 + K4 * y4
    dy[3] =  K3 * y1 * y3 - K4 * y4
    return dy

# -------------------------
# 3. Initial conditions and time
# -------------------------
y0 = [0.00176, 0, 0, 0]
t_start = 0.0
t_end = 1e6
t_step = 2.5e3

# Ensure t_eval is within t_span
n_points = int(np.ceil((t_end - t_start) / t_step)) + 1
t_eval = np.linspace(t_start, t_end, n_points)

# -------------------------
# 4. Solve ODE
# -------------------------
sol = solve_ivp(
    e5,
    (t_start, t_end),
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=1e-6,
    atol=1e-6
)

print("Success:", sol.success)
print("Message:", sol.message)

# -------------------------
# 5. Save solution to CSV
# -------------------------
data = np.column_stack((sol.t, sol.y.T))
col_names = ['t', 'y1', 'y2', 'y3', 'y4']

df = pd.DataFrame(data, columns=col_names)
filename = f"E5-using-{METHOD}.csv"
df.to_csv(filename, index=False)

print(f"Solution saved to '{filename}'")
 */
export function printE5(): void {
  printSolution(e5, ros34prw, SEPARATOR);
}

/** Print solution of the POLL benchmark problem
 * @example
 * ```ts
 * // TypeScript usage example
 * printPollution();
 * ```
 *
 * ```python
 * # Python reference calculation for comparisonimport numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# -------------------------
# 0. Solver method and tolerance
# -------------------------
METHOD = "Radau"  # stiff solver
RTOL = 1e-6
ATOL = 1e-6

# -------------------------
# 1. Reaction constants
# -------------------------
k1 = 0.35
k2 = 26.6
k3 = 1.23e4
k4 = 8.6e-4
k5 = 8.2e-4
k6 = 1.5e4
k7 = 1.3e-4
k8 = 2.4e4
k9 = 1.65e4
k10 = 9e3
k11 = 0.022
k12 = 1.2e4
k13 = 1.88
k14 = 1.63e4
k15 = 4.8e6
k16 = 3.5e-4
k17 = 0.0175
k18 = 1e8
k19 = 4.44e11
k20 = 1240
k21 = 2.1
k22 = 5.78
k23 = 0.0474
k24 = 1780
k25 = 3.12

# -------------------------
# 2. Initial conditions
# -------------------------
y0 = [
    0,     # y1 NO2
    0.2,   # y2 NO
    0,     # y3 O3P
    0.04,  # y4 O3
    0,     # y5 HO2
    0,     # y6 OH
    0.1,   # y7 HCHO
    0.3,   # y8 CO
    0.01,  # y9 ALD
    0,     # y10 MEO2
    0,     # y11 C2O3
    0,     # y12 CO2
    0,     # y13 PAN
    0,     # y14 CH3O
    0,     # y15 HNO3
    0,     # y16 O1D
    0.007, # y17 SO2
    0,     # y18 SO4
    0,     # y19 NO3
    0      # y20 N2O5
]

# -------------------------
# 3. Define POLL ODE system
# -------------------------
def poll(t, y):
    y = np.array(y)
    r1  = k1*y[0]
    r2  = k2*y[1]*y[3]
    r3  = k3*y[4]*y[1]
    r4  = k4*y[6]
    r5  = k5*y[6]
    r6  = k6*y[6]*y[5]
    r7  = k7*y[8]
    r8  = k8*y[8]*y[5]
    r9  = k9*y[10]*y[1]
    r10 = k10*y[10]*y[0]
    r11 = k11*y[12]
    r12 = k12*y[9]*y[1]
    r13 = k13*y[13]
    r14 = k14*y[0]*y[5]
    r15 = k15*y[2]
    r16 = k16*y[3]
    r17 = k17*y[3]
    r18 = k18*y[15]
    r19 = k19*y[15]
    r20 = k20*y[16]*y[5]
    r21 = k21*y[18]
    r22 = k22*y[18]
    r23 = k23*y[0]*y[3]
    r24 = k24*y[18]*y[0]
    r25 = k25*y[19]

    dy = np.zeros(20)
    dy[0]  = -(r1 + r10 + r14 + r23 + r24) + (r2 + r3 + r9 + r11 + r12 + r22 + r25)
    dy[1]  = -r2 - r3 - r9 - r12 + r1 + r21
    dy[2]  = -r15 + r1 + r17 + r19 + r22
    dy[3]  = -r2 - r16 - r17 - r23 + r15
    dy[4]  = -r3 + 2*r4 + r6 + r7 + r13 + r20
    dy[5]  = -r6 - r8 - r14 - r20 + r3 + 2*r18
    dy[6]  = -r4 - r5 - r6 + r13
    dy[7]  = r4 + r5 + r6 + r7
    dy[8]  = -r7 - r8
    dy[9]  = -r12 + r7 + r9
    dy[10] = -r9 - r10 + r8 + r11
    dy[11] = r9
    dy[12] = -r11 + r10
    dy[13] = -r13 + r12
    dy[14] = r14
    dy[15] = -r18 - r19 + r16
    dy[16] = -r20
    dy[17] = r20
    dy[18] = -r21 - r22 - r24 + r23 + r25
    dy[19] = -r25 + r24

    return dy

# -------------------------
# 4. Time array
# -------------------------
t_start = 0
t_end = 60  # final time
t_step = 0.1
n_points = int(np.ceil((t_end - t_start)/t_step)) + 1
t_eval = np.linspace(t_start, t_end, n_points)

# -------------------------
# 5. Solve ODE
# -------------------------
sol = solve_ivp(
    poll,
    (t_start, t_end),
    y0,
    method=METHOD,
    t_eval=t_eval,
    rtol=RTOL,
    atol=ATOL
)

print("Success:", sol.success)
print("Message:", sol.message)

# -------------------------
# 6. Save solution to CSV
# -------------------------
col_names = ['t'] + [f'y{i+1}' for i in range(20)]
data = np.column_stack((sol.t, sol.y.T))
df = pd.DataFrame(data, columns=col_names)
filename = f"POLL-using-{METHOD}.csv"
df.to_csv(filename, index=False)

print(f"Solution saved to '{filename}'")
 */
export function printPollution(): void {
  printSolution(pollution, ros34prw, SEPARATOR);
}

printHires();
