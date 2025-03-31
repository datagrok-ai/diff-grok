
---
title: 'Diff Studio: Ecosystem for Interactive Modeling by Ordinary Differential Equations'
tags:
  - Differential equations
  - JavaScript
authors:
  - name: Viktor Makarichev   
    orcid: 0000-0000-0000-0000
    affiliation: "1" # (Multiple affiliations must be quoted)
    corresponding: true # (This is how to denote the corresponding author)
  - name: Larisa Bankurova
    affiliation: "1"
  - name: Gennadii Zakharov
    affiliation: "1,3"
  - name: Leonid Stolbov
    affiliation: "1"
  - name: Steven Mehrman
    affiliation: "2"
  - name: Dan Skatov
    affiliation: "2"
  - name: Jeffrey Cohen
    affiliation: "2"
  - name: Paul Sass
    affiliation: "2"
  - name: Davit Rizhinashvili 
    affiliation: "1"
  - name: Andrew Skalkin
    affiliation: "1"
affiliations:
 - name: Datagrok Inc, USA
   index: 1
   ror: 00hx57361
 - name: Johnson & Johnson Inc, USA
   index: 2
 - name: Wellcome Sanger Institute, UK
   index: 3
date: 03 April 2025
bibliography: paper.bib
---

# Summary

Differential equations play a crucial role in modeling complex systems
and phenomena, with applications ranging from pharmacology and drug
manufacturing to financial modeling and environmental studies.

Diff Studio is a TypeScript-based application designed for solving
initial value problems in ordinary differential equations (ODEs)
directly within web browsers. 
It enables users without programming experience 
define models declaratively using a simple domain-specific language 
and offers real-time visualization of solutions. 
The integrated editor includes example models, syntax
highlighting, and error detection. Diff Studio is highly performant,
capable of solving stiff equations significantly faster than many other
tools.

Diff Studio consists of two components. 
The first component is the
[Diff Grok](https://github.com/datagrok-ai/DiffGrok)
that implements numerical methods and formula parsing tools. 
The second component is 
[Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
an application that integrates these tools with 
[Datagrok](https://datagrok.ai/), 
a next-generation scientific computing platform.

Diff Studio creates an ecosystem within the Datagrok platform that
supports the quick development of specialized scientific applications
based on differential equation models. 
Also, it provides features for model design and exploration. 
It democratize research by acting as a central repository for ODE models 
and guaranteeing the reproducibility and accessibility of results.

# Statement of need

Scientific modeling of complex processes and phenomena relies heavily on
ordinary differential equations. 
Their applications span diverse fields,
including physical processes [@chicone2006ordinary], 
biochemical kinetics [@ingalls2013mathematical], 
drug delivery systems [@mircioiu2019mathematical], 
cloud computing [@jafarnejad2019applying], 
and population dynamics [@hastings2013population].

While analytic methods provide exact solutions, 
they can be applied only to a limited class of ODEs. 
Moreover, in many cases, the use of analytic
solutions proves impractical due to their complexity [@hairer2008solving1].
Consequently, numerical methods, which compute approximate solutions,
are often preferred.

Numerous methods for solving ODEs have been recently developed 
[@hairer2008solving1; @hairer2002solving2]. 
These methods have been implemented in various software
tools, including libraries and packages for programming languages and
scientific computing environments. 
Notable examples include 
SUNDIALS [@gardner2022sundials; @hindmarsh2005sundials], 
Julia Differential Equations package [@rackauckas2017differentialequations], 
SciPy [@2020SciPyNMeth], 
Maple [@maple2025], 
Mathematica [@Mathematica2024],
Matlab [@MATLAB], 
and deSolve [@soetaert2010solving].

Though most existing tools require significant technical expertise,
shifting emphasis from applied research to software development. 
The goal of this project is to develop an ecosystem 
providing comprehensive capabilities for in-browser ODE modeling and analysis.

# The solution: Diff Studio

In this paper, we introduce **Diff Studio** – 
a JavaScript-based application designed to provide capabilities 
for rapid in-browser ODE modeling and analysis for users 
without significant experience in software development.

Diff Studio facilitates rapid development and debugging of models
defined by ODEs. 
It enables the efficient design of complex processes
and phenomena simulators with user-friendly interfaces.

The key features of **Diff Studio** include:

- Browser-based numerical solving of initial value problems for systems of ODEs;
- Declarative problem specification implementing a "no-code" approach;
- Support for solving both stiff and non-stiff ODE systems;
- Automatic generation of user interfaces;
- Interactive visualization and model exploration tools;
- Enables sensitivity analysis and parameter optimization functionality;
- Development of standalone applications directly from ODE’s system.
- Capabilities for sharing models and computational results;

The project architecture comprises two main components:

- [Diff Grok](https://github.com/datagrok-ai/DiffGrok),
  is a library providing numerical methods and features for
  automatically generating JavaScript code to solve problems specified
  through declarative syntax.

- **Diff Studio**, is a production-ready modeling environment integrated
  within the [Datagrok](https://datagrok.ai/), scientific computing platform.

## DiffGrok library

The Diff Grok library (DG-lib) consists of two main components:

- **Solving tools:** numerical methods providing approximate solutions
  of Initial Value Problems (IVPs);
- **Scripting tools:** methods for automatic generation of JavaScript
  code that solves problems specified in declarative form.

**Solving tools** incorporate the following numerical methods:

- `mrt` - Modified Rosenbrock triple (MRT) [@Shampine1997]
- `ros3prw` - the ROS3PRw method [@rang2021rosenbrock]
- `ros34prw` - the ROS34PRw method [@rang2015improved]

To solve
\begin{equation}\label{eq:diffeq}
\begin{split}
dy/dt = f(t, y), \\
y(t_0) = y_0
\end{split}
\end{equation}

on the interval $[t0, t1]$ with step size $h$, 
define an ODEs object. 
This object specifies the independent variable ($t$) 
and its range ($[t0, t1]$), 
solution grid step size ($h$), 
initial conditions ($y0$),
right-hand side of the ODEs, 
tolerance, 
and names of dependent variables. 
Subsequently, a selected method (`mrt`, `ros3prw` or `ros34prw`) 
is applied to this object. 
The output consists of an array of `Float64Array`
arrays containing the values of the independent variable and the
corresponding approximate solutions.

For example, consider
\begin{equation}\label{eq:ivp}
\begin{split}
dx/dt = x + y - t, \\
dy/dt = xy + t \\
x(0) = 1, y(0) = -1 \\
t \in [0, 2], h = 0.001
\end{split}
\end{equation}

In this case, the ODEs object is defined as follows:

```javascript
const task: ODEs = {
    name: 'Example',
    arg: {
        name: 't',
        start: 0,
        finish: 2,
        step: 0.001,
    },
    initial: [1, -1],
    func: (t: number, y: Float64Array, output: Float64Array) => {
        output[0] = y[0] + y[1] - t;
        output[1] = y[0] * y[1] + t;
    },
    tolerance: 1e-7,
    solutionColNames: ['x', 'y'],
};
```

The following code solves the given problem:
```javascript
const solution = mrt(task);
```

The solution contains three items:

- `solution[0]` - values of `t`, i.e. the range `0..2` with the step `0.001`;
- `solution[1]` - values of `x(t)` at the points of this range;
- `solution[2]` - values of `y(t)` at the same points.

The library demonstrates outstanding computational performance. 
The subsequent section presents a performance analysis based on solving a
set of classical benchmark problems.

The scripting tools enable specification of IVPs in a declarative form
known as the Diff Studio model (DS-model), which employs an intuitive
block-structured syntax. 
For example, the problem defined in \autoref{eq:ivp} can be
expressed as follows:

```
#name: Example

#equations:
dx/dt = x + y - t
dy/dt = x \* y + t

#argument: t

initial = 0
final = 2
step = 0.001

#inits:
x = 1
y = -1

#tolerance: 1e-7
```

Diff Studio model corresponding to \autoref{eq:ivp}

The method getIVP parses strings of a DS-model and produces an IVP
object specifying a problem. If a model contains invalid expressions, an
error is raised. And the method getJScode generates JavaScript code,
involving an appropriate ODEs object, that can be applied for
**in-browser** solving differential equations. 
Besides, an intermediate
IVP object can be used for the automatic creation of user interfaces. 
We have implemented these features when developing the Diff Studio
application for the Datagrok platform.

The getIVP method parses DS-model strings and generates an IVP object
that specifies the problem. 
If a model contains invalid expressions, an error is raised. 
The getJScode method generates JavaScript code 
that incorporates an appropriate ODEs object, 
enabling **in-browser** solution of differential equations. 
Furthermore, the intermediate IVP
object can be utilized for the automatic generation of user interfaces.
These features have been implemented in the development of the Diff
Studio application for Datagrok.

### Performance

DG-lib ensures efficient integration of both stiff and non-stiff ODEs. 
The **Table 1** compares the computational performance of MRT, ROS3PRw, and
ROS34PRw in solving these problems.
To evaluate performance, we employ the following collection of classical
benchmark problems:

- **Robertson:** A stiff system describing the kinetics of an
  autocatalytic reaction [@robertson1966solution];
- **HIRES:** A stiff system explaining photomorphogenesis phytochrome
  responses through a chemical reaction involving eight reactants 
  [@schafer1975new];
- **VDPOL:** A system of ODEs describing nonlinear vacuum tube circuit
  behavior [@vanderpol1926relaxation];
- **OREGO:** A stiff system simulating the Belousov-Zhabotinskii
  reaction [@hairer2002solving2];
- **E5:** A stiff system of nonlinear ODEs representing a chemical
  pyrolysis model [@hairer2002solving2];
- **Pollution:** A stiff system of nonlinear equations describing
  chemical reactions in an air pollution model developed at the Dutch
  National Institute of Public Health and Environmental Protection
  [@verwer1994gauss].

**Table 1.** Computational time comparison: MRT vs. ROS3PRw vs.
ROS34PRw on AMD Ryzen 5 5600H 3.30 GHz CPU \label{table:performance}

|  Problem  | Equations |     Segment      | Points | Tolerance | MRT, ms | ROS3PRw, ms | ROS34PRw, ms |
|:---------:|:---------:|:----------------:|:------:|:---------:|:-------:|:-----------:|:------------:|
| Robertson |     3     |  `[0, 10E+11]`   |  40K   |   1E-7    |   103   |     446     |     285      |
|   HIRES   |     8     | `[0, 321.8122]` |  32K   |   1E-10   |   222   |     362     |     215      |
|   VDPOL   |     2     |   `[0, 2000]`   |  20K   |   1E-12   |   936   |    1576     |     760      |
|   OREGO   |     3     |   `[0, 360]`    |  36K   |   1E-8    |   381   |     483     |     199      |
|    E5     |     4     |  `[0, 10E+13]`  |  40K   |   1E-6    |   14    |     17      |      8       |
| Pollution |    20     |    `[0, 60]`    |  30K   |   1E-6    |   36    |     50      |      23      |

The efficient computation allows users to see the modeling results in
near-real time, thereby facilitating interactive model exploration and
improvement.

## Diff Studio application

The Diff Studio application (**DS-app**) provides production-ready
capabilities for solving initial value problems (IVPs) directly within a
web browser. It integrates DG-lib with the 
[Datagrok JavaScript API](https://datagrok.ai/js-api/).

DS-app has a model editor (\autoref{fig:diffstudio}), 
where users can define mathematical expressions, 
such as those shown in Fig. 2. 
The application then parses these expressions 
and automatically generates the following components:

- **JavaScript script (JS-script):** The JS-script invokes methods from
  DS-tools, and the Datagrok platform executes it to compute the
  numerical solution of the given IVP;
- **User interface (UI):** The UI provides interactive controls for
  users to modify model inputs (\autoref{fig:autoui}). Whenever a user updates the
  inputs, Datagrok reruns the JS-script to recompute the solution. The
  results are displayed using a grid and line charts.

![The Diff Studio application: 
the equation editor, numerical solution of the problem \autoref{eq:ivp} 
and its visualization.\label{fig:diffstudio}](./images/diffstudio.png)

![The Diff Studio application, Autogenerated UI: Diffstudio
creates input entries for all variables listed in the equation editor.
Each time model inputs are changed, a solution is computed and
displayed.\label{fig:autoui}](./images/diffstudio_autogenerated_ui.png)

DS-app generates the user interface (UI) based on the model
specification. To enhance usability, each model input can be annotated
with options that define the caption, category, measurement units,
minimum and maximum values, and tooltips (see eq:ivp). Diff Studio
provides the UI shown in \autoref{fig:autoui}, 
where inputs are displayed with their
specified captions and units. 
Additionally, inputs belonging to the same
category are grouped together. 
When a user hovers over an element, 
a tooltip appears, providing additional context.

Applying input annotations results in a self-explanatory UI. 
If minimum and maximum values are specified, 
a slider is generated, allowing users to adjust input values efficiently. 
Combined with the high computational performance of DS-app, 
this feature enables real-time, 
interactive model exploration.

![The correspondence of input annotation from \autoref{fig:autoui} and UI
elements from \autoref{fig:diffstudio}
\label{fig:annotations2ui}](./images/annotations-to-ui.png)

DS-app highlights errors and invalid expressions \autoref{fig:errorhighlight}, 
which is particularly important when creating and debugging complex models.

![Error indication in the Diff Studio application.
\label{fig:errorhighlight}](./images/error_highlighting.png)

DS-models can be stored on the Datagrok platform or downloaded to local
storage as `.ivp` text files, where `ivp` denotes an initial value
problem. Additionally, computation results and their visualizations can
be saved as CSV and PNG files, respectively.

Thus, DS-app serves as a comprehensive scientific computing environment,
facilitating both the development and application of models defined by
systems of ordinary differential equations (ODEs).

## Integration with Datagrok

An autogenerated JS-script is a 
[Datagrok function](https://datagrok.ai/help/datagrok/concepts/functions/),
which constitutes a key feature for several reasons.

First, the platform executes the JS-script to compute a numerical
solution for a given problem. The output of these computations is a
dataframe, which is displayed in a grid and visualized using line charts. 
Users can further customize the output using 
[Datagrok viewers](https://datagrok.ai/help/visualize/viewers/).
This also enables seamless sharing of the unique URI and history of
results.

Second, the Datagrok platform offers advanced function analysis
capabilities, which can be applied to any numerical 
[platform function](https://datagrok.ai/help/compute/function-analysis).
Using 
[sensitivity analysis](https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis)
you can explore 
how changes in input parameters affect your model outputs.
[Parameter optimization](https://datagrok.ai/help/compute/function-analysis#parameter-optimization)
solves the reversed task, allowing you to
find input parameters that generate a specified model output. 
These tools can be accessed directly from Diff Studio, enabling comprehensive
model analysis and parameter fitting.

Third, the JS-script is modifiable. 
It can be extended to incorporate non-elementary and special functions, 
allowing for the simulation of a broader range of processes. 
Additionally, other Datagrok functions can be called within the script, 
further enhancing its functionality.

Finally, JS-script can be saved in Datagrok scripts. 
Hence, it can be called by other functions. 
Such a capability makes it possible to create
building blocks for more complex applications.

Furthermore, **Diff Studio** supports the handling of IVP files within Datagrok. 
Any file with the `.ivp` extension is automatically opened in
Diff Studio, triggering the computation process. 
Users can execute models simply by clicking on IVP files stored on the platform.
Additionally, they can share links to these files with other users,
enabling seamless collaboration. 
This feature allows researchers to efficiently share their results with colleagues, 
positioning Diff Studio as a centralized hub 
for models defined by ordinary differential equations (ODEs).

In summary, the comprehensive capabilities of Diff Studio establish an
interactive modeling ecosystem for ODE-based simulations.

## Availability

The Diff Grok library is available on
[GitHub](https://github.com/datagrok-ai/DiffGrok) 
and includes both solving and scripting tools, along with appropriate examples.

The Diff Studio application code is also accessible on
[GitHub](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
while its documentation can be found at
[Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio).
Sample models are also available at
[Datagrok help pages](https://datagrok.ai/help/compute/models).

To run the application, open the
[DiffStudio application](https://public.datagrok.ai/apps/DiffStudio).
First-time users can log in using a Google account. 
Also, an interactive
[tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations)
provides a step-by-step guide to using Diff Studio.

Moving forward, our efforts will focus on enhancing Diff Studio,
including the integration of GPU-accelerated computations.

# Acknowledgements 

# Conflicts of interest
Authors declare not conflicts of interest.

# References



