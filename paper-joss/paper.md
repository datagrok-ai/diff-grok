---
title: 'Diff Studio: Ecosystem for Interactive Modeling by Ordinary Differential Equations'
tags:
  - Differential equations
  - JavaScript
authors:
  - name: Viktor Makarichev   
    orcid: 0000-0003-1481-9132
    affiliation: "1" 
    corresponding: true 
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
 - name: Johnson & Johnson Inc, USA
   index: 2
 - name: Wellcome Sanger Institute, UK
   index: 3
date: 03 April 2025
bibliography: paper.bib
---

# Summary

Differential equations play a crucial role in modeling complex systems and phenomena, 
with applications ranging from pharmacology and drug manufacturing 
to financial modeling and environmental studies.

**Diff Studio** is a high-performance TypeScript application 
for solving initial value problems in ordinary differential equations (ODEs)
directly within web browsers. 
It enables users without programming expertise to declaratively define models 
using a simple domain-specific language 
and provides real-time solution visualization,
syntax highlighting, error detection, 
and one-click creation of standalone applications directly from an ODE's system. 

Diff Studio consists of two components. 
The **Diff Grok library** (<https://github.com/datagrok-ai/diff-grok>)
implements numerical methods and formula parsing tools. 
The **Diff Studio application** (<https://datagrok.ai/help/compute/diff-studio>)
integrates Diff Grok tools with 
**Datagrok** (<https://datagrok.ai>), 
a next-generation scientific computing platform.

Diff Studio creates an ecosystem
for rapid development of scientific applications based on ODE models. 
It democratizes research by providing a central repository for ODE models 
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

Analytic methods providing exact solutions
can be applied only to a limited class of ODEs. 
The use of analytic
solutions often proves impractical due to their complexity [@hairer2008solving1].
Numerical methods computing approximate solutions
are often preferred.

Many methods for solving ODEs have been recently developed 
[@hairer2008solving1; @hairer2002solving2]. 
These methods have been implemented in various software tools, 
including libraries and packages for programming languages and
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

We introduce **Diff Studio** – 
a TypeScript-based application 
that enables rapid in-browser ODE modeling and analysis, 
even for users without significant software development experience. 
Diff Studio facilitates the rapid development of scientific applications, 
providing a user-friendly interface
for efficient development of ODE-based simulators for complex processes and phenomena.

The key features of **Diff Studio** include:

- Numerical solving of initial value problems for systems of ODEs directly in a web-browser;
- Declarative problem specification implementing a "no-code" approach;
- Solving both stiff and non-stiff ODE systems;
- Automatic generation of user interfaces;
- Interactive visualization and model exploration tools;
- Sensitivity analysis and parameter optimization functionality;
- Automatic creation of standalone applications directly from an ODE's system;
- Sharing models and computational results;

The application has two components:

- The **Diff Grok library** [@diffgrok] provides numerical methods and 
  automatic generation of JavaScript code from a declarative problem specification.
- The **Diff Studio package** [@diffstudio] is a production-ready modeling environment integrated
  within **Datagrok** [@datagrok], a scientific computing platform.

## Diff Grok library

The Diff Grok library consists of two main parts:

- **Solving tools:** numerical methods providing approximate solutions
  of Initial Value Problems (IVPs);
- **Scripting tools:** methods for automatic generation of JavaScript
  code that solves problems specified in the declarative form.

### Solving tools

Solving tools incorporate the following numerical methods:

- `mrt` - Modified Rosenbrock triple (MRT) [@Shampine1997]
- `ros3prw` - the ROS3PRw method [@jax2021]
- `ros34prw` - the ROS34PRw method [@rang2015improved]

To solve
\begin{equation}\label{eq:diffeq}
\begin{split}
dy/dt = f(t, y), \\
y(t_0) = y_0
\end{split}
\end{equation}

on the interval $[t0, t1]$, 
define an ODEs object. 
This object specifies the independent variable ($t$), 
its range ($[t0, t1]$), 
solution grid step size ($h$), 
initial conditions ($y_0$),
right-hand side of the ODEs, 
tolerance, 
and names of dependent variables. 
Next, apply a selected method (`mrt`, `ros3prw` or `ros34prw`) 
to this object. 
The output consists of a list of `float64`
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

- `solution[0]` - values of `t`, i.e., the range `0..2` with the step `0.001`;
- `solution[1]` - values of `x(t)` at the points of this range;
- `solution[2]` - values of `y(t)` at the same points.

The Diff Grok library demonstrates outstanding computational performance. 
The later section presents a performance analysis based on solving a
set of classical benchmark problems.

### Scripting tools

Scripting tools enable specification of IVPs in a declarative form
known as the **Diff Studio model**, 
which uses an intuitive block-structured syntax. 
For example, the problem defined in \autoref{eq:ivp} can be
expressed as shown in \autoref{fig:ivp}.

![Diff Studio model corresponding to the \autoref{eq:ivp}.\label{fig:ivp}](./images/DiffStudio_example_IVP.png){ height=6cm }

The method `getIVP()` parses a Diff Studio model and produces the IVP
object specifying the problem. 
If a model contains invalid expressions, an error is raised. 
The method `getJScode()` generates JavaScript code,
involving an appropriate ODEs object, that can be applied for
**in-browser** solving differential equations. 
Besides, an intermediate
IVP object can be used for the automatic creation of user interfaces. 

### Performance

Diff Grok efficiently solves both stiff and non-stiff ODEs. 
**Table 1** compares the computational performance of `MRT`, `ROS3PRw`, and
`ROS34PRw` in solving these problems.
For performance evaluation, we used the following collection of classical
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

**Table 1.** Computational time comparison: `MRT` vs. `ROS3PRw` vs.
`ROS34PRw` on AMD Ryzen 5 5600H 3.30 GHz CPU \label{table:performance}

|  Problem  | Equations |    Segment    | Points | Tolerance | MRT, ms | ROS3PRw, ms | ROS34PRw, ms |
|:---------:|:---------:|:-------------:|:------:|:---------:|:-------:|:-----------:|:------------:|
| Robertson |     3     | `[0, 10E+11]` |  40K   |   1E-7    |   103   |     446     |     285      |
|   HIRES   |     8     | `[0, 321.81]` |  32K   |   1E-10   |   222   |     362     |     215      |
|   VDPOL   |     2     |  `[0, 2000]`  |  20K   |   1E-12   |   936   |    1576     |     760      |
|   OREGO   |     3     |  `[0, 360]`   |  36K   |   1E-8    |   381   |     483     |     199      |
|    E5     |     4     | `[0, 10E+13]` |  40K   |   1E-6    |   14    |     17      |      8       |
| Pollution |    20     |   `[0, 60]`   |  30K   |   1E-6    |   36    |     50      |      23      |

Diff Grok allows users to obtain the modeling results in near-real time, 
which facilitates interactive model exploration and improvement.

## Diff Studio package

The Diff Studio package provides production-ready
capabilities for solving initial value problems (IVPs) 
directly within a web browser. 
It integrates Diff Grok with the **Datagrok** platform [@datagrok] 
using [Datagrok JavaScript API](https://datagrok.ai/js-api/).

Diff Studio has a model editor (\autoref{fig:diffstudio}), 
where users can edit ODEs system and all other equations 
defining IVP. 

![The Diff Studio application: 
the equation editor, numerical solution of the problem \autoref{eq:ivp}, 
and its visualization.\label{fig:diffstudio}](./images/diffstudio.png)

The application parses these expressions 
and automatically generates the following components:

- **JavaScript script** that calls methods from Diff Studio tools.
  The Datagrok platform executes it to compute the
  numerical solution of the given IVP;
- **User interface (UI)**, providing interactive controls for 
  users to modify model inputs.  
  Datagrok can automatically recompute the solution when the user updates the inputs. 
  The results are displayed using a grid and line chart.

Diff Studio generates the UI based on the model specification. 
Each model input can use options that define the caption, category, measurement units,
minimum and maximum values, and tooltips. 
Inputs from the same category are grouped together. 
When a user hovers over an element, 
a tooltip appears, providing additional context.
The example of Diff Studio autogenerated UI with
captions and units for inputs
is shown on the \autoref{fig:autoui}. 

![The Diff Studio application, Autogenerated UI: Diffstudio
creates input entries (highlighted) for all variables listed in the equation editor.
Each time model inputs are changed, a solution is computed and
displayed.\label{fig:autoui}](./images/diffstudio_autogenerated_ui_highlighted.png)

Applying input annotations results in a self-explanatory UI. 
If minimum and maximum values are specified, 
a slider is generated, allowing users to adjust input values efficiently. 
Combined with the high computational performance of Diff Grok, 
this feature enables real-time, 
interactive model exploration. 
The correspondence of input annotation and UI elements
is shown on the \autoref{fig:annotations2ui}

![The correspondence of input annotation from \autoref{fig:diffstudio} and UI
elements from \autoref{fig:autoui}
\label{fig:annotations2ui}](./images/annotations-to-ui.png){ height=6cm }

Diff Studio highlights errors and invalid expressions (\autoref{fig:errorhighlight}), 
which is particularly important when creating and debugging complex models.

![Error indication in the Diff Studio application.
\label{fig:errorhighlight}](./images/error_highlighting.png)

Diff Studio models can be stored on the Datagrok platform or downloaded to local
storage as `.ivp` text files, where "`ivp`" denotes an Initial Value Problem. 
Computation results and their visualizations can be saved as CSV and PNG files.

Thus, Diff Studio serves as a comprehensive scientific computing and modeling environment,
facilitating both the development and application of models 
defined by systems of ordinary differential equations (ODEs).

## Integration with the Datagrok platform

An autogenerated JS script solving ODEs is a 
[Datagrok function](https://datagrok.ai/help/datagrok/concepts/functions/),
which constitutes a key feature for several reasons.

First, Diff Studio returns calculation results in a Datagrok Dataframe. 
Users can use all Datagrok's interactive capabilities to [transform](https://datagrok.ai/help/transform/) 
and [visualize](https://datagrok.ai/help/visualize/viewers/) calculation results. 

Second, Datagrok offers advanced 
[function analysis capabilities]((https://datagrok.ai/help/compute/function-analysis)), 
which can be applied to any numerical Datagrok function.
Using 
[sensitivity analysis](https://datagrok.ai/help/compute/function-analysis#sensitivity-analysis)
the user can explore 
how changes in input parameters affect model outputs.
[Parameter optimization](https://datagrok.ai/help/compute/function-analysis#parameter-optimization)
solves the reversed task, allowing the user to
find an input parameters combination that generates a specified model output. 
These tools can be accessed directly from Diff Studio, enabling comprehensive
model analysis and parameter fitting.

Third, the user can use all enterprise Datagrok capabilities:
seamless sharing of the model and history of computations via unique URI,
access control, and audit.

Fourth, the JS script is modifiable. 
The user can incorporate non-elementary and special functions,
call other Datagrok functions,
and customize default visualizations via 
[function annotations](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation), 
further enhancing script functionality.

Finally, JS script can be saved in Datagrok scripts and called by other Datagrok functions. 
This capability makes it possible to use the script as a
building block for more complex applications.

Furthermore, Diff Studio supports the handling of IVP files within Datagrok. 
Users can execute models either by clicking on IVP files stored on the platform 
or by dragging and dropping local `.ivp` files into the Datagrok window. 
Any file with the `.ivp` extension is automatically opened in Diff Studio, 
triggering the computation process.

Additionally, the user can share links to these files with other users,
enabling seamless collaboration. 
This feature allows researchers to efficiently share their results with colleagues, 
positioning Diff Studio as a centralized hub 
for models defined by ordinary differential equations (ODEs).

In summary, the comprehensive capabilities of Diff Studio establish an
interactive modeling ecosystem for ODE-based simulations.

## Availability

The Diff Grok library is available on
[GitHub](https://github.com/datagrok-ai/diff-grok) 
and includes both solving and scripting tools, along with appropriate examples.

The Diff Studio application code is also accessible on
[GitHub](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
while its documentation can be found at
[Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio).
Sample models are also available at
[Datagrok help pages](https://datagrok.ai/help/compute/models).

To run the application, open the
[Diff Studio application](https://public.datagrok.ai/apps/DiffStudio).
First-time users can log in using a Google account. 
Also, an interactive
[tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations)
provides a step-by-step guide to using Diff Studio.

Moving forward, our efforts will focus on enhancing Diff Studio,
including the integration of GPU-accelerated computations.

# Acknowledgements 
The authors express sincere gratitude to the 
all **Datagrok Inc** team for their invaluable support, 
and to the **JnJ ModelHub** project team for their contributions, 
and feedback, which significantly improved the quality 
and capabilities of **Diff Studio**.

# Conflicts of interest
The authors declare no conflict of interest.

# References



