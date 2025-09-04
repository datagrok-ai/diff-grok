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
date: 04 September 2025
bibliography: paper.bib
---

# Summary

Ordinary differential equations (ODEs) are very important in modeling complex systems and phenomena. Their applications range from pharmacology and drug manufacturing 
to financial modeling and environmental studies.

**Diff Studio** is a high-performance TypeScript application 
for solving initial value problems for ODEs directly within web browsers. It enables non-programmers to define models via a domain-specific language with real-time visualization, syntax support, error detection, and one-click generation of standalone ODE applications.

Diff Studio consists of two components. 
The **Diff Grok library** (<https://github.com/datagrok-ai/diff-grok>)
implements numerical methods and formula parsing tools. 
The **Diff Studio application** (<https://datagrok.ai/help/compute/diff-studio>)
integrates Diff Grok tools with **Datagrok** (<https://datagrok.ai>),  a scientific computing platform free for personal and academic use.

Diff Studio constitutes an ecosystem for rapid development of ODE-based applications that ensure reproducibility and accessibility through a shared model repository.


# Statement of need

Scientific modeling of complex processes and phenomena uses ordinary differential equations. They are widely applied in diverse fields, including physical processes [@chicone2006ordinary], 
biochemical kinetics [@ingalls2013mathematical], drug delivery systems [@mircioiu2019mathematical], cloud computing [@jafarnejad2019applying], and population dynamics [@hastings2013population].

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

We introduce **Diff Studio** - a solution that lets users define complex models using the special declarative notation and interactively explore model behavior in a web browser. It consists of two components: the **Diff Grok library** [@diffgrok] and the **Diff Studio package** [@diffstudio].

## Diff Grok library

This library provides numerical methods and automatic generation of JavaScript code from a declarative problem specification. It includes:

- **Solving tools:** numerical methods for solving Initial Value Problems (IVPs);
- **Scripting tools:** methods for automatic generation of JavaScript
  code that solves problems specified in the declarative form.

Solving tools implement:

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

The Diff Grok library demonstrates outstanding computational performance [@diffgrokperformance] that is evaluated using the following collection of classical
benchmark problems: **Robertson** [@robertson1966solution], **HIRES** [@schafer1975new], **VDPOL** [@vanderpol1926relaxation], **OREGO** [@hairer2002solving2], **E5** [@hairer2002solving2], and **Pollution** [@verwer1994gauss]. It allows users to obtain the modeling results in near-real time, which facilitates interactive model exploration and improvement.

Scripting tools enable specification of IVPs in a declarative form, which uses an intuitive block-structured syntax. For example, the problem defined in \autoref{eq:ivp} can be
expressed as shown in \autoref{fig:ivp}.

![Diff Studio model corresponding to the \autoref{eq:ivp}.\label{fig:ivp}](./images/DiffStudio_example_IVP.png){ height=6cm }

The method `getIVP()` parses a Diff Studio model and produces the IVP object specifying the problem. The method `getJScode()` generates JavaScript code, involving an appropriate ODEs object, that can be applied for **in-browser** solving differential equations. Besides, an intermediate IVP object can be used for the automatic creation of user interfaces.

## Diff Studio package

The Diff Studio package integrates Diff Grok tools with the Datagrok platform [@datagrok]. It implements an application with a model editor (\autoref{fig:diffstudio}), where users can edit equations utilizing a "no-code" approach. Errors and invalid expressions are highlighted (\autoref{fig:errorhighlight}). 

![The Diff Studio application: 
the equation editor, numerical solution of the problem \autoref{eq:ivp}, 
and its visualization.\label{fig:diffstudio}](./images/diffstudio.png)

![Error indication in the Diff Studio application.
\label{fig:errorhighlight}](./images/error_highlighting.png)

Diff Studio automatically generates the user interface ( \autoref{fig:autoui}). Each model input can be annotated using self-explanatory options (\autoref{fig:annotations2ui}).

![The Diff Studio application, Autogenerated UI: Diffstudio
creates input entries (highlighted) for all variables listed in the equation editor.
Each time model inputs are changed, a solution is computed and
displayed.\label{fig:autoui}](./images/diffstudio_autogenerated_ui_highlighted.png)

![The correspondence of input annotation from \autoref{fig:diffstudio} and UI
elements from \autoref{fig:autoui}
\label{fig:annotations2ui}](./images/annotations-to-ui.png){ height=6cm }

Datagrok performs in-browser computations and displays the results using a grid and line chart.

Other features of **Diff Studio** include:

- Sensitivity analysis and parameter optimization functionality;
- Storing and sharing models and computational results.

Thus, Diff Studio serves as a comprehensive scientific computing and modeling environment.

## Availability

The Diff Grok library is available on [GitHub](https://github.com/datagrok-ai/diff-grok) 
and includes both solving and scripting tools, along with appropriate examples.

The Diff Studio application code is also accessible on [GitHub](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
while its documentation can be found at [Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio). Sample models are also available at
[Datagrok help pages](https://datagrok.ai/help/compute/models).

To run the application, open the [Diff Studio application](https://public.datagrok.ai/apps/DiffStudio). First-time users can log in using a Google account. Also, an interactive [tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations) provides a step-by-step guide to using Diff Studio.

Moving forward, our efforts will focus on enhancing Diff Studio, including the integration of GPU-accelerated computations.

# Acknowledgements 

The authors express sincere gratitude to the 
all **Datagrok Inc** team for their invaluable support, 
and to the **JnJ ModelHub** project team for their contributions, 
and feedback, which significantly improved the quality 
and capabilities of **Diff Studio**.

# Conflicts of interest

The authors declare no conflict of interest.

# References
