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

Ordinary differential equations (ODEs) are crucial in modeling complex systems and phenomena. Their applications range from pharmacology and drug manufacturing 
to financial modeling and environmental studies.

**Diff Studio** is a high-performance TypeScript application for solving initial value problems (IVPs) for ODEs directly within web browsers. It consists of two components. The **Diff Grok library** implements numerical methods and formula parsing tools. The **Diff Studio application** integrates Diff Grok tools with **Datagrok**, a scientific computing platform free for personal and academic use.

Diff Studio provides an ecosystem for rapid development of ODE-based applications with reproducible and accessible models.


# Statement of need

Scientific modeling of complex processes and phenomena uses ODEs. They are widely applied in diverse fields, including physical processes [@chicone2006ordinary], 
biochemical kinetics [@ingalls2013mathematical], drug delivery systems [@mircioiu2019mathematical], cloud computing [@jafarnejad2019applying], and population dynamics [@hastings2013population].

Analytic methods providing exact solutions can be applied only to a limited class of ODEs. 
The use of analytic solutions often proves impractical due to their complexity [@hairer2008solving1]. Numerical methods computing approximate solutions are often preferred.

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

Most tools require expertise, shifting focus from research to development. 
The goal of this project is to develop an ecosystem providing a combination of a "no-code" approach with comprehensive capabilities for in-browser modeling and analysis.

# The solution: Diff Studio

## Diff Grok library

This library provides numerical methods and automatic generation of JavaScript code from a declarative problem specification. It includes:

- **Solving tools:** numerical methods for solving IVPs;
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

on the interval $[t_0, t_1]$, 
define an ODEs object. 
This object specifies the independent variable ($t$), 
its range ($[_t0, t_1]$), 
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

Diff Grok delivers outstanding computational performance [@diffgrokperformance], benchmarked on **Robertson** [@robertson1966solution], **HIRES** [@schafer1975new], **VDPOL** [@vanderpol1926relaxation], **OREGO** [@hairer2002solving2], **E5** [@hairer2002solving2], and **Pollution** [@verwer1994gauss]. It allows users to obtain the modeling results in near-real time.

Scripting tools enable specification of IVPs declaratively using an intuitive syntax. For example, the problem \autoref{eq:ivp} can be expressed as shown in \autoref{fig:ivp}.

![Diff Studio model corresponding to the \autoref{eq:ivp}.\label{fig:ivp}](./images/DiffStudio_example_IVP.png){ height=6cm }

The method `getIVP()` parses a model and produces the IVP object specifying the problem. The method `getJScode()` generates JavaScript code, involving an appropriate ODEs object, that can be applied for solving equations.

## Diff Studio package

The Diff Studio package integrates Diff Grok with Datagrok [@datagrok]. It implements an application with a model editor (\autoref{fig:diffstudio}).

![The Diff Studio application: 
the equation editor, numerical solution of the problem \autoref{eq:ivp}, 
and its visualization.\label{fig:diffstudio}](./images/diffstudio.png)

Diff Studio automatically generates the user interface ( \autoref{fig:autoui}). Each model input can be annotated using self-explanatory options (\autoref{fig:annotations2ui}).

![The Diff Studio application, Autogenerated UI: Diffstudio
creates input entries (highlighted) for all variables listed in the equation editor.
Each time model inputs are changed, a solution is computed and
displayed.\label{fig:autoui}](./images/diffstudio_autogenerated_ui_highlighted.png)

![The correspondence of input annotation from \autoref{fig:diffstudio} and UI
elements from \autoref{fig:autoui}
\label{fig:annotations2ui}](./images/annotations-to-ui.png){ height=6cm }

Datagrok provides in-browser computations and visualizations. Other features include:

- Sensitivity analysis and parameter optimization functionality;
- Storing and sharing computations.

Thus, Diff Studio serves as a comprehensive modeling environment.

## Availability

Diff Grok is available on [GitHub](https://github.com/datagrok-ai/diff-grok).

The Diff Studio package is accessible on [GitHub](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
while its documentation can be found at [Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio).

Run Diff Studio online [here](https://public.datagrok.ai/apps/DiffStudio), or complete an interactive [tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations).

Moving forward, our efforts will focus on enhancing Diff Studio, including the integration of GPU-accelerated computations.

# Acknowledgements 

The authors are grateful to the entire **Datagrok Inc** team and to the **JnJ ModelHub** project team for their contributions and feedback, which significantly improved the project.

# Conflicts of interest

The authors declare no conflict of interest.

# References
