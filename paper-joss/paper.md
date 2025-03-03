
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
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: G. Zakharov
    affiliation: 3
affiliations:
 - name: Datagrok Inc
   index: 1
   ror: 00hx57361
 - name: JnJ
   index: 2
 - name: Wellcome Sanger Institute, UK
   index: 3
date: 13 August 2017
bibliography: paper.bib

# Summary
Differential equations play a crucial role in modeling complex systems and phenomena, 
with applications ranging from pharmacology and drug manufacturing 
to financial modeling and environmental studies.
Diff Studio is a TypeScript-based application designed for solving initial value problems 
in ordinary differential equations (ODEs) directly within web browsers, 
featuring real-time solution visualization. 
By allowing users to define models declaratively without the need for programming knowledge, 
it streamlines both model generation and exploration.
Diff Studio consists of two components. 
The first component is the [Diff Grok](https://github.com/datagrok-ai/DiffGrok) 
that implements numerical methods and formula parsing tools. 
The second component is [Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
an application that integrates these tools with [Datagrok](https://datagrok.ai/), 
a next-generation scientific computing platform. 

Diff Studio creates an ecosystem within the Datagrok platform 
that supports the quick development of specialized scientific applications 
while offering premises for model design and exploration. 
It enables the democratization of research by acting as a central repository for ODE models 
and guaranteeing the reproducibility and accessibility of results.

# Statement of need
Scientific modeling of complex processes and phenomena relies heavily on ordinary differential equations. 
Their applications span diverse fields, including physical processes [Carmen Chicone], 
biochemical kinetics [Ingalls Brian], drug delivery systems, [Mircioiu, C.], 
cloud computing [Einollah Jafarnejad Ghomi], and population dynamics [Alan Hastings].

While analytic methods provide exact solutions, they can be applied only to a limited class of ODEs. 
Moreover, in many cases, the use of analytic solutions proves impractical due to their complexity [Ernst Hairer]. 
Consequently, numerical methods, which compute approximate solutions, are often preferred.

Numerous methods for solving ODEs have been recently developed [Ernst Hairer 1, 2]. 
These methods have been implemented in various software tools, 
including libraries and packages for programming languages and scientific computing environments. 
Notable examples include SUNDIALS [gardner2022sundials, hindmarsh2005sundials], 
Julia Differential Equations package [rackauckas2017differentialequations], 
SciPy [2020SciPy-NMeth], Maple [Maple], Mathematica [Mathematica], 
Matlab [MATLAB], and deSolve [Soetaert K]. 

Though most existing tools require significant technical expertise, 
shifting emphasis from applied research to software development. 
The goal of this project is to develop  an ecosystem providing comprehensive capabilities 
for in-browser ODE modeling and analysis.

# The solution: Diff Studio
In this paper, we introduce Diff Studio – a JavaScript-based application 
designed to provide capabilities for in-browser ODE modeling and analysis 
for user without significant experience in software development.

Diff Studio facilitates rapid development and debugging of models defined by ODEs. 
It enables the efficient design of complex process and phenomena simulators with user-friendly interfaces.

The key features of Diff Studio include:
* Browser-based numerical solving of initial value problems for systems of ODEs;
* Declarative problem specification implementing a "no-code" approach;
* Support for solving both stiff and non-stiff ODE systems;
* Automatic generation of user interfaces;
* Interactive visualization and model exploration tools;
* Sensitivity analysis and parameter optimization functionality;
* Development of standalone applications directly from ODE’s system.
* Capabilities for sharing models and computational results;


The project architecture comprises two main components:
[Diff Grok](https://github.com/datagrok-ai/DiffGrok) 
is a library providing numerical methods and features 
for automatically generating JavaScript code to solve problems specified through declarative syntax. 
[Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio) 
is a production-ready modeling environment 
integrated within the Datagrok [https://datagrok.ai/], scientific computing platform.

## DiffGrok library

The Diff Grok library (DG-lib) consists of two main components:
* Solving tools: numerical methods providing approximate solutions of Initial Value Problems (IVPs);
* Scripting tools: methods for automatic generation of JavaScript code that  solves problems specified in declarative form.

**Solving tools** incorporate the following numerical methods:
* `mrt` - Modified Rosenbrock triple (MRT) [Shampine, Lawrence F.];
* `ros3prw` - the ROS3PRw method [Joachim Rang];
* `ros34prw` - the ROS34PRw method [Joachim Rang].

To solve
$$dy /dt = f(t, y)
y(t_0) = y_0$$

on the interval `[t0, t1]` with step size `h`, one must define an ODEs object. 
This object specifies the independent variable (t) and its range (`[t0, t1]`), 
solution grid step size (`h`), initial conditions ($y_0$), right-hand side of the ODEs, 
tolerance, and names of dependent variables. 
Subsequently, a selected method (`mrt`, `ros3prw` or `ros34prw`) is applied to this object. 
The output consists of an array of Float64Array arrays 
containing the values of the independent variable and the corresponding approximate solutions.

For example, consider 
$$
dx/dt = x + y - t                                                (1)
dy/dt = xy + t
x(0) = 1, y(0) = -1
t∊[0, 2], h = 0.001
$$

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

The following code:
```javascript
const solution = mrt(task);
```

solves the given problem. Here, solution contains three items:
`solution[0]` - values of t, i.e. the range `0..2` with the step `0.001`;
`solution[1]` - values of $x(t)$ at the points of this range;
`solution[2]` - values of $y(t)$ at the same points.

The library demonstrates outstanding computational performance. 
The subsequent section presents a performance analysis based on solving a set of classical benchmark problems.

The scripting tools enable specification of IVPs in a declarative form 
known as the **Diff Studio model** (DS-model), 
which employs an intuitive block-structured syntax. 
For example, the problem defined in (1) can be expressed as follows:

```
#name: Example
#equations:
    dx/dt = x + y - t
    dy/dt = x * y + t
#argument: t
    initial = 0
    final =  2
    step = 0.001
#inits:
    x = 1
    y = -1
#tolerance: 1e-7
```



