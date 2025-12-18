---
title: 'Diff Studio: Ecosystem for Interactive Modeling by Ordinary Differential Equations'
tags:
  - Differential equations
  - TypeScript
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

Ordinary differential equations (ODEs) are crucial in modeling complex systems and phenomena. Their applications range from pharmacology and drug manufacturing to financial modeling and environmental studies.

**Diff Studio** is a high-performance, interactive TypeScript application for solving initial value problems (IVPs) for ODEs directly within web browsers. It consists of two components. The **Diff Grok library** implements numerical methods and formula parsing tools. The **Diff Studio application** integrates Diff Grok tools with **Datagrok**, a scientific computing platform free for personal and academic use.

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

Current ODE modeling tools have significant drawbacks. Commercial systems are costly, blocking access for underfunded researchers, small labs, and institutions. Sharing simulations requires recipients to own the same software, complicating collaboration and reproducibility. Traditional platforms offer static interfaces. Development of interactive applications requires expertise, shifting focus from research to development.

This project proposes a web-based approach that addresses the stated problems. The goal is to develop an ecosystem providing a combination of a "no-code" approach with comprehensive capabilities for in-browser modeling and analysis.

# The solution: Diff Studio

The Diff Grok library provides numerical methods for solving problems given in a declarative form. It includes:

- **Solving tools:** numerical methods and computational pipelines.
- **Scripting tools:** features for specifying problems in declarative form.

Solving tools implement: the modified Rosenbrock triple (MRT) [@Shampine1997], the ROS3PRw [@jax2021], and the ROS34PRw [@rang2015improved] methods that
provide solving of both stiff and non-stiff systems. The performance is benchmarked on **Robertson** [@robertson1966solution], **HIRES** [@schafer1975new], **VDPOL** [@vanderpol1926relaxation], **OREGO** [@hairer2002solving2], **E5** [@hairer2002solving2], and **Pollution** [@verwer1994gauss]. Diff Grok allows users to obtain modeling results in near-real time (see \autoref{fig:peformance}).

![Diff Grok performance: computational time comparison.\label{fig:peformance}](./images/dg-performance.png)

Computational pipelines support multi-stage modeling and solving of IVPs in web workers, enabling parallel computations and analyses, including parameter optimization and sensitivity studies.

Scripting tools provide the ability to define IVPs as a set of text strings containing equations and model input annotations.
Their combination with solver tools makes Diff Grok a foundation for creating web applications for ODE-based modeling.

The Diff Studio package integrates Diff Grok with the Datagrok platform [@datagrok]. It implements a web application with a model editor (\autoref{fig:dseditor}) and autogenerated user interface (\autoref{fig:autoui}).

![Pharmacokinetic-pharmacodynamic simulation with Diff Studio: the equation editor, numerical solution, and its visualization.\label{fig:dseditor}](./images/ds-editor.png)

![Diff Studio: autogenerated UI.\label{fig:autoui}](./images/ds-ui.png)

Datagrok provides in-browser computations and visualizations. Other features include:

- Sensitivity analysis and parameter optimization;
- Storing and sharing computations via URI.

Thus, Diff Studio serves as a comprehensive modeling environment.

## Existing Solutions

WebAssembly [@wasm2025], Pyodide [@pyodide2025], and pure implementation with JavaScript or TypeScript are the most promising approaches for web-based simulation with ODEs.

WebAssembly ensures near-native performance. One can implement numerical methods in C/C++ or Rust, and
then compile them into WebAssembly, making it possible to use existing solvers.

Each time the user updates the equations, recompilation is required, which is impractical when designing complex models. A WebAssembly-based distribution of Python is provided by Pyodide. It enables the application of NumPy, SciPy, and other scientific libraries. However, it introduces large package sizes and performance overhead, as well as less seamless integration with browser APIs compared to native JavaScript implementations. Thus, the use of a pure JavaScript/TypeScript implementation of in-browser ODE solvers is the most promising.

There exists a set of libraries that provide ODE solvers: Math.js [@mathjs], odex-js [@odexjs], and others. Their distinctive feature is the requirement of proficiency in the corresponding programming language, which contrasts with a no-code approach to model description.

The Diff Grok library and Diff Studio package democratize scientific computing and make numerical modeling more accessible to a broader audience providing:

- Declarative model specification.
- Ability to solve both non-stiff and stiff ODEs.
- Parallel computing.
- Seamless integration of computations with visualization tools.
- Modeling directly in the web browser.

## Availability

Diff Grok is available on [GitHub](https://github.com/datagrok-ai/diff-grok).

The Diff Studio package is accessible on [GitHub](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio),
while its documentation can be found at [Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio).

Run Diff Studio online [here](https://public.datagrok.ai/apps/DiffStudio), or complete an interactive [tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations).

# Acknowledgements

The authors are grateful to the entire **Datagrok Inc** team and to the **JnJ ModelHub** project team for their contributions and feedback, which significantly improved the project.

# Conflicts of interest

The authors declare no conflict of interest.

# References
