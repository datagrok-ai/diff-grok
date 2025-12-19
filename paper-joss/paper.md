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

Diff Studio is an ecosystem for ODE-based modeling within Datagrok. It enables rapid development of applications with reproducible computations and supports interactive exploration of model behavior. The platform facilitates collaboration through sharing, reuse, and community-driven development.

# Statement of need

Scientific modeling of complex processes and phenomena often uses ODEs. They are widely applied in diverse fields, including physical processes [@chicone2006ordinary],
biochemical kinetics [@ingalls2013mathematical], drug delivery systems [@mircioiu2019mathematical], cloud computing [@jafarnejad2019applying], and population dynamics [@hastings2013population].

Analytic methods providing exact solutions can be applied only to a limited class of ODEs. The use of analytic solutions often proves impractical due to their complexity [@hairer2008solving1]. Numerical methods computing approximate solutions are often preferred. Many such methods have been developed [@hairer2008solving1; @hairer2002solving2] and implemented in a variety of software tools, including libraries and packages for programming languages and scientific computing environments. Notable examples include SUNDIALS [@gardner2022sundials; @hindmarsh2005sundials], the Julia Differential Equations package [@rackauckas2017differentialequations], SciPy [@2020SciPyNMeth], Maple [@maple2025], Mathematica [@Mathematica2024], Matlab [@MATLAB], and deSolve [@soetaert2010solving].

For web-based ODE simulation, the following approaches are promising: WebAssembly [@wasm2025], Pyodide [@pyodide2025], and pure JavaScript/TypeScript implementations. WebAssembly provides near-native performance by allowing numerical methods written in C/C++ or Rust to be compiled for the browser, but it requires recompilation whenever equations are updated, limiting its flexibility for iterative model design. Pyodide offers a WebAssembly-based Python distribution that enables the use of NumPy, SciPy, and other scientific libraries, but it incurs large package sizes, performance overhead, and less seamless integration with browser APIs. In contrast, pure JavaScript/TypeScript implementations offer a flexible and performant solution, with libraries such as Math.js [@mathjs] and odex-js [@odexjs], although they require programming proficiency.

Current ODE modeling tools also face other limitations. Commercial platforms are costly that restricts access for underfunded researchers, small labs, and institutions. Sharing simulations often necessitates that recipients have the same software, complicating collaboration and reproducibility. Traditional platforms typically provide static interfaces, and developing interactive applications demands programming expertise, which can shift the focus from scientific exploration to software development.

This project proposes a web-based approach that addresses these challenges. The goal is to develop an ecosystem combining a "no-code" interface with comprehensive capabilities for interactive, in-browser modeling and analysis of ODE-based systems.

# The solution: Diff Studio

The Diff Grok library provides numerical methods for solving problems specified in a declarative form. It includes:

- **Solving tools:** numerical methods and computational pipelines. The library implements the modified Rosenbrock triple (MRT) [@Shampine1997], the ROS3PRw [@jax2021], and the ROS34PRw [@rang2015improved] methods that provide solving of both stiff and non-stiff systems. The performance is benchmarked on **Robertson** [@robertson1966solution], **HIRES** [@schafer1975new], **VDPOL** [@vanderpol1926relaxation], **OREGO** [@hairer2002solving2], **E5** [@hairer2002solving2], and **Pollution** [@verwer1994gauss]. Diff Grok allows users to obtain modeling results in near-real time (see \autoref{fig:peformance}). Computational pipelines support multi-stage modeling and solving of IVPs in web workers, enabling parallel computations and analyses, including parameter optimization and sensitivity studies.
- **Declarative modeling language:** a domain-specific language for specifying problems. It provides the ability to define IVPs as a text containing equations and model input annotations.

![Diff Grok performance: computational time comparison.\label{fig:peformance}](./images/dg-performance.png)

Diff Grok is a foundation for creating web applications for ODE-based modeling. Diff Studio integrates this library with the Datagrok platform [@datagrok]. This web application has an equations editor (\autoref{fig:dseditor}) and provides a model analysis via an autogenerated user interface (\autoref{fig:autoui}).

![Pharmacokinetic-pharmacodynamic simulation with Diff Studio: the equation editor, numerical solution, and its visualization. The view displays the declarative model specification with input annotations.\label{fig:dseditor}](./images/ds-editor.png)

The Datagrok platform provides in-browser computations and visualizations, while interactivity in exploring ODE-based models is a key feature of Diff Studio. Each time the user adjusts sliders or modifies any inputs, the application automatically recalculates the results. The high-performance Diff Grok library performs these computations almost in real time. Another important feature is the ability to collaborate on models and share computations via URLs. The platform also provides sensitivity analysis and parameter optimization features.

Thus, Diff Studio serves as a comprehensive modeling environment and a centralized hub for ODE-based models.

![Diff Studio in model exploration mode, showing the autogenerated UI. Adjusting sliders or modifying any inputs automatically updates the results.\label{fig:autoui}](./images/ds-ui.png)

We plan to implement automatic switching between stiff and non-stiff integration modes. This will allow the solver to apply the most appropriate method: using Rosenbrock for stiff regions and an explicit method for non-stiff regions.

## Availability

Links:

- [Diff Grok](https://github.com/datagrok-ai/diff-grok)
- [Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio)
- [Datagrok Help pages](https://datagrok.ai/help/compute/diff-studio)

Run Diff Studio online [here](https://public.datagrok.ai/apps/DiffStudio), or complete an interactive [tutorial](https://public.datagrok.ai/apps/tutorials/Tutorials/Scientificcomputing/Differentialequations).

# Acknowledgements

The authors are grateful to the entire **Datagrok** team and to the **JnJ ModelHub** project team for their contributions and feedback, which significantly improved the project.

# Conflicts of interest

The authors declare no conflict of interest.

# References
