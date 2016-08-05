# qdyn
## Welcome to QDYN
## A Quasi-DYNmatic earthquake simulator

--------------------------------
## [Access to QDYN Wiki](https://github.com/ydluo/qdyn/wiki)
--------------------------------

## Downloads

[User's manual] (https://docs.google.com/document/d/11FkfCGAwbjfrn9U07EkIZRNuXmvY6PJt67DFNQU_mW8/pub) online (most up-to-date version) or [PDF] (https://drive.google.com/file/d/0B13IeWPRjYq8WFFTUk1FemRSdzA/view?usp=sharing).

[QDYN software and examples] (https://drive.google.com/folderview?id=0B13IeWPRjYq8c2h5SGFsWDVyZ1k&usp=sharing).

--------------------------------

## News 

August 2016: We submitted the following paper, in which we used *QDYN* to simulate earthquake cycles, with events spanning a broad range of magnitudes, to study earthquake scaling relations:  
  Y. Luo, J. P. Ampuero and K. Irikura  
  *Free surface effects on earthquake moment-area scaling relations*  
  PAGEOPH, Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*

*QDYN* was previously hosted on [GoogleCode] (https://code.google.com/p/qdyn/) and is now on [Github] (http://ydluo.github.io/qdyn). Both [SVN] (https://subversion.apache.org) and [GIT] (https://git-scm.com) version control systems are supported. 

A *QDYN* tutorial session was offered at the [School on "Earthquakes: nucleation, triggering and relationship with aseismic processes"] (http://earthquakes.sciencesconf.org/) at Cargese, Corsica, 3 - 10 November 2014 

--------------------------------


## Summary
*QDYN* is a boundary element software package to simulate earthquake cycles (tectonic fault slip) under the quasi-dynamic approximation (quasi-static elasticity with radiation damping).  

The code implements adaptive time stepping and shared-memory parallelization to simulate earthquake cycles including seismic and aseismic slip. QDYN includes various forms of rate-and-state friction and state evolution laws. It handles non-planar fault geometries in 3D and 2D elastic media, as well as spring-block simulations. It has a user-friendly matlab interface and graphical output.

--------------------------------


## Features

  * rate-and-state friction, with velocity cut-offs, aging and slip laws
  * heterogeneous frictional properties
  * slow and fast, aseismic and seismic slip transients
  * non-planar faults (currently limited to variable dip, rectangular elements)
  * 3D, 2D and 1D (spring-block)
  * tectonic and transient loads
  * matlab wrapper and graphic output display utilities
  * parallelized for shared memory systems (OpenMP)
  * normal stress coupling
  * fully coupled with SPECFEM3D via QSB (QDYN-SPECFEM Bridge)




-------------------------

## Introduction Poster

![](https://lh4.googleusercontent.com/-OjKBE5_Ipf8/T9wk2GtVRXI/AAAAAAAAABg/a1diUWu7tFU/s763/Poster_QDYN.jpg)

[Download Poster] (http://code.google.com/p/qdyn/downloads/detail?name=Poster_QDYN.pdf) 

-------------------------


## Featured Simulations

### [Simulation_Tohoku](https://github.com/ydluo/qdyn/wiki/Simulation_Tohoku)
![] (https://lh5.googleusercontent.com/-JPaTpBXo5eA/USdSArzQ0QI/AAAAAAAAKew/9wnVu30Lhf4/s900/Tohoku_cycle_logo.gif)

### [Simulation_Cascadia_Tremor](https://github.com/ydluo/qdyn/wiki/Simulation_Cascadia_Tremor)
![] (https://lh5.googleusercontent.com/-a_2MRxcUgf8/T-v2JCjmxBI/AAAAAAAAAB8/NlQTwfra4fY/s900/Tremor_3D_Cascadia.gif)


------------------------
## Developers

QuakeID Group, Caltech Seismo Lab:

Yingdi Luo [http://www.seismolab.caltech.edu/luo_y.html] (http://www.seismolab.caltech.edu/luo_y.html)

Jean-Paul Ampuero [http://www.seismolab.caltech.edu/ampuero_jp.html] (http://www.seismolab.caltech.edu/ampuero_jp.html)

Bryan Riel [http://www.seismolab.caltech.edu/riel_b.html] (http://www.seismolab.caltech.edu/riel_b.html)

For questions, feedback or suggestions, please contact the developers via the [Issue] (https://github.com/ydluo/qdyn/issues) on GITHUB.

-------------------------

## Suggested References:

Y. Luo, J. P. Ampuero (2011), Numerical Simulation of Tremor Migration Triggered by Slow Slip and Rapid Tremor Reversals, AGU Fall Meeting 2011 Abstract S33C-02

Y. Luo, J. P. Ampuero (2012), Simulation of Complex Tremor Migration Patterns, AGU Fall Meeting 2012 Abstract S44B-02


