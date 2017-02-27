# QDYN
## A Quasi-DYNamic earthquake simulator

--------------------------------

## News 

*Feb 2017* | **QDYN v1.1 has been released** and published online via [zenodo] (https://zenodo.org/record/322459#.WLNq3BiZNE4) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459). You can cite it as:  
  Y. Luo, J. P. Ampuero, P. Galvez, M. van den Ende and B. Idini (2017)  
  *QDYN: a Quasi-DYNamic earthquake simulator (v1.1)*  
  Zenodo. doi:10.5281/zenodo.322459

*Jan 2017* | In the following paper we used *QDYN* to simulate earthquake cycles and obtained events spanning a broad range of magnitudes to study earthquake scaling relations:  
  Y. Luo, J. P. Ampuero, K. Miyakoshi and K. Irikura (2017)  
  *Surface effects on earthquake moment-area scaling relations*   
  PAGEOPH, Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*  
  doi:10.1007/s00024-017-1467-4  [PDF] (https://rdcu.be/oOL9)

*Nov 2014* | A *QDYN* tutorial session was offered at the [School on "Earthquakes: nucleation, triggering and relationship with aseismic processes"] (http://earthquakes.sciencesconf.org/) at Cargese, Corsica, 3 - 10 November 2014 

*QDYN* was previously hosted on [GoogleCode] (https://code.google.com/p/qdyn/) and is now on [Github] (http://ydluo.github.io/qdyn). Both [SVN] (https://subversion.apache.org) and [GIT] (https://git-scm.com) version control systems are supported. 


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
  * MPI parallelization
  * normal stress coupling
  * Faults surrounded by damaged zones
  * fully coupled with SPECFEM3D via QSB (QDYN-SPECFEM Bridge)


--------------------------------

## Downloads, documentation and support

The primary documentation for QDYN is the [User's manual](https://github.com/ydluo/qdyn/blob/master/doc/QDYN_man_GIT.pdf). 

To download and install QDYN please follow the instructions in _Section 2_ of the manual.

[QDYN Wiki](https://github.com/ydluo/qdyn/wiki)

Submit questions, feedback or suggestions via our [issue tracking system] (https://github.com/ydluo/qdyn/issues).


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

[Yingdi Luo](http://www.seismolab.caltech.edu/luo_y.html) (Caltech Seismolab, USA)

[Jean-Paul Ampuero](http://www.seismolab.caltech.edu/ampuero_jp.html) (Caltech Seismolab, USA)

Percy Galvez (AECOM, Switzerland)

Martijn van den Ende (Utrecht University, The Netherlands)

Benjamin Idini (University of Chile)



-------------------------

## Suggested References:

Y. Luo, J. P. Ampuero (2011), Numerical Simulation of Tremor Migration Triggered by Slow Slip and Rapid Tremor Reversals, AGU Fall Meeting 2011 Abstract S33C-02

Y. Luo, J. P. Ampuero (2012), Simulation of Complex Tremor Migration Patterns, AGU Fall Meeting 2012 Abstract S44B-02

Y. Luo, J. P. Ampuero, K. Miyakoshi and K. Irikura (2017)
  *Surface effects on earthquake moment-area scaling relations*  
  PAGEOPH, Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*
  [doi:10.1007/s00024-017-1467-4] (https://link.springer.com/article/10.1007/s00024-017-1467-4)  
  [PDF] (https://rdcu.be/oOL9)
  
Y. Luo, J. P. Ampuero, P. Galvez, M. Ende and B. Idini. (2017). 
*QDYN: a Quasi-DYNamic earthquake simulator (v1.1) [Data set]*. Zenodo. doi:10.5281/zenodo.322459  
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459)
