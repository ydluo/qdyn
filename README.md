# QDYN
## A Quasi-DYNamic earthquake simulator

--------------------------------

## News 

*21 August 2019* | QDYN version 2.1 has been released! See the [release notes](https://github.com/ydluo/qdyn/releases/tag/qdyn_2.1) for new features and (minor) bug fixes.

*11 December 2018* | QDYN version 2.0 has been released!

*6 December 2018* | Results of QDYN earthquake cycle simulations on faults surrounded by damaged zones were presented at:  
* the AGU Fall Meeting 2017 (San Francisco, 11-15 December 2017): B. Idini and J. P. Ampuero, *Rupture complexity promoted by damaged fault zones in earthquake cycle models* (doi:10.1002/essoar.10500080.1, [poster](https://www.essoar.org/doi/abs/10.1002/essoar.10500080.1))  
* the [10th ACES International Workshop](http://quaketm.bosai.go.jp/~shiqing/ACES2018/index_aces.html) (Japan, September 25-28, 2018), J. P. Ampuero et al., *Rupture complexity promoted by damaged fault zones in earthquake cycle models* ([slides](http://quaketm.bosai.go.jp/~shiqing/ACES2018/abstracts/aces_abstract_ampuero.pdf)) 
* the [Workshop on Modeling Earthquake Source Processes](http://www.seismolab.caltech.edu/workshop.html) (Caltech, October 8-10, 2018). 

*29 November 2018* | The QDYN team participated in  in the two first [benchmarks](https://scecdata.usc.edu/cvws/seas/index.html) and [workshops](https://www.scec.org/workshops/2018/seas) of the SCEC SEAS project (Southern California Earthquake Center, Sequences of Earthquakes and Aseismic Slip). Yingdi Luo presented work using QDYN on *A heterogenous fault model of episodic tremor and slow-slip events with spatial-temporal variability* ([slides](https://files.scec.org/s3fs-public/2018_SEAS_Workshop_1115_Luo_A_heterogenous_fault_model_of_episodic_tremor_and_slow-slip_event_with_spatial-temporal_variability.pdf)).

*11 September 2018* | Results based on QDYN were highlighted at the [2018 WEGENER conference](https://wegener2018.sciencesconf.org/) in the keynote presentation *The spectrum of slip behaviors emerging from the interactions between seismic and aseismic slip* by J. P. Ampuero ([slides](https://wegener2018.sciencesconf.org/data/pages/Jean_Paul_Ampuero_Wegener_2018.pdf)).

*28 May 2018* | In order to address a problem in the QDYN history tree, the version history of the QDYN repository was rewritten. All developers that have cloned the QDYN repository before this date should follow the steps described in [this manual](doc/git_fix_2018-05-28.pdf). Users of QDYN that did not make any code changes can simply delete the old repository and clone the latest version to get a clean QDYN repository.

*May 2018* | Two papers with simulation results based on *QDYN*, published in Tectonophysics - Special Issue on "[Physics of Earthquake Rupture Propagation](https://www.sciencedirect.com/journal/tectonophysics/vol/733/suppl/C)":  
  M. van den Ende, J. Chen, J. P. Ampuero and A. Niemeijer
  *A comparison between rate-and-state friction and microphysical models, based on numerical simulations of fault slip* (doi:[10.1016/j.tecto.2017.11.040](https://doi.org/10.1016/j.tecto.2017.11.040))  
  Y. Luo and J. P. Ampuero  
  *Stability and effective friction of faults with heterogeneous friction properties and fluid pressure* (doi:[10.1016/j.tecto.2017.11.006](https://doi.org/10.1016/j.tecto.2017.11.006))

*Apr 2017* | Earthquake cycle simulations with a micro-physics model of granular flow and ductile creep of fault gouges based on *QDYN* presented at the EGU General Assembly 2017, session on *Earthquakes: from slow to fast, from the field to the laboratory*:    
  M. van den Ende, J. Chen, J. P. Ampuero and A. Niemeijer (2017)  
  *Earthquake and slow-slip nucleation investigated with a micro-physics based seismic cycle simulator*  
  Geophysical Research Abstracts, Vol. 19, EGU2017-7249  
  Abstract [PDF](http://meetingorganizer.copernicus.org/EGU2017/EGU2017-7249.pdf) -- Presentation [PPT](http://presentations.copernicus.org/EGU2017-7249_presentation.pptx)  

*Feb 2017* | [**QDYN v1.1 has been released**](https://github.com/ydluo/qdyn/releases/tag/qdyn_1.1) and published online via [zenodo](https://zenodo.org/record/322459#.WLNq3BiZNE4). This new version introduces MPI parallelization for 3D simulations in HPC clusters and faults surrounded by damaged zones in 2D. You can cite it as:  
  Y. Luo, J. P. Ampuero, P. Galvez, M. van den Ende and B. Idini (2017)  
  *QDYN: a Quasi-DYNamic earthquake simulator (v1.1)*  
  Zenodo. doi:10.5281/zenodo.322459  
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459)

*Jan 2017* | In the following paper we used *QDYN* to simulate earthquake cycles and obtained events spanning a broad range of magnitudes to study earthquake scaling relations:  
  Y. Luo, J. P. Ampuero, K. Miyakoshi and K. Irikura (2017)  
  *Surface effects on earthquake moment-area scaling relations*     
  PAGEOPH, Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*  [doi:10.1007/s00024-017-1467-4](https://link.springer.com/article/10.1007/s00024-017-1467-4)  
  [PDF](https://rdcu.be/oOL9)  

*Nov 2014* | A *QDYN* tutorial session was offered at the [School on "Earthquakes: nucleation, triggering and relationship with aseismic processes"](http://earthquakes.sciencesconf.org/) at Cargese, Corsica, 3 - 10 November 2014 

*QDYN* was previously hosted on [GoogleCode](https://code.google.com/p/qdyn/) and is now on [Github](http://ydluo.github.io/qdyn). Both [SVN](https://subversion.apache.org) and [GIT](https://git-scm.com) version control systems are supported. 


--------------------------------

## Summary

*QDYN* is a boundary element software to simulate earthquake cycles (tectonic fault slip) under the quasi-dynamic approximation (quasi-static elasticity with radiation damping).  

*QDYN* has adaptive time stepping to simulate earthquake cycles including seismic and aseismic slip, on non-planar faults in 3D and 2D elastic media and on spring-block systems. It includes various forms of rate-and-state friction and state evolution laws, a microphysical model implementation, dynamic weakening (thermal pressurization), parallelization for shared and distributed memory systems, and a user-friendly matlab or python interface for setting inputs and visualizing outputs.

--------------------------------

## Features

  * rate-and-state friction, with velocity cut-offs, aging and slip laws
  * microphysically based frictional model (Chen-Niemeijer-Spiers model)
  * heterogeneous frictional properties
  * slow and fast, aseismic and seismic slip transients
  * dynamic weakening (thermal pressurization)
  * non-planar faults (currently limited to variable dip, rectangular elements)
  * 3D, 2D and 1D (spring-block)
  * tectonic and transient loads
  * normal stress coupling
  * faults surrounded by damaged zones
  * matlab and python wrappers, and graphic output display utilities
  * parallelized for shared memory systems (OpenMP)
  * parallelized for distributed memory systems (MPI)
  * fully coupled with SPECFEM3D via QSB (QDYN-SPECFEM Bridge)


--------------------------------

## Downloads, documentation and support

The primary documentation for QDYN is the [User's manual](https://github.com/ydluo/qdyn/blob/master/doc/QDYN_man_GIT.pdf). 

Stable versions of QDYN can be downloaded from the [release page](https://github.com/ydluo/qdyn/releases). To download development versions please follow the instructions in _Section 2.2_ of the manual.

To install QDYN please follow the instructions in _Section 2.3_ of the manual. 

Questions, feedback or suggestions can be submitted via our [issue tracking system](https://github.com/ydluo/qdyn/issues).


-------------------------

## Introduction Poster

![](https://lh4.googleusercontent.com/-OjKBE5_Ipf8/T9wk2GtVRXI/AAAAAAAAABg/a1diUWu7tFU/s763/Poster_QDYN.jpg)

[Download Poster](http://code.google.com/p/qdyn/downloads/detail?name=Poster_QDYN.pdf) 

-------------------------


## Featured Simulations

### [Simulation_Tohoku](https://github.com/ydluo/qdyn/wiki/Simulation_Tohoku)
![](https://lh5.googleusercontent.com/-JPaTpBXo5eA/USdSArzQ0QI/AAAAAAAAKew/9wnVu30Lhf4/s900/Tohoku_cycle_logo.gif)

### [Simulation_Cascadia_Tremor](https://github.com/ydluo/qdyn/wiki/Simulation_Cascadia_Tremor)
![](https://lh5.googleusercontent.com/-a_2MRxcUgf8/T-v2JCjmxBI/AAAAAAAAAB8/NlQTwfra4fY/s900/Tremor_3D_Cascadia.gif)


------------------------
## Developers

[Yingdi Luo](https://science.jpl.nasa.gov/people/YLuo/) (UCLA / NASA JPL, USA)

[Jean-Paul Ampuero](http://www.seismolab.caltech.edu/ampuero_jp.html) (Caltech Seismolab, USA and IRD / UCA, Géoazur, France)

[Percy Galvez](https://smi.kaust.edu.sa/Pages/People-Galvez.aspx) (KAUST, Saudi Arabia; AECOM, Switzerland)

[Martijn van den Ende](https://www.linkedin.com/in/martijnvandenende) (Université Côte d'Azur, Géoazur, France)

[Benjamin Idini](http://www.seismolab.caltech.edu/idini_b.html) (Caltech Seismolab, USA)


-------------------------

## Suggested References

#### For all uses of the QDYN software

Luo, Y., Ampuero, J. P., Galvez,  P., van den Ende, M., & Idini, B. (2017). 
QDYN: a Quasi-DYNamic earthquake simulator (v1.1) [Data set]. Zenodo. doi:10.5281/zenodo.322459  
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459)

#### For simulations on heterogeneous faults

Luo, Y., & Ampuero, J. P. (2018). 
Stability of faults with heterogeneous friction properties and effective normal stress. 
Tectonophysics, 733, 257-272, doi:[10.1016/j.tecto.2017.11.006](https://doi.org/10.1016/j.tecto.2017.11.006)

Luo, Y., Ampuero, J. P., Miyakoshi, K., & Irikura, K. (2017). Surface rupture effects on earthquake moment-area scaling relations. Pure and Applied Geophysics, 174(9), 3331-3342, doi:[10.1007/s00024-017-1467-4](https://link.springer.com/article/10.1007/s00024-017-1467-4).
In Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*.  
[PDF](https://rdcu.be/oOL9)

Luo, Y., & Ampuero, J. P. (2012), Simulation of Complex Tremor Migration Patterns, AGU Fall Meeting 2012 Abstract S44B-02

Luo, Y., & Ampuero, J. P. (2011), Numerical Simulation of Tremor Migration Triggered by Slow Slip and Rapid Tremor Reversals, AGU Fall Meeting 2011, abstract S33C-02

#### For the microphysically based (CNS) simulations

van den Ende, M. P. A., Chen, J., Ampuero, J. P., & Niemeijer, A. R. (2018).
A comparison between rate-and-state friction and microphysical models, based on numerical simulations of fault slip.
Tectonophysics, 733, 273-295, doi:[10.1016/j.tecto.2017.11.040](https://doi.org/10.1016/j.tecto.2017.11.040)
 
#### For simulations on faults surrounded by damaged zones
  
B. Idini, B., & Ampuero, J. P. (2017).
Rupture complexity promoted by damaged fault zones in earthquake cycle models.
AGU Fall Meeting 2017, abstract T41C-0632.  
Poster [PDF](https://www.essoar.org/doi/abs/10.1002/essoar.10500080.1), doi:10.1002/essoar.10500080.1
