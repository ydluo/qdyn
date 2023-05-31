# QDYN <img src="docs/img/qdyn_logo_small.jpeg" alt="QDYN logo" align="right" />

## A Quasi-DYNamic earthquake simulator

[![GitHub Action Status](https://github.com/ydluo/qdyn/actions/workflows/tests/badge.svg)](https://github.com/ydluo/qdyn/actions) [![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://ydluo.github.io/qdyn/) [![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://perso.crans.org/besson/LICENSE.html)

[**Quickstart**](#quickstart)
| [**Documentation**](hhttps://ydluo.github.io/qdyn/)
| [**Installation**](https://ydluo.github.io/qdyn/getting_started.html)
| [**Suggested references**](#suggested-references)

------

## What is QDYN?

*QDYN* is a boundary element software to simulate earthquake cycles (seismic and aseismic slip on tectonic faults) under the quasi-dynamic approximation (quasi-static elasticity combined with radiation damping) on faults governed by rate-and-state friction and embedded in elastic media.

*QDYN* includes various forms of rate-and-state friction and state evolution laws, and handles non-planar fault geometry in 3D and 2D media, as well as spring-block simulations. Loading is controlled by remote displacement, steady creep or oscillatory load. In 3D it handles free surface effects in a half-space, including normal stress coupling. The medium surrounding the fault is linear, isotropic and elastic, and may be uniform or (in 2D) contain a damaged layer.

*QDYN* implements adaptive time stepping, shared-memory parallelization, and can deal with multi-scale earthquake cycle simulations with fine details in both time and space. It is equipped with user-friendly MATLAB and Python interfaces and graphical output utilities.

## Main features

- Rate-and-state friction, with velocity cut-offs, aging and slip laws
- Microphysically based frictional model (*Chen-Niemeijer-Spiers* model)
- Heterogeneous frictional properties
- Slow and fast, aseismic and seismic slip transients
- Dynamic weakening (thermal pressurization)
- Multiple and non-planar faults (rectangular elements)
- 3D, 2D and 1D (spring-block)
- Tectonic and transient loads
- Normal stress coupling (free surface effects and fault interactions)
- Faults surrounded by damaged zones
- Python wrapper and graphic output display utilities
- Parallelized for shared memory systems (OpenMP)
- Parallelized for distributed memory systems (MPI)


## Quickstart

Download and install the latest version of QDYN directly from GitHub:
```
git clone git@github.com:ydluo/qdyn.git
```

Install the Python wrapper with pip:
```
cd qdyn
pip install -e .
```
This will add the QDYN wrapper to your Python path, and can be imported as:
```
from qdyn import qdyn
```

Compile the Fortran code:
```
cd qdyn
make clean && make
```

Try out one of the example Jupyter notebooks in the `/examples/notebooks` directory.

For a more detailed installation guide (including requirements), see the [*Getting started*](https://ydluo.github.io/qdyn/getting_started.html) section in the documentation.

Previous (stable) versions of QDYN can be downloaded from the [release page](https://github.com/ydluo/qdyn/releases). Development versions are available as separate [branches](https://github.com/ydluo/qdyn/branches) following the naming convention `release/x.x.x`.



Questions, feedback or suggestions can be submitted via our [issue tracking system](https://github.com/ydluo/qdyn/issues).

## Core developers

*(listed alphabetically)*

[Jean-Paul Ampuero](http://www.seismolab.caltech.edu/ampuero_jp.html) (IRD/UCA, Géoazur, France; Caltech Seismolab, USA)

[Martijn van den Ende](https://www.linkedin.com/in/martijnvandenende) (Université Côte d'Azur, Géoazur, France)

[Percy Galvez](https://smi.kaust.edu.sa/Pages/People-Galvez.aspx) (KAUST, Saudi Arabia; AECOM, Switzerland)

[Benjamin Idini](http://www.seismolab.caltech.edu/idini_b.html) (Caltech Seismolab, USA)

[Yingdi Luo](https://science.jpl.nasa.gov/people/YLuo/) (NASA JPL, USA)

## Suggested References

### For all uses of the QDYN software

Luo, Y., Ampuero, J. P., Galvez,  P., van den Ende, M., & Idini, B. (2017). 
QDYN: a Quasi-DYNamic earthquake simulator (v1.1) [Data set]. Zenodo. doi:10.5281/zenodo.322459  
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459)

### For simulations on heterogeneous faults

Luo, Y., & Ampuero, J. P. (2018). 
Stability of faults with heterogeneous friction properties and effective normal stress. 
Tectonophysics, 733, 257-272, doi:[10.1016/j.tecto.2017.11.006](https://doi.org/10.1016/j.tecto.2017.11.006)

Luo, Y., Ampuero, J. P., Miyakoshi, K., & Irikura, K. (2017). Surface rupture effects on earthquake moment-area scaling relations. Pure and Applied Geophysics, 174(9), 3331-3342, doi:[10.1007/s00024-017-1467-4](https://link.springer.com/article/10.1007/s00024-017-1467-4).
In Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*.  
[PDF](https://rdcu.be/oOL9)

Luo, Y., & Ampuero, J. P. (2012), Simulation of Complex Tremor Migration Patterns, AGU Fall Meeting 2012 Abstract S44B-02

Luo, Y., & Ampuero, J. P. (2011), Numerical Simulation of Tremor Migration Triggered by Slow Slip and Rapid Tremor Reversals, AGU Fall Meeting 2011, abstract S33C-02

### For the microphysically based (CNS) simulations

van den Ende, M. P. A., Chen, J., Ampuero, J. P., & Niemeijer, A. R. (2018).
A comparison between rate-and-state friction and microphysical models, based on numerical simulations of fault slip.
Tectonophysics, 733, 273-295, doi:[10.1016/j.tecto.2017.11.040](https://doi.org/10.1016/j.tecto.2017.11.040)

### For simulations on faults surrounded by damaged zones

Idini, B., & Ampuero, J. P. (2017).
Rupture complexity promoted by damaged fault zones in earthquake cycle models.
AGU Fall Meeting 2017, abstract T41C-0632.  
Poster [PDF](https://www.essoar.org/doi/abs/10.1002/essoar.10500080.1), doi:[10.1002/essoar.10500080.1](https://dx.doi.org/10.1002/essoar.10500080.1)

### For rate-and-state friction simulations including viscosity

Beall, A., van den Ende, M. P. A., Ampuero, J. P., Capitanio, F. A., Fagereng, A. (2022).
Linking Earthquake Magnitude‐Frequency Statistics and Stress in Visco‐Frictional Fault Zone Models.
Geophysical Research Letters 49(20), doi:[10.1029/2022GL099247](https://doi.org/10.1029/2022GL099247)