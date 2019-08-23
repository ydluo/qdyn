---
layout: default
title: Summary
---



# Overview

QDYN is a boundary element software to simulate earthquake cycles (seismic and aseismic slip on tectonic faults) under the quasi-dynamic approximation (quasi-static elasticity combined with radiation damping) on faults governed by rate-and-state friction and embedded in elastic media.

QDYN includes various forms of rate-and-state friction and state evolution laws, and handles non-planar fault geometry in 3D and 2D media, as well as spring-block simulations. Loading is controlled by remote displacement, steady creep or oscillatory load. In 3D it handles free surface effects in a half-space, including normal stress coupling. The medium surrounding the fault is linear, isotropic and elastic, and may be uniform or (in 2D) contain a damaged layer.

QDYN implements adaptive time stepping, shared-memory parallelization, and can deal with multi-scale earthquake cycle simulations with fine details in both time and space. It is equipped with user-friendly MATLAB and Python interfaces and graphical output utilities.



<video width="650" autoplay loop>
    <source src="img/double_asperity.mp4" type="video/mp4">
</video>

<div style="width: 100%; text-align: center;"><em>(evolution of slip velocity on a fault with two seismogenic asperities)</em></div>


## Main features

- Rate-and-state friction, with velocity cut-offs, aging and slip laws
- Microphysically based frictional model (*Chen-Niemeijer-Spiers* model)
- Heterogeneous frictional properties
- Slow and fast, aseismic and seismic slip transients
- Dynamic weakening (thermal pressurization)
- Non-planar faults (currently limited to variable dip, rectangular elements)
- 3D, 2D and 1D (spring-block)
- Tectonic and transient loads
- Normal stress coupling
- Faults surrounded by damaged zones
- MATLAB and python wrappers, and graphic output display utilities
- Parallelized for shared memory systems (OpenMP)
- Parallelized for distributed memory systems (MPI)
- Fully coupled with SPECFEM3D via QSB (QDYN-SPECFEM Bridge)



## Support

The QDYN development team offers online support to users who report bugs, installation problems, documentation issues, feature requests, or questions about QDYN usage by submitting "issues" via the [GitHub issue system](https://github.com/ydluo/qdyn/issues). Please do not contact the QDYN developers directly by email. 

Before submitting an issue please make sure that:

- you have read the QDYN documentation
- you are running the most recent stable version of QDYN
- your problem has not been treated in previous issues. You can browse and search the list of closed issues

When submitting a new issue, please include all information needed to reproduce your problem: input files, operating system, compiler, QDYN version (git hash).



## Acknowledgements

### Code contributions

QDYN started from a 2D code written by Allan Rubin (Princeton University) in the early 2000s (*Rubin & Ampuero*, [2005](http://dx.doi.org/10.1029/2005JB003686), [2009](http://dx.doi.org/10.1029/2009JB006529); *Ampuero & Rubin*, [2008](http://dx.doi.org/10.1029/2007JB005082)). The main developers of QDYN are Jean-Paul Ampuero, Martijn van den Ende, Percy Galvez, Benjamin Idini and Yingi Luo. Bryan Riel contributed the double-FFT version.

The subroutines implementing Okada’s formulas were provided by Shinichi Miyazaki (Kyoto University). They include subroutines written by Y. Okada. The FFT subroutines are based on the  General Purpose FFT Package written by Takuya Ooura (Kyoto University).



### Funding

The development of QDYN has been supported by the US National Science Foundation, the Southern California Earthquake Center, Japan’s Nuclear Regulation Authority (formerly Japan Nuclear Energy Safety Organization), and the European Research Council, and the French National Research Agency (UCAJEDI Investments in the Future project ANR-15-IDEX-01).



## License

This software is freely available for academic research purposes. If you use QDYN please include proper attributions to its authors and cite one of the references in section 1.7 in your scientific papers and reports.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


