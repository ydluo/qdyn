---
layout: default
title: Developer notes
---

# Developer notes

The QDYN team welcomes contributions from the community, which is facilitated by GitHub through [pull requests](https://help.github.com/en/articles/about-pull-requests). The optimal workflow for this is as follows:

- [Fork](https://help.github.com/en/articles/fork-a-repo) the QDYN repository
- Modify or extend the code
- Push the commits to your forked repository
- Create a pull request to merge your code with that in the QDYN repository

When the pull request is made, a GitHub Actions workflow will be triggered to assess the validity of the code by compiling the code and running the testing suite. The core developers will ensure that the pull request will be merged with the appropriate branch (e.g. `release/2.3.4`).

For large modifications or new features, please contact the QDYN team or [open an issue](https://github.com/ydluo/qdyn/issues) on GitHub to discuss the implementation strategy.



## Work in progress

The following items are under active development and will be included in the stable branch in the future.

### Physics

- Surface deformation to compare with GPS data
- Slip-dependent friction law
- 3D kernel for faults in infinite media
- Variable strike (including free surface)
- Layered media: EDKS (Luis Rivera), Relax (Sylvain Barbot)
- Heterogeneous media: import kernel from Relax, Pylith, or SPECFEM3D?
- Triangular mesh elements: e.g. Meade ([2007](https://doi.org/10.1016/j.cageo.2006.12.003)), Gimbutas et al. ([2012](https://doi.org/10.1785/0120120127)), Pan et al. ([2014](https://doi.org/10.1785/0120140161)), Nikkhoo & Walter ([2015](https://doi.org/10.1093/gji/ggv035))

### Code engineering and optimisation

- QDYN-SPECFEM3D (QSB) bridge
  - [ ] Update I/O management for compatibility with the stable QDYN version
  - [ ] Design benchmark/integration test
  - [ ] Create examples/tutorials using QSB
  - [ ] Update documentation
  - [ ] Incorporate QSB into the Python wrapper
- [Intel MKL](https://software.intel.com/en-us/mkl) optimised libraries (e.g. for FFT)
- Hierarchical matrix techniques ("H-matrices")

### Input/Output

- Parallel MPI I/O or parallel I/O library (HDF5, NETCDF4, ADIOS)
- Binary output




## The QDYN team

*(listed alphabetically)*

[Jean-Paul Ampuero](http://www.seismolab.caltech.edu/ampuero_jp.html) (IRD/UCA, Géoazur, France; Caltech Seismolab, USA)

[Martijn van den Ende](https://www.linkedin.com/in/martijnvandenende) (Université Côte d'Azur, Géoazur, France)

[Percy Galvez](https://smi.kaust.edu.sa/Pages/People-Galvez.aspx) (KAUST, Saudi Arabia; AECOM, Switzerland)

[Benjamin Idini](http://www.seismolab.caltech.edu/idini_b.html) (Caltech Seismolab, USA)

[Yingdi Luo](https://science.jpl.nasa.gov/people/YLuo/) (NASA JPL, USA)