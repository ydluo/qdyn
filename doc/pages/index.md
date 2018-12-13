---
layout: default
title: Summary
---

QDYN is a boundary element software to simulate earthquake cycles (seismic and aseismic slip on tectonic faults) under the quasi-dynamic approximation (quasi-static elasticity combined with radiation damping) on faults governed by rate-and-state friction and embedded in elastic media.

QDYN includes various forms of rate-and-state friction and state evolution laws, and handles non-planar fault geometry in 3D and 2D media, as well as spring-block simulations. Loading is controlled by remote displacement, steady creep or oscillatory load. In 3D it handles free surface effects in a half-space, including normal stress coupling. The medium surrounding the fault is linear, isotropic and elastic, and may be uniform or (in 2D) contain a damaged layer.

QDYN implements adaptive time stepping, shared-memory parallelization, and can deal with multi-scale earthquake cycle simulations with fine details in both time and space. It is equipped with a user-friendly Matlab interface and graphical output utilities.

Test



### Main features

- Rate-and-state friction, with velocity cut-offs, aging and slip laws
- Arbitrarily heterogeneous frictional properties
- Slow and fast, aseismic and seismic slip transients (adaptive timestep)
- Non-planar faults (currently limited to variable dip, rectangular elements)
- 3D, 2D and 1D (spring-block)
- Steady and oscillatory loads
- Normal stress coupling 
- Faults surrounded by damaged zones
- Matlab wrapper and graphic output display utilities
- Parallelized for shared memory systems (OpenMP)
- Parallelized for distributed memory systems (MPI)
- Fully coupled with SPECFEM3D via QSB (QDYN-SPECFEM Bridge)

