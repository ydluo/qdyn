---
layout: default
title: Optimizing performance
---

# Optimizing performance

## Running simulations outside the wrapper

To run simulations outside the MATLAB/Python environment, e.g. when computing on an HPC cluster: run first the wrapper only to generate the input file `qdyn.in`, then run the `qdyn` executable outside the wrapper.

## Managing parallel computing

### OpenMP

For 3D simulations on 2D faults (`MESHDIM=2` or `4`), QDYN is parallelized for shared memory multi-processor systems with OpenMP. Before compiling the code, you should set the specific compiler optimisation flags that enable OpenMP, as described in the  Makefile (see step 2.b in section 2.3). Before performing parallel simulations, you should set the following environment variable:

```
setenv OMP_NUM_THREADS 8
```

This command allows QDYN to run on 8 threads, which will roughly speed up calculations by a factor of 8. The number of threads should be set according to demand. In general, set this value to the maximum number of threads available on your shared memory system.

### MPI

For 3D simulations with `MESHDIM=2`, QDYN can run in parallel in distributed memory clusters with MPI. The number of processors must be set in the variable `NPROCS`.

### Managing outputs of large simulations

QDYN by default outputs results as a single ASCII file (`fort.19`). In most multi-cycle 3D simulations this file is very large. It is then helpful to set `OX_SEQ = 1` in the wrapper, to generate separate `ox` files outputs for each snapshot (`fort.1001, ...`). Also, setting `OX_DYN = 1` will automatically detect seismic events
(according to parameters `DYN_TH_ON` and `DYN_TH_OFF`) and generate 3 snapshots for each event.