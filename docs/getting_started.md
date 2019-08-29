---
layout: default
title: Getting started
---
# Getting started

## Requirements

- [GNU Make](https://www.gnu.org/software/make/)

- A Fotran compiler

- MPI (e.g. [MPICH](https://www.mpich.org/) or [Open MPI](https://www.open-mpi.org/)) linked to your Fortran compiler (`mpif90`)

- Python 3+, MATLAB, or Octave. The Python wrapper relies on [NumPy](http://www.numpy.org/) and [Pandas](https://pandas.pydata.org/). The Python test suite (optional but recommended for developers) additionally requires [matplotlib](https://matplotlib.org/) and [termcolor](https://pypi.org/project/termcolor/). These dependencies are conveniently acquired through [`pip`](https://pypi.org/project/pip/):<br />
  ```
  pip install numpy pandas matplotlib termcolor
  ```
  or through [Anaconda](https://www.anaconda.com/):<br />
  ```
  conda create -n QDYN python=3.7 numpy matplotlib pandas termcolor
  conda activate QDYN
  ```

- For Windows 10 users, Linux tools can be acquired and run natively through a [Linux subsystem](https://docs.microsoft.com/en-us/windows/wsl/install-win10)



## Downloading QDYN

QDYN is hosted on [GitHub](https://github.com/ydluo/qdyn). To download for the first time the stable version of QDYN, execute the following git command:

```
git clone https://github.com/ydluo/qdyn qdyn-read-only
```

This creates a directory  `qdyn-read-only` which contains the whole QDYN package. You can create a directory with a different name. The code contained by the `master` branch (default) is tested and stable, but other development branches may be available. Consult the GitHub repository for the availability of
development code.



## Installing QDYN

1. Navigate to the `src` directory
2. Modify the section "User Settings" of the  `Makefile` following the instructions and examples therein:
    *  In section 1, set the variable `EXEC_PATH = [target path to your executable file]`. If you set the default value (recommended) the executable file `qdyn` is placed in the `src` directory. If you change this variable, you must set the `EXEC_PATH` input variable accordingly when calling `qdyn` from the wrappers (`qdyn.m` or `pyqdyn.py`).
    *  In section 2, adjust your Fortran compiler settings: set the variables `F90 = [your compiler]`, `OPT = [your compiler optimization flags]` and `PREPROC = [your compiler preprocessing flags]`. Settings for several commonly used compilers are provided. Note that the specific optimization flags need to be set to enable parallelization through OpenMP.
3. Set the parameters in the section "User Settings" of `constants.f90` following the instructions therein
4. Run `make`




## Keeping QDYN up-to-date

After the first-time checkout you can update the QDYN package by executing the following command in your QDYN directory:

```
git pull origin master
```

Git automatically detects conflicts and attempts to resolve them. In case of unresolvable conflicts you have to fix them manually following the instructions in the [GitHub help pages](https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/).



## Additional notes for Windows 10 users

As of 2017, Windows 10 officially supports a bash command line environment by installing a Linux subsystem (as of writing, Ubuntu and OpenSUSE are currently offered in the Windows Store). Within a subsystem, Unix-compiled executables can be run natively, and the user has access to the Canonical software repository (`apt-get install [package]`). Running QDYN in a Linux subsystem is done as follows:
1. Install your preferred Linux subsystem, see [this instruction page](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

2. Install `make`, `gfortran`, and Open MPI as:  
    ```bash
    sudo apt-get install make gfortran libopenmpi-dev
    ```

3. Download QDYN as instructed above. Note that Windows does not have access to the Linux file system, so in order to exchange files between the subsystem and Windows, it is recommended to download QDYN to (and run simulations from) a local Windows directory (e.g. C:\Users\bob\qdyn). The Windows file system can be accessed in the Linux subsystem as:  `cd /mnt/c/Users/bob/qdyn`

4. Navigate to the QDYN `src` directory and compile QDYN as described above

5. In the case that the required Python or command line MATLAB/Octave packages are installed on the Linux subsystem, QDYN can be called directly from a wrapper. If none of these software packages are available, generate a `qdyn.in` file in Windows (through a wrapper), navigate within the subsystem to the location of  `qdyn.in` (e.g.  `cd /mnt/c/Users/bob/test_simulation`) and run: `/mnt/c/Users/bob/qdyn/src/qdyn`

6. QDYN should now be running within the Linux subsystem, creating output files in C:\Users\bob\test_simulation that can be accessed by Windows for further processing.

7. The Python wrapper (`pyqdyn.py`) also has built-in functionalities to call the subsystem directly from a Windows 10 environment. In order to set-up and run QDYN simulations from the Python wrapper, set
    `qdyn.W10_bash = True`. When doing so, the wrapper will automatically switch between the Windows and Linux environments.