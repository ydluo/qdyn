"""
                    .: PyQDYN: Python wrapper for QDYN :.

Author: Martijn van den Ende
Modification date: 2017/07/24

This wrapper can be imported and called from a different Python script file to
do the heavy lifting regarding the generation and population of the mesh,
writing of the qdyn input file, and reading of the qdyn output files. See the
examples directory for usage.

As opposed to the MATLAB wrapper, this Python wrapper heavily uses dictionaries
to contain the input parameters. These parameters are set to default values,
but should be overridden in the calling script file.
Like in the MATLAB wrapper, all parameters are capitalised.

TODO:
- Rendering of 2D mesh for 3D simulations is not tested

"""

from __future__ import print_function

import gzip
import os
import pickle
import re
from subprocess import call

import numpy as np
import pandas as pd
from pandas import read_csv


# import antigravity	# xkcd.com/353/


class qdyn:

    mesh_rendered = False
    qdyn_input_written = False
    mesh_dict = {}

    # Directory path to qdyn executable
    # Default: same directory as this python script
    qdyn_path = os.path.abspath(
        os.path.join(os.path.realpath(__file__), os.pardir)
    )
    # Working directory can be kept empty, except for special cases
    work_dir = ""
    # Flag for using the bash environment in Windows 10
    W10_bash = False

    # Initialise qdyn wrapper
    def __init__(self):
        set_dict = {}

        set_dict["FRICTION_MODEL"] = "RSF"	# Model for fault rheology: "RSF" = rate-and-state friction, "CNS" = Chen-Niemeijer-Spiers microphysical model

        # Fault properties and features
        set_dict["MESHDIM"] = 0				# Dimensionality of fault: 0 = spring-slider, 1 = 1D fault in 2D medium, 2 = 2D fault in 3D medium

        # Along strike boundary conditions:
        # 0 = periodic
        # 1 = fixed-length asperity surrounded by steady creep
        # 2 = periodic and symmetric around the first element
        # 3 = fixed-length and symmetric around the first element
        set_dict["FINITE"] = 1

        # ODE solver mode
        # 1 = Bulirsch-Stoer
        # 2 = Runge-Kutta-Fehlberg
        set_dict["SOLVER"] = 1

        # Fault loading geometry
        #  1 = strike-slip (right-lateral)
        # -1 = strike-slip (left-lateral)
        #  2 = thrust
        # -2 = normal
        set_dict["FAULT_TYPE"] = 1

        set_dict["V_PL"] = 1e-9             # Tectonic loading rate [m/s]
        set_dict["L"] = 1					# Length of the fault if MESHDIM = 1. Stiffness is MU/L if MESHDIM = 0
        set_dict["W"] = 1					# Distance between load point and fault if MESHDIM = 1 and FINITE = 0
        set_dict["SIGMA"] = 50e6            # Normal stress [Pa]
        set_dict["SIGMA_CPL"] = 0			# Normal stress coupling (only for dipping faults)
        set_dict["VS"] = 3000				# Shear wave speed [m/s]. If VS = 0 radiation damping is off
        set_dict["MU"] = 30e9				# Shear modulus [Pa]
        set_dict["LAM"] = 30e9				# Elastic modulus for 3D simulations [Pa]
        set_dict["TPER"] = 31536000			# Period of additional time-dependent oscillatory shear stress loading [s]
        set_dict["APER"] = 0				# Amplitude of additional time-dependent oscillatory shear stress loading [Pa]
        set_dict["DIP_W"] = 0				# Fault dip in 3D
        set_dict["Z_CORNER"] = 0			# 3D

        # Optional simulation features
        set_dict["FEAT_STRESS_COUPL"] = 0	# Normal stress coupling
        set_dict["FEAT_TP"] = 0				# Thermal pressurisation
        set_dict["FEAT_LOCALISATION"] = 0	# Gouge zone localisation of strain (CNS only)

        # Rate-and-state friction parameters
        set_dict["SET_DICT_RSF"] = {
            "RNS_LAW": 0,						# Type of rate-and-state friction law (0: original, 1: with cut-off velocities)
            "THETA_LAW": 1,						# State evolution law (0: aging in no-healin approximation, 1: ageing law, 2: slip law)
            "A": 0.0010,						# Direct effect parameter
            "B": 0.0011,						# Evolution effect parameter
            "DC": 1e-5,							# Characteristic slip distance
            "V1": 0.01,							# Cut-off velocity of direct effect (when RNS_LAW = 1) [m/s]
            "V2": 1e-7,							# Cut-off velocity of evolution effect (when RNS_LAW = 1) [m/s], Controls transition from weakening to strengthening. V2 should be <= V1
            "MU_SS": 0.6,						# Reference steady-state friction
            "V_SS": 1e-6,						# Reference steady-state slip velocity [m/s]
            "CO": 0,							# Cohesion [Pa]
            "V_0": 1.01e-5,						# Initial slip velocity [m/s]
            "TH_0": 1.0,						# Initial state [s]
        }

        # CNS microphysical model parameters
        set_dict["SET_DICT_CNS"] = {
            "THICKNESS": 1e-3,					# Total thickness of gouge zone
            # Granular flow parameters
            "A_TILDE": 0.006,					# Coefficient of logarithmic rate-dependence of grain-boundary friction
            "MU_TILDE_STAR": 0.4,				# Reference grain-boundary friction at Y_GR_STAR
            "Y_GR_STAR": 30e-10,				# Reference strain rate corresponding with MU_TILDE_STAR
            "H": 0.577,							# Dilatancy geometric factor
            # Porosity parameters
            "PHI_C": 0.32,						# Critical state porosity
            "PHI0": 0.03,						# Lower cut-off velocity (e.g. percolation threshold)
            "PHI_INI": 0.25,					# Simulation initial porosity
            # Creep mechanism parameters (grouped in list)
            "N_CREEP": -1,                      # Number of creep mechanisms (must be set by render_mesh() )
            "A": [0.0],				            # (T-dependent) kinetic parameters for creep mechanisms
            "N": [1.0],                         # Stress exponents for creep mechanisms
            "M": [2.0],				            # Porosity exponents for creep mechanisms
            # State of stress
            "TAU": 2e7,							# Initial shear stress [Pa]
        }

        # Localisation settings
        set_dict["SET_DICT_LOCALISATION"] = {
            "LOCALISATION": 1.0,				# Degree of localisation: 1.0 = entire gouge layer is active, 0.0 = entire gouge zone is passive
            "PHI_INI_BULK": 0.2,				# Initial porosity passive gouge zone
            # Creep mechanism parameters for passive/bulk zone (grouped in list)
            "A": [0.0],				            # (T-dependent) kinetic parameters for creep mechanisms
            "N": [1.0],                         # Stress exponents for creep mechanisms
            "M": [2.0],				            # Porosity exponents for creep mechanisms
        }

        # Thermal pressurisation model
        set_dict["SET_DICT_TP"] = {
            "RHOC": 2.16e6,						# Density times specific heat capacity (rho*cp) [J/k/m3]
            "HALFW": 0.5e-3,					# Half-width of the fault zone [m]
            "BETA": 2e-9,						# Bulk compressibility (fluid + pore) [1/Pa]
            "ETA": 2e-4,						# Fluid dynamic viscosity [Pa s]
            "LAM": 2.78e-4,						# Nett thermal expansion coefficient (fluid - solid) [1/K]
            "K_T": 2.0,							# Thermal conductivity [J/s/K/m]
            "K_P": 1e-16,						# Intrinsic hydraulic permeability [m2]
            "P_A": 0.0,							# Ambient fluid pressure [Pa]
            "T_A": 293.0,						# Ambient temperature [K]
            "DILAT_FACTOR": 0.0,                # Factor > 0 to control amount of dilatancy hardening
        }

        # Benjamin Idini's damage model
        set_dict["D"] = 0
        set_dict["HD"] = 0

        # Discretisation and accuracy parameters
        set_dict["N"] = -1					# Number of mesh elements if MESHDIM = 1
        set_dict["NX"] = 1					# Number of fault elements along-strike in 3D
        set_dict["NW"] = -1					# Number of fault elements along-dip in 3D
        set_dict["DW"] = 1					# Length of segment along-strike in 3D
        set_dict["TMAX"] = 1000				# Total simulation time [s]
        set_dict["NSTOP"] = 0				# Stopping criterion (0: stop at t=TMAX, 1: stop end localisation phase, 2: stop first slip rate peak, 3: stop v > TMAX)
        set_dict["DTTRY"] = 1e-1			# First trial step [s]
        set_dict["DTMAX"] = 0				# Maximum time step [s] (0 = unrestricted)
        set_dict["ACC"] = 1e-7				# Solver accuracy

        # Output control parameters
        set_dict["V_TH"] = 1e-2				# Threshold velocity for seismic event
        set_dict["NTOUT_OT"] = 1			# Temporal interval (number of time steps) for time series output
        set_dict["NTOUT"] = 1				# Temporal interval (number of time steps) for snapshout output
        set_dict["NXOUT"] = 1				# Spatial interval (number of elements) for snapshot output
        set_dict["OX_SEQ"] = 0				# Type of snapshot outputs (0: all snapshots in single file, 1: one file per snapshot)
        set_dict["OX_DYN"] = 0				# Output specific snapshots of dynamic events defined by thresholds on peak slip velocity DYN_TH_ON and DYN_TH_OFF
        set_dict["NXOUT_DYN"] = 1			# Spatial interval (number of elements) for dynamic snapshot outputs
        set_dict["DYN_TH_ON"] = 1e-3		# peak slip rate threshold defining beginning of dynamic event
        set_dict["DYN_TH_OFF"] = 1e-4		# peak slip rate threshold defining end of dynamic event
        set_dict["IC"] = 0					# Index of selected elements for time series output (starting at 0)
        set_dict["IOT"] = 0					# Indices of elements for additional time series outputs
        set_dict["IASP"] = 0				# Flags for elements (identification only)
        set_dict["SUFFIX"] = ""

        # SPECFEM3D parameters
        set_dict["DYN_FLAG"] = 0			# Integration with dynamic code
        set_dict["DYN_SKIP"] = 0			# Number of dynamic events to skip (warm-up cycles)
        set_dict["DYN_M"] = 1e18			# Target seismic moment of dynamic event

        set_dict["NPROC"] = 1				# Number of processors, default 1 = serial (no MPI)


        self.set_dict = set_dict

        pass

    # Receive and store settings dict given by script file
    def settings(self, set_dict):
        for key in set_dict:
            val = set_dict[key]
            if type(val) == dict:
                self.set_dict[key].update(val)
            else:
                self.set_dict[key] = val
        return True

    # Render the mesh, and populate it uniformly with the input parameters
    # After render_mesh() was called, the self.mesh_dict can be modified
    # in the script file to introduce heterogeneity in the parameters
    # See examples/PyQDYN/Tse+Rice/tse_rice.py for an example
    def render_mesh(self):

        settings = self.set_dict
        dim = settings["MESHDIM"]
        N = settings["N"]
        self.set_dict["NW"] = N

        # Determine number of creep mechanisms
        A = np.array(settings["SET_DICT_CNS"]["A"])
        N_creep = len(A)
        self.set_dict["SET_DICT_CNS"]["N_CREEP"] = N_creep

        if dim == 0:
            N = 1
            self.set_dict["N"] = N
            self.set_dict["NW"] = N

        if dim == 2:
            N = self.set_dict["NX"]*self.set_dict["NW"]
            self.set_dict["N"] = N

        if N < 1:
            print("Number of mesh elements needs to be set before rendering mesh, unless MESHDIM = 0 (spring-block)")
            exit()

        mesh_params_general = ("IOT", "IASP", "V_PL", "SIGMA")
        mesh_params_RSF = ("V_0", "TH_0", "A", "B", "DC", "V1", "V2", "MU_SS", "V_SS", "CO")
        mesh_params_CNS = ("TAU", "PHI_INI", "A_TILDE", "MU_TILDE_STAR", "Y_GR_STAR", "H", "PHI_C", "PHI0", "THICKNESS")
        mesh_params_creep = ("A", "N", "M")
        mesh_params_localisation = ("LOCALISATION", "PHI_INI_BULK")
        mesh_params_TP = ("RHOC", "BETA", "ETA", "HALFW", "K_T", "K_P", "LAM", "P_A", "T_A", "DILAT_FACTOR")

        mesh_dict = {}

        # Populate mesh with general settings
        for param in mesh_params_general:
            mesh_dict[param] = np.ones(N)*settings[param]

        # Check for friction model
        if settings["FRICTION_MODEL"] == "RSF":
            for param in mesh_params_RSF:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_RSF"][param]
        elif settings["FRICTION_MODEL"] == "CNS":
            for param in mesh_params_CNS:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_CNS"][param]
            for i in range(N_creep):
                for param in mesh_params_creep:
                    param_no = "%s%i" % (param, i)
                    mesh_dict[param_no] = np.ones(N)*settings["SET_DICT_CNS"][param][i]
        else:
            print("FRICTION_MODEL '%s' not recognised. Supported models: RSF, CNS" % (settings["FRICTION_MODEL"]))
            exit()

        # Populate mesh with localisation parameters
        if settings["FEAT_LOCALISATION"] == 1:
            for param in mesh_params_localisation:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_LOCALISATION"][param]
            for i in range(N_creep):
                for param in mesh_params_creep:
                    param_bulk = "%s%i_bulk" % (param, i)
                    mesh_dict[param_bulk] = np.ones(N)*settings["SET_DICT_LOCALISATION"][param][i]

        # Populate mesh with thermal pressurisation parameters
        if settings["FEAT_TP"] == 1:
            for param in mesh_params_TP:
                mesh_dict[param] = np.ones(N)*settings["SET_DICT_TP"][param]

        # Mesh XYZ coordinates (and dip angle)
        mesh_dict["X"] = np.zeros(N)
        mesh_dict["Y"] = np.zeros(N)
        mesh_dict["Z"] = np.zeros(N)
        mesh_dict["DIP_W"] = np.zeros(N)

        if dim == 1:
            mesh_dict["X"] = np.linspace(-N/2, N/2, N)*settings["L"]/N
            mesh_dict["Y"] = np.ones(N)*settings["W"]

        if dim == 2:
            print("Warning: 2D faults currently require constant dip angle")
            theta = settings["DIP_W"]*np.pi/180.0
            mesh_dict["DW"] = np.ones(settings["NW"])*settings["DW"]
            mesh_dict["DIP_W"] = np.ones(N)*settings["DIP_W"]
            dx = 1.0*settings["L"]/settings["NX"]
            mesh_dict["X"] = (np.arange(N)%settings["NX"] + 0.5)*dx
            mesh_dict["Y"] = (np.floor(np.arange(N)/settings["NX"]) + 0.5)*settings["DW"]*np.cos(theta)
            mesh_dict["Z"] = mesh_dict["Y"]*np.tan(theta)

        self.mesh_dict.update(mesh_dict)
        self.mesh_rendered = True

        return True

    def write_input(self):

        if self.mesh_rendered == False:
            print("The mesh has not yet been rendered!")
            exit()

        # Optionally, output can be created to an external working directory
        # This is still experimental...
        if self.work_dir != "":
            if not os.path.isdir(self.work_dir):
                print("Creating %s" % (self.work_dir))
                call(["mkdir", self.work_dir])
            print("Switching working directory to %s" % (self.work_dir))
            os.chdir(self.work_dir)

        settings = self.set_dict
        mesh = self.mesh_dict
        Nprocs = settings["NPROC"]
        delimiter = "    "

        # Define chunk size for each processor
        nwLocal = (settings["NW"]//Nprocs)*np.ones(Nprocs, dtype=int)
        nwLocal[Nprocs-1] += settings["NW"] % Nprocs
        nnLocal = 0

        # Number of creep mechanisms
        N_creep = settings["SET_DICT_CNS"]["N_CREEP"]

        # Loop over processor nodes
        for iproc in range(Nprocs):

            nloc = nwLocal[iproc]*settings["NX"]	# number of fault elements hosted on node
            iloc = np.zeros(nloc, dtype=int)		# buffer for indices of those elements

            if Nprocs == 1: filename = "qdyn.in"
            else: filename = "qdyn%06i.in" % (iproc)

            # Start building the contents of our QDYN input file (input_str)
            input_str = ""
            input_str += "%u%s meshdim\n" % (settings["MESHDIM"], delimiter)

            # Input specific to 3D faults
            if settings["MESHDIM"] == 2:
                input_str += "%u %u%s NX, NW\n" % (settings["NX"], nwLocal[iproc], delimiter)
                input_str += "%.15g %.15g %.15g%s L, W, Z_CORNER\n" % (settings["L"], settings["W"], settings["Z_CORNER"], delimiter)
                # TODO: MPI the stuff below
                for i in range(settings["NW"]):
                    input_str += "%.15g %.15g \n" % (mesh["DW"][i], mesh["DIP_W"][i])
            else:
                input_str += "%u%s NN\n" % (settings["N"], delimiter)
                input_str += "%.15g %.15g%s L, W\n" % (settings["L"], settings["W"], delimiter)

            # Input specific to 2D faults
            if settings["MESHDIM"] == 1:
                input_str += "%u%s finite\n" % (settings["FINITE"], delimiter)

            # Input specific to fault rheological model (CNS or RSF)
            if settings["FRICTION_MODEL"] == "RSF":
                input_str += "%u%s itheta_law\n" % (settings["SET_DICT_RSF"]["THETA_LAW"], delimiter)
                input_str += "%u%s i_rns_law\n" % (settings["SET_DICT_RSF"]["RNS_LAW"], delimiter)
            else: # CNS model
                input_str += "%u%s itheta_law\n" % (1, delimiter)
                input_str += "%u%s i_rns_law\n" % (3, delimiter)

            # Some more general settings
            # Note that i_sigma_law is replaced by the feature stress_coupling, but is kept for compatibility
            input_str += "%u%s i_sigma_law\n" % (settings["SIGMA_CPL"], delimiter)
            # If the CNS model is used, define the number of creep mechanisms
            if settings["FRICTION_MODEL"] == "CNS":
                # Raise an error when the default value has not been overwritten
                if N_creep == -1:
                    print("The number of creep mechanisms was not set properly!")
                    print("This value should be determined in render_mesh()")
                    exit()
                input_str += "%u%s N_creep\n" % (N_creep, delimiter)
            input_str += "%u %u %u%s stress_coupling, thermal press., localisation\n" % (settings["FEAT_STRESS_COUPL"], settings["FEAT_TP"], settings["FEAT_LOCALISATION"], delimiter)
            input_str += "%u %u %u %u %u %u %u%s ntout_ot, ntout_ox, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n" % (settings["NTOUT_OT"], settings["NTOUT"], settings["IC"]+1, settings["NXOUT"], settings["NXOUT_DYN"], settings["OX_SEQ"], settings["OX_DYN"], delimiter)
            input_str += "%.15g %.15g %.15g %.15g %.15g %.15g%s beta, smu, lambda, v_th\n" % (settings["VS"], settings["MU"], settings["LAM"], settings["D"], settings["HD"], settings["V_TH"], delimiter)
            input_str += "%.15g %.15g%s Tper, Aper\n" % (settings["TPER"], settings["APER"], delimiter)
            input_str += "%.15g %.15g %.15g %.15g%s dt_try, dtmax, tmax, accuracy\n" % (settings["DTTRY"] ,settings["DTMAX"] ,settings["TMAX"] ,settings["ACC"] , delimiter)
            input_str += "%u%s nstop\n" % (settings["NSTOP"], delimiter)
            input_str += "%u %u%s DYN_FLAG, DYN_SKIP\n" % (settings["DYN_FLAG"], settings["DYN_SKIP"], delimiter)
            input_str += "%.15g %.15g %.15g%s M0, DYN_th_on, DYN_th_off\n" % (settings["DYN_M"], settings["DYN_TH_ON"], settings["DYN_TH_OFF"], delimiter)
            input_str += "%i %i%s FAULT_TYPE, SOLVER\n" % (settings["FAULT_TYPE"], settings["SOLVER"], delimiter)

            # Loop over all fault segments that are hosted on this processor node
            for i in range(nwLocal[iproc]):
                for j in range(settings["NX"]):
                    n = i*settings["NX"] + j
                    iloc[n] = j + i*settings["NX"] + nnLocal

            # Check for fault rheology model
            if settings["FRICTION_MODEL"] == "CNS":
                for i in range(nloc):
                    # Basic parameters (initial conditions, granular flow, porosity)
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g " % (mesh["SIGMA"][iloc[i]], mesh["TAU"][iloc[i]], mesh["PHI_INI"][iloc[i]], mesh["V_PL"][iloc[i]], mesh["A_TILDE"][iloc[i]], mesh["MU_TILDE_STAR"][iloc[i]], mesh["Y_GR_STAR"][iloc[i]], mesh["H"][iloc[i]], mesh["PHI_C"][iloc[i]], mesh["PHI0"][iloc[i]])
                    # Creep parameters (3 for each creep mechanism)
                    for j in range(N_creep):
                        input_str += "%.15g %.15g %.15g " % (mesh["A%i" % j][iloc[i]], mesh["N%i" % j][iloc[i]], mesh["M%i" % j][iloc[i]])
                    # Last couple of parameters
                    input_str += "%.15g %.15g %.15g\n" % (mesh["THICKNESS"][iloc[i]], mesh["IOT"][iloc[i]], mesh["IASP"][iloc[i]])
            else: # RSF model
                for i in range(nloc):
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["SIGMA"][iloc[i]], mesh["V_0"][iloc[i]], mesh["TH_0"][iloc[i]], mesh["A"][iloc[i]], mesh["B"][iloc[i]], mesh["DC"][iloc[i]], mesh["V1"][iloc[i]], mesh["V2"][iloc[i]], mesh["MU_SS"][iloc[i]], mesh["V_SS"][iloc[i]], mesh["IOT"][iloc[i]], mesh["IASP"][iloc[i]], mesh["CO"][iloc[i]], mesh["V_PL"][iloc[i]])

            # Check if localisation is requested
            if settings["FEAT_LOCALISATION"] == 1:
                if settings["FRICTION_MODEL"] != "CNS":
                    print("Localisation is compatible only with the CNS friction model")
                    exit()
                for i in range(nloc):
                    # Basic parameters (degree of localisation, initial bulk porosity)
                    input_str += "%.15g %.15g " % (mesh["LOCALISATION"][iloc[i]], mesh["PHI_INI_BULK"][iloc[i]])
                    # Creep parameters (3 for each creep mechanism)
                    for j in range(N_creep):
                        input_str += "%.15g %.15g %.15g " % (mesh["A%i_bulk" % j][iloc[i]], mesh["N%i_bulk" % j][iloc[i]], mesh["M%i_bulk" % j][iloc[i]])
                    # End the line
                    input_str += "\n"

            # Check if thermal pressurisation is requested
            if settings["FEAT_TP"] == 1:
                for i in range(nloc):
                    input_str += "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n" % (mesh["RHOC"][iloc[i]], mesh["BETA"][iloc[i]], mesh["ETA"][iloc[i]], mesh["HALFW"][iloc[i]], mesh["K_T"][iloc[i]], mesh["K_P"][iloc[i]], mesh["LAM"][iloc[i]], mesh["P_A"][iloc[i]], mesh["T_A"][iloc[i]], mesh["DILAT_FACTOR"][iloc[i]])

            # Add mesh grid location information
            for i in range(nloc):
                input_str += "%.15g %.15g %.15g %.15g\n" % (mesh["X"][iloc[i]], mesh["Y"][iloc[i]], mesh["Z"][iloc[i]], mesh["DIP_W"][iloc[i]])

            nnLocal += settings["NX"]*nwLocal[iproc]

            # Write QDYN input file
            with open(filename, "w") as f:
                f.write(input_str)

        self.qdyn_input_written = True

        return True

    # Run QDYN
    def run(self, test=False, unit=False):

        if not self.qdyn_input_written and unit is False:
            print("Input file has not yet been written. First call write_input() to render QDyn input file")
            exit()

        # Executable file (including path)
        qdyn_exec = os.path.join(self.qdyn_path, "qdyn")

        # If we're on a Windows 10 + Unix subsystem, we need to adjust executable path
        if self.W10_bash is True:
            # Replace Win-style backward slashes with Unix forward slashes
            qdyn_exec_bash = qdyn_exec.replace("\\", "/")
            # Find drive letter (usually C:, but not always)
            drive = re.match(r"([A-Z](?=:))", qdyn_exec_bash).group()
            # Replace drive with mounted partition
            qdyn_exec_bash = qdyn_exec_bash.replace("%s:" % (drive), "/mnt/%s" % (drive.lower()))

        # If serial (NPROC = 1)
        if self.set_dict["NPROC"] == 1:

            # If we're on a Windows 10 + Unix subsystem, call bash
            if self.W10_bash is True: cmd = ["bash", "-c", "\"%s\"" % (qdyn_exec_bash)]
            # If we're on Unix, simply call the qdyn executable directly
            else: cmd = [qdyn_exec]

            # If unit testing is requested, append test argument
            if unit:
                cmd.append("test")
            if not test:
                # Run command
                return call(cmd)
            else:
                # Run and suppress stdout
                with open(os.devnull, "w") as output:
                    call(cmd, stdout=output)

        else: # MPI parallel

            # If we're on a Windows 10 + Unix subsystem, call bash
            if self.W10_bash is True:
                cmd = ["bash", "-c",
                "\"/usr/local/bin/mpirun -np %i %s\"" % (self.set_dict["NPROC"], qdyn_exec_bash)]
            # If we're on Unix, simply call the mpi executable directly
            else:
                cmd = ["/usr/local/bin/mpirun", "-np", "%i" % (self.set_dict["NPROC"]), qdyn_exec]

            # If unit testing is requested, append test argument
            if unit:
                cmd.append("test")
            # Run command
            if not test:
                # Run command
                call(cmd)
            else:
                # Run and suppress stdout
                with open(os.devnull, "w") as output:
                    call(cmd, stdout=output)

        # If a suffix is requested, rename output files
        suffix = self.set_dict["SUFFIX"]
        if suffix != "":
            for filename in os.listdir("."):
                if filename.startswith("fort") and "_" not in filename: os.rename(filename, "%s_%s" % (filename, suffix))

        return True


    # Read QDYN output data
    def read_output(self, mirror=False, read_ot=True, filename_ot="fort.18", read_ox=True, filename_ox="fort.19"):

        # If a suffix is specified (optional), change file names
        suffix = self.set_dict["SUFFIX"]
        if suffix != "":
            filename_ot = "%s_%s" % (filename_ot, suffix)
            filename_ox = "%s_%s" % (filename_ox, suffix)

        # Default number of header lines
        nheaders = 6

        # Output file contents depends on the requested features
        if self.set_dict["FEAT_LOCALISATION"] == 1:
            #cols = ("t", "loc_size", "crack_size", "potcy", "pot_rate", "v", "theta", "theta2", "v_theta_dc", "tau", "x", "x_max", "v_max", "theta_max", "theta2_max", "omeg_max", "tau_max", "x_max", "sigma_max")
            cols = ("t", "loc_size", "crack_size", "potcy", "pot_rate", "v", "theta", "v_theta_dc", "tau", "slip", "x_max", "v_max", "theta_max", "omeg_max", "tau_max", "slip_max", "sigma_max")
        elif self.set_dict["FEAT_TP"] == 1:
            cols = ("t", "loc_size", "crack_size", "potcy", "pot_rate", "v", "theta", "v_theta_dc", "tau", "slip", "x_max", "v_max", "theta_max", "omeg_max", "tau_max", "slip_max", "sigma_max", "P", "T", "P_max", "T_max")
            nheaders = 8
        else:
            cols = ("t", "loc_size", "crack_size", "potcy", "pot_rate", "v", "theta", "v_theta_dc", "tau", "slip", "x_max", "v_max", "theta_max", "omeg_max", "tau_max", "slip_max", "sigma_max")

        # If time series data is requested
        if read_ot:
            data = read_csv(filename_ot, header=nheaders, names=cols, delim_whitespace=True)
            self.ot = data	# Store data in self.ot

            # Check if time series output was requested for other fault elements
            N_iot = int(np.sum(self.mesh_dict["IOT"]))
            self.N_iot = N_iot
            if N_iot > 0:
                self.iot = [self.read_other_output(i) for i in range(N_iot)]
        else:
            self.ot = None
            self.N_iot = 0

        if read_ox:

            # Snapshot output template depends on requested features
            if self.set_dict["FEAT_TP"] == 1:
                cols_ox = ("x", "t", "v", "theta", "dtau", "tau_dot", "slip", "sigma", "P", "T")
            else:
                cols_ox = ("x", "t", "v", "theta", "dtau", "tau_dot", "slip", "sigma")

            # Read snapshot output
            data_ox = read_csv(filename_ox, header=None, names=cols_ox, delim_whitespace=True, comment="#")

            # Store snapshot data in self.ox
            self.ox = data_ox

            # Sanitise output (check for near-infinite numbers, etc.)
            self.ox = self.ox.apply(pd.to_numeric, errors="coerce")

            # If free surface was generated manually (i.e. without FINITE = 2 or 3),
            # take only half data set (symmetric around first element)
            if mirror == True:
                data_ox = data_ox.loc[(data_ox["x"] > 0)]
                self.ox = data_ox
        else:
            self.ox = None

        return True

    # Read other time series data, when they are available
    def read_other_output(self, i):
        cols = ("t", "v", "theta", "tau", "x", "sigma")
        file_i = 10000+i+1
        filename = "fort.%i" % file_i
        suffix = self.set_dict["SUFFIX"]
        if suffix != "":
            filename = "%s_%s" % (filename, suffix)
        data = read_csv(filename, header=2, names=cols, delim_whitespace=True)
        return data

    # Optionally, export dicts with all of the simulation settings. It is good
    # practice to call this function before running each simulation, so that
    # the simulation parameters are permanently stored (i.e. independent of any
    # changes made to the script file)
    def export_dicts(self):
        print("Exporting settings/mesh dict...")
        with gzip.GzipFile("set_dict.tar.gz", "w") as f:
            pickle.dump(self.set_dict, f, pickle.HIGHEST_PROTOCOL)
        with gzip.GzipFile("mesh_dict.tar.gz", "w") as f:
            pickle.dump(self.mesh_dict, f, pickle.HIGHEST_PROTOCOL)






