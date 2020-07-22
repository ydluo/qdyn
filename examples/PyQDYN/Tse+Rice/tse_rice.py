# Usage:
# 
# Run simulation with 'python tse_rice.py run'
# Omit the run argument to execute script file without
# re-running the simulation (e.g. for plotting data)

# Importing some required modules
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# Go up in the directory tree
upup = [os.pardir]*4
qdyn_dir = os.path.join(*upup)
# Get QDYN src directory
src_dir = os.path.abspath(
    os.path.join(
        os.path.join(__file__, qdyn_dir), "src")
)
# Append src directory to Python path
sys.path.append(src_dir)
# Import QDYN wrapper
from pyqdyn import qdyn

# QDYN class object
p = qdyn()

# Number of fault segments
N = np.power(2, 11)

# Number of secs in year
t_yr = 3600*24*365.0

# Plate velocity = 35 mm/yr
V = 35e-3 / t_yr

# Python dictionary with general settings
set_dict = {
    "N": N,   					# number of fault segments
    "NXOUT": np.power(2,1),   	# snapshot output spacing
    "NTOUT": 100,   			# snapshot output frequency
    "ACC": 1e-7,   				# Solver accuracy
    "MU": 30e9,   				# Shear modulus
    "DTTRY": 1e-8,   			# First time step (needs to be small)
    "TMAX": 200*t_yr,   		# Run simulation 200 years
    "MESHDIM": 1,   			# One-dimensional fault
    "VS": 3000,   				# Shear wave velocity
    "L": 30.0e3,   				# Fault length (depth)
    "FINITE": 3,   				# Finite fault with free surface
    "IC": N//4,					# Location for time-series output (7.5 km depth)
    "SOLVER": 1,
    "V_PL": V,
}

# Python dictionary with rate-and-state friction parameters
set_dict_RSF = {
    "RNS_LAW": 0,		# using classical rate-and-state
    "THETA_LAW": 1,		# using the ageing law
    "DC": 40e-3,		# Dc = 40 mm
    "V_0": V,			# Initial velocity
    "V_SS": V,			# Reference velocity
}

# Set state variable near steady-state value
set_dict_RSF["TH_0"] = 0.99*set_dict_RSF["DC"]/set_dict_RSF["V_0"]
set_dict["SET_DICT_RSF"] = set_dict_RSF

"""
Compute depth-dependent parameters
"""

# Polynomial coefficients for temperature data of Lachenbruch & Sass (1973)
T_poly_coeffs = [3.96094744e-06, -6.94566205e-04, 4.37987070e-02, -1.38543340e+00, 3.69151304e+01, 8.90082578e-01]
T_int = np.poly1d(T_poly_coeffs)

# Depth vector
z = np.linspace(0, 30e3, N)

# Temperature vector
T = T_int(z*1e-3)

# RSF parameters			
a = 3.28e-5*T - 9.288e-3
a_min_b = a
a = np.clip(a, a_min=0.004, a_max=10)
a_min_b = np.clip(a_min_b, a_min=-0.0029, a_max=10)
b = -(a_min_b - a)

# Effective normal stress
sigma = 18e3*z + 1e7

print("Number of fault elements: %i \t Element size: %.2f m" % (N, np.max(z)/N))

# Set-up simulation
def setup():

    # Feed our settings to QDYN and render mesh
    p.settings(set_dict)
    p.render_mesh()

    # Populate mesh dict with depth-dependent parameters
    p.mesh_dict["SIGMA"][:] = sigma[:]
    p.mesh_dict["A"][:] = a[:]
    p.mesh_dict["B"][:] = b[:]

    # Write input file
    p.write_input()

    # Run the simulation when the "run" argument is provided
    try:
     if sys.argv[1] == "run":
        # Export settings dicts (optional, but good practice)
        p.export_dicts()
        # Run simulation
        p.run()
    except IndexError:
        pass

    # Read model output
    p.read_output()


setup()

#Plot some output
plt.plot(p.ot_vmax["t"]/t_yr, p.ot_vmax["v"])
plt.yscale("log")
plt.ylabel(r"$V_{max}$ [m/s]")
plt.xlabel("time [yr]")
plt.tight_layout()
plt.show()