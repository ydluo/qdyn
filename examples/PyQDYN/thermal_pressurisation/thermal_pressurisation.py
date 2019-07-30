# Usage:
#
# Run simulation with 'python thermal_pressurisation.py run'

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

# Load-point velocity [m/s]
V = 1e-6

# Python dictionary with general settings
set_dict = {
    "FRICTION_MODEL": "RSF",
    "ACC": 1e-8,
    "MU": 30e9,
    "DTTRY": 1e-6,
    "MESHDIM": 0,
    "NTOUT": 1000,
    "FEAT_TP": 1,
    "L": 300,
    "SIGMA": 5e6,
    "SOLVER": 2,
    "V_PL": V,
    "TMAX": 1e4,
}

set_dict_RSF = {
    "A": 1e-3,
    "B": 5e-3,
    "DC": 1e-5,
    "V_0": V,
    "V_SS": V,
    "TH_0": 1e-5 / V,
}

set_dict_TP = {
    "HALFW": 0.5e-3,
    "RHOC": 2.16e6,
    "BETA": 2e-9,
    "ETA": 2e-4,
    "LAM": 2.78e-4,
    "K_T": 2.0,
    "K_P": 1e-16,
    "P_A": 1e6,
    "T_A": 273.0,
}

set_dict["SET_DICT_RSF"] = set_dict_RSF
set_dict["SET_DICT_TP"] = set_dict_TP

# Feed our settings to QDYN and render mesh
p.settings(set_dict)
p.render_mesh()
p.write_input()

# Run the simulation
p.run()

# Read model output
p.read_output()

# Plot some output
plt.subplot(311)
plt.plot(p.ot["t"], p.ot["tau"]/5e6)
plt.ylabel("tau / sigma [-]")
plt.twinx()
plt.plot(p.ot["t"], p.ot["P"]*1e-6-1, "g-")
plt.ylabel("dP [MPa]")

plt.subplot(312)
plt.plot(p.ot["t"], p.ot["T"]-273.0, "r")
plt.ylabel("dT [K]")

plt.subplot(313)
plt.plot(p.ot["t"], p.ot["v"])
plt.ylabel("V [m/s]")
plt.xlabel("time [s]")
plt.yscale("log")

plt.tight_layout()
plt.show()
