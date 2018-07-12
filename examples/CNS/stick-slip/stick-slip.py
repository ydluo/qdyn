# Path where QDYN executable and wrapper are located
qdyn_path = "/home/martijn/QDyn/src"

# Importing some required modules
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(qdyn_path)
from pyqdyn import qdyn

# QDYN class object
p = qdyn()

# Define where the QDYN executable is located
p.qdyn_path = qdyn_path

# Python dictionary with general settings
set_dict = {
    "FRICTION_MODEL": "CNS",
    "ACC": 1e-10,
    "SOLVER": 2,
    "MU": 2e10,
    "TMAX": 1e4,
    "DTTRY": 1e-6,
    "MESHDIM": 0,
    "NTOUT": 10000,
    "SIGMA": 5e6,
    "V_PL": 1e-6,
    "L": 1e1,
}

# Python dictionary with CNS parameters
set_dict_CNS = {
    "IPS_CONST_DISS": 1e-10,
    "A_TILDE": 0.02,
    "MU_TILDE_STAR": 0.4,
    "Y_GR_STAR": 1e-6,
    "PHI_INI": 0.25,
    "THICKNESS": 1e-4,
    "TAU": 3e6,
}

# Add CNS dictionary to QDYN dictionary
set_dict["SET_DICT_CNS"] = set_dict_CNS

p.settings(set_dict)
p.render_mesh()

# Write input file
p.write_input()

# Run simulation
p.run()

# Get our results
p.read_output()

# Plot results
plt.subplot(211)
plt.plot(p.ot["t"], p.ot["tau"] / p.set_dict["SIGMA"])
plt.ylabel("friction [-]")
plt.subplot(212)
plt.plot(p.ot["t"], p.ot["theta"] * 100)
plt.xlabel("time [s]")
plt.tight_layout()
plt.show()
