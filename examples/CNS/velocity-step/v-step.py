# Importing some required modules
import gzip
import pickle

import matplotlib.pyplot as plt
import numpy as np

# Import QDYN wrapper
from qdyn.pyqdyn import qdyn
from numpy.testing import assert_allclose

# QDYN class object
p = qdyn()

Z_ps = 1e-10

# Python dictionary with general settings
set_dict = {
    "FRICTION_MODEL": "CNS",
    "ACC": 1e-10,
    "SOLVER": 2,
    "MU": 2e10,
    "TMAX": 20,
    "DTTRY": 1e-6,
    "MESHDIM": 0,
    "NTOUT": 10000,
    "VS": 0,
    "SIGMA": 5e6,
}

# Python dictionary with CNS parameters
set_dict_CNS = {
    "A": [0.5*Z_ps, 0.5*Z_ps, 0.0*Z_ps],
    "N": [1, 1, 1],
    "M": [1, 1, 1],
    "A_TILDE": 0.02,
    "MU_TILDE_STAR": 0.4,
    "Y_GR_STAR": 1e-6,
    "PHI_INI": 0.25,
    "THICKNESS": 1e-4,
    "TAU": 3e6,
}

# Add CNS dictionary to QDYN dictionary
set_dict["SET_DICT_CNS"] = set_dict_CNS

# Velocity step sequence 10 > 15 > 10 micron/s
Vs = [10e-6, 15e-6, 10e-6]

# Define some starting points for the first simulation
tau_final = 3e6
phi_final = 0.25
t_final = 0

# Total slip distance per simulation
x_ss = 500e-6

t_all = np.array([])
tau_all = np.array([])
phi_all = np.array([])

# Loop over all velocity steps
for i, V in enumerate(Vs):
    # Set load-point velocity
    set_dict["V_PL"] = V
    # Set initial values from previous step
    set_dict["SET_DICT_CNS"]["PHI_INI"] = phi_final
    set_dict["SET_DICT_CNS"]["TAU"] = tau_final
    # Set simulated time
    set_dict["TMAX"] = x_ss/V

    # Feed our settings to QDYN and render mesh
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
    plt.plot(p.ot[0]["t"]-p.ot[0]["t"].iloc[0]+t_final, p.ot[0]["tau"]/p.set_dict["SIGMA"])
    plt.subplot(212)
    plt.plot(p.ot[0]["t"]-p.ot[0]["t"].iloc[0]+t_final, p.ot[0]["theta"]*100)

    t_all = np.hstack([t_all, p.ot[0]["t"]-p.ot[0]["t"].iloc[0]+t_final])
    tau_all = np.hstack([tau_all, p.ot[0]["tau"]])
    phi_all = np.hstack([phi_all, p.ot[0]["theta"]])

    # Set starting point for next simulation
    tau_final = p.ot[0]["tau"].values[-1]
    phi_final = p.ot[0]["theta"].values[-1]
    t_final += p.ot[0]["t"].values[-1]

plt.subplot(211)
plt.ylabel("friction [-]")
plt.subplot(212)
plt.ylabel("porosity [%]")
plt.xlabel("time [s]")
plt.tight_layout()
plt.show()

with gzip.GzipFile("benchmark.tar.gz", "r") as f:
    benchmark = pickle.load(f)

plt.figure(2)
plt.plot(t_all, tau_all*1e-6, label="Current")
plt.plot(benchmark["t"], benchmark["tau"]*1e-6, "k--", label="Benchmark")
plt.xlabel("time [s]")
plt.ylabel("shear stress [MPa]")
plt.legend(loc=4, ncol=2)
plt.tight_layout()
plt.show()

print("Performing benchmark comparison")
assert_allclose(benchmark["t"], t_all, rtol=1e-6)
assert_allclose(benchmark["tau"], tau_all, rtol=1e-6)
assert_allclose(benchmark["phi"], phi_all, rtol=1e-6)
print("Benchmark comparison OK")
