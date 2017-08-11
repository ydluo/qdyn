# Path where QDYN executable and wrapper are located
qdyn_path = "/home/martijn/QDyn/src"

# Importing some required modules
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(qdyn_path)
from pyqdyn import qdyn

# CNS kinetic parameters

d = 10e-6			# Nominal grain size
O = 2.69e-5			# Molar volume
R = 8.3144			# Universal gas constant
T = 293				# Absolute temperature
DCS0 = 2.79e-15		# Kinetic pre-factor
dH = 2.45e4			# Activation energy
DCS = DCS0*np.exp(-dH/(R*T))			# Kinetic parameter at T
Z_ps = (24/np.pi)*DCS*O/(R*T*d**3)		# Combined (lumped) kinetic constant

# QDYN class object
p = qdyn()

# Define where the QDYN executable is located
p.qdyn_path = qdyn_path

# Python dictionary with general settings
set_dict = {
	"FRICTION_MODEL": "CNS",
	"ACC": 1e-10,
	"MU": 2e10,
	"TMAX": 20,
	"DTTRY": 1e-6,
	"MESHDIM": 0,
	"NTOUT": 10000,
	"VS": 0,
}

# Python dictionary with CNS parameters
set_dict_CNS = {
	"IPS_CONST_DIFF": Z_ps,
	"A_TILDE": 0.02,
	"MU_TILDE_STAR": 0.4,
	"Y_GR_STAR": 1e-6,
	"PHI_INI": 0.25,
	"THICKNESS": 1e-4,
	"TAU": 3e6,	
	"SIGMA": 5e6,
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

# Loop over all velocity steps
for i, V in enumerate(Vs):
	# Set load-point velocity
	set_dict["SET_DICT_CNS"]["V_LP"] = V
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
	plt.plot(p.ot["t"]-p.ot["t"].iloc[0]+t_final, p.ot["tau"]/p.set_dict["SET_DICT_CNS"]["SIGMA"])
	plt.subplot(212)
	plt.plot(p.ot["t"]-p.ot["t"].iloc[0]+t_final, p.ot["theta"]*100)

	# Set starting point for next simulation
	tau_final = p.ot["tau"].values[-1]
	phi_final = p.ot["theta"].values[-1]
	t_final += p.ot["t"].values[-1]

plt.subplot(211)
plt.ylabel("friction [-]")
plt.subplot(212)
plt.ylabel("porosity [%]")
plt.xlabel("time [s]")
plt.tight_layout()
plt.show()

