import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn
from bisect import bisect_left
from datetime import datetime

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

# Instantiate QDYN wrapper
p = qdyn()
# Get current date
now = datetime.now()
# Define output directory; create if needed
data_dir = "BP_output"
if not os.path.isdir(data_dir):
    os.mkdir(data_dir)

# Output header files
static_header = ""
static_header += "# problem=SEAS Benchmark No.1\n"
static_header += "# code=QDYN\n"
static_header += "# version=2.0\n"
static_header += "# modeler=Luo et al.\n"
static_header += "# date=%s\n" % now.strftime("%Y/%m/%d")


def write_timeseries(data, IOT, dz):
    """
    Helper routine to write time-series data to file
    """

    # If data size exceeds limit: downsample
    Nt = len(data[0])
    stride = 1
    if Nt > int(1e5):
        stride = 10
        Nt = Nt // stride

    # Loop over all depth intervals
    for i, z in enumerate(IOT):
        # Generate headers
        header = "" + static_header
        header += "# element_size=%.3fm\n" % dz
        header += "# location= on fault, %.3fkm depth\n" % (z * 1e-3)
        header += "# num_time_steps=%d\n" % Nt
        header += "# Column #1 = Time (s)\n"
        header += "# Column #2 = Slip (m)\n"
        header += "# Column #3 = Slip rate (log10 m/s)\n"
        header += "# Column #4 = Shear stress (MPa)\n"
        header += "# Column #5 = State (log10 s)\n"
        header += "#\n"
        header += "t slip slip_rate shear_stress state\n"

        # Filename
        filename = "fltst_dp%03d_dz=%.1f" % (1e-2 * z, dz)
        file = os.path.join(data_dir, filename)

        # Write headers
        with open(file, "w") as f:
            f.write(header)

        # Extract output data (and downsample)
        output = data[i][["t", "x", "v", "tau", "theta"]][::stride]
        # Convert slip_rate, shear_stress, theta
        output["v"] = np.log10(output["v"])
        output["tau"] *= 1e-6
        output["theta"] = np.log10(output["theta"])
        # Write data to file (append)
        with open(file, "a") as f:
            output.to_csv(f, index=False, header=False, sep=" ")

    print("Time-series output successful")
    pass


def write_slip(ox, z, dz):
    """
    Helper routine to write snapshot data to file
    """

    mask = np.isfinite(ox["x"])
    x = ox["x"][mask].unique()
    x_order = np.argsort(x)
    t_vals = np.sort(ox["t"].unique())

    Nx = len(x)
    Nt = len(t_vals)
    t_vals = t_vals[:-1]

    slip = ox["slip"][:Nx * (Nt - 1)].values.reshape((Nt - 1, Nx))
    slip = slip.T[x_order].T
    v = ox["v"][:Nx * (Nt - 1)].values.reshape((Nt - 1, Nx))
    v = v.T[x_order].T
    v_max = np.array([np.nanmax(v[i]) for i in range(Nt - 1)])
    v_max = np.log10(v_max)

    # First row contains depth intervals
    first_row = np.hstack([0, 0, z])
    # Create output buffer
    output = np.hstack([t_vals.reshape(-1, 1), v_max.reshape(-1, 1), slip])
    output = np.vstack([first_row, output])
    output = pd.DataFrame(output)

    # Create headers
    header = "" + static_header
    header += "# element_size=%.3f\n" % dz
    header += "# Row #1 = Depth (m) with two zeros first\n"
    header += "# Column #1 = Time (s)\n"
    header += "# Column #2 = Max slip rate (log10 m/s)\n"
    header += "# Columns #3-%d = Slip (m)\n" % (slip.shape[1] + 2)
    header += "z\n"
    header += "t max_slip_rate slip\n"

    # Filename
    filename = "devol_dz=%.1f.dat" % dz
    file = os.path.join(data_dir, filename)

    # Write headers
    with open(file, "w") as f:
        f.write(header)

    # Write data to file (append)
    with open(file, "a") as f:
        output.to_csv(f, index=False, header=False, sep=" ")

    print("Slip output successful")
    pass


def slip_distribution_profile(ox, t_step, t_step_subseismic, t_step_seismic, slip_ref, depth=15e3):
    """
    Helper routine to plot the snapshot data (slip contours)
    """
    mask = np.isfinite(ox["x"])
    x = ox["x"][mask].unique()
    x_order = np.argsort(x)
    t_vals = np.sort(ox["t"].unique())

    z0 = 0
    z_max = depth*1e-3
    z = np.linspace(0, z_max, len(x)) + z0

    Nx = len(x)
    Nt = len(t_vals)
    t_vals = t_vals[:-1]
    slip = ox["slip"][:Nx*(Nt-1)].values.reshape((Nt-1, Nx))
    slip = slip.T[x_order].T
    v = ox["v"][:Nx*(Nt-1)].values.reshape((Nt-1, Nx))
    v = v.T[x_order].T

    v_max = np.array([np.nanmax(v[i]) for i in range(Nt-1)])

    v_subseismic = 1e-7
    v_seismic = 1e-3

    t_prev = 0
    inds_seismic = (v_max >= v_seismic)
    inds_subseismic = (v_max >= v_subseismic) & (v_max < v_seismic)

    ref_ind = np.where(slip[:,0] > slip_ref)[0][0]
    slip_ref = slip[ref_ind,:]

    fig = plt.figure(figsize=(15,8), facecolor="white")
    colours = seaborn.color_palette("deep", 5)
    colours[0] = "b"
    colours[1] = "r"
    colours[2] = "b"

    ax = fig.add_subplot(111)

    for i in range(Nt-1):
        if inds_seismic[i]:
            if t_vals[i] > t_prev + t_step_seismic:
                plt.plot(slip[i]-slip_ref, z, ls="--", c=colours[1], lw=0.8)
                t_prev = t_vals[i]
        elif inds_subseismic[i]:
            if t_vals[i] > t_prev + t_step_subseismic:
                plt.plot(slip[i]-slip_ref, z, ls="-", c=colours[2], lw=1.0)
                t_prev = t_vals[i]
        else:
            if t_vals[i] > t_prev + t_step:
                plt.plot(slip[i]-slip_ref, z, ls="-", c=colours[0], lw=1.5)
                t_prev = t_vals[i]

    t_day = 24*3600.0
    t_yr = 365*t_day
    plt.plot([np.nan]*2, [np.nan]*2, "-", c=colours[0], label="Interseismic (%.0f yr)" % (t_step/t_yr))
    plt.plot([np.nan]*2, [np.nan]*2, "-", c=colours[2], label="Subseismic (%.1f day)" % (t_step_subseismic/t_day))
    plt.plot([np.nan]*2, [np.nan]*2, "--", c=colours[1], label="Coseismic (%.1f sec)" % (t_step_seismic))
    plt.legend(bbox_to_anchor=(0.0, 1.1, 1.0, .102), loc="center", ncol=3, borderaxespad=0.0)

    plt.ylim((np.min(z), np.max(z)))
    plt.xlim((0, np.max(slip)-np.max(slip_ref)))
    plt.ylabel("depth [km]")
    plt.xlabel("accumulated slip [m]")
    plt.gca().invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()


""" Simulation parameters """

t_yr = 3600*24*365.0        # Seconds in 1 year [s]
t_max = 1200*t_yr           # Simulation time [s]

rho = 2670.0                # Rock density [kg/m3]
Vs = 3464.0                 # Shear wave speed [m/s]
mu = rho * Vs**2            # Shear modulus [Pa]
eta = mu / (2 * Vs)
sigma = 50e6                # Effective normal stress [Pa]
V_pl = 1e-9                 # Plate rate [m/s]
V_ini = V_pl                # Initial velocity [m/s]
H = 15e3                    # Depth extent of VW region [m]
h = 3e3                     # Width of VW-VS transition region [m]
L = 40e3                    # With of frictional fault [m]

# Target element sizes: {25, 50, 100, 200, 400, 800} [m]
dz = 25.0                   # Target cell size
N_approx = L / dz           # Approximate number of fault elements

# Compute the nearest power of 2
logN = np.log2(N_approx)
N = int(np.power(2, np.round(logN)))
dz2 = L / N
print("Target dz: %.1f \t Current dz: %.1f" % (dz, dz2))

dz_slip = 100.0
dN_slip = 1
if dz2 < dz_slip:
    logdN_slip = np.log2(dz_slip / dz2)
    dN_slip = int(np.power(2, np.round(logdN_slip)))

""" Rate-and-state friction parameters """

a0 = 0.010          # Min a
a_max = 0.025       # Max a
b = 0.015           # Constant b
Dc = 0.004          # Critical slip distance
f_ref = 0.6         # Reference friction value
V_ref = 1e-6        # Reference velocity

z = np.arange(dz2 / 2, L, dz2)      # Depth vector
a = a0 + (a_max - a0) * (z - H)/h   # Depth-dependent values of a
a[z < H] = a0                       # Shallow cutoff
a[z >= H + h] = a_max               # Deep cutoff

# Construct initial stress state
tau_exp = np.exp( (f_ref + b*np.log(V_ref / V_ini)) / a_max )
tau_sinh = np.arcsinh(V_ini * tau_exp / (2 * V_ref))
tau_ini = sigma * a_max * tau_sinh + eta * V_ini

# Construct initial "state" state
theta_sinh = np.sinh((tau_ini - eta * V_ini) / (a * sigma))
theta_log = np.log(theta_sinh * 2 * V_ref / V_ini)
theta_ini = (Dc / V_ref) * np.exp((a / b) * theta_log - f_ref / b)

# Requested timeseries output depths
z_IOT = np.array([0.0, 2.4, 4.8, 7.2, 9.6, 12.0, 14.4, 16.8, 19.2, 24.0, 28.8, 36.0]) * 1e3
# Corresponding fault element indices
IOT = [bisect_left(z, zi) for zi in z_IOT]
# Get actual output depth intervals
z_IOT_true = z[IOT]

# Python dictionary with general settings
set_dict = {
    "N": N,   					# number of fault segments
    "NXOUT": dN_slip,   	    # snapshot output spacing
    "NTOUT": 100,   			# snapshot output frequency
    "ACC": 1e-7,   				# Solver accuracy
    "MU": mu,   				# Shear modulus
    "DTTRY": 1e-8,   			# First time step (needs to be small)
    "TMAX": t_max,   		    # Run simulation 200 years
    "MESHDIM": 1,   			# One-dimensional fault
    "VS": Vs,   				# Shear wave velocity
    "L": L,        				# Fault length (depth)
    "FINITE": 3,   				# Finite fault with free surface
    "IC": 0,					# Location for time-series output (7.5 km depth)
    "SOLVER": 2,                # Runge-Kutta solver
    "V_PL": V_pl,               # Plate velocity
    "FAULT_TYPE": 1,            # Strike-slip fault (right-lateral)
}

# Python dictionary with rate-and-state friction parameters
set_dict_RSF = {
    "RNS_LAW": 2,		# Using regularised RSF
    "THETA_LAW": 1,		# Using the ageing law
    "DC": Dc,		    # Dc = 40 mm
    "V_0": V_ini,		# Initial velocity
    "V_SS": V_ref,	    # Reference velocity
    "MU_SS": f_ref,     # Reference friction
    "A": a0,            # Direct effect
    "B": b,             # Evolution effect
    "SIGMA": sigma,     # Effective normal stress
}

set_dict["TH_0"] = 0.99*set_dict_RSF["DC"]/set_dict_RSF["V_0"]
set_dict["SET_DICT_RSF"] = set_dict_RSF

p.settings(set_dict)
p.render_mesh()

p.mesh_dict["A"][:] = a[:]
p.mesh_dict["TH_0"][:] = theta_ini[:]
p.mesh_dict["IOT"][IOT] = 1

print("Writing input file")
p.write_input()
print("Running simulation")
p.run()
print("Reading output")
p.read_output(read_ot=True, read_ox=True)

print("Writing results")
write_timeseries(data=p.iot, IOT=z_IOT_true, dz=dz2)
write_slip(ox=p.ox, z=z[::set_dict["NXOUT"]], dz=dz2)

print("Plotting results")
slip_distribution_profile(
    ox=p.ox, t_step=5*t_yr, t_step_subseismic=5*t_yr, t_step_seismic=1.0,
    slip_ref=0.0, depth=set_dict["L"]
)
