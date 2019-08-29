# Importing some required modules
import os
import sys
from termcolor import colored
import warnings
warnings.filterwarnings("ignore")
# Go up in the directory tree
upup = [os.pardir]*2
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
from vstep import TestVstep
from stickslip import TestStickSlip
from singleasperity import TestSingleAsperity
from tserice import TestTseRice

# Width of output
msg_width = 66

# General simulation parameters
set_dict = {
    "FAULT_TYPE": 1,
    "ACC": 1e-10,
    "SOLVER": 2,
    "MU": 3e10,
    "TMAX": 20,
    "DTTRY": 1e-6,
    "MESHDIM": 0,
    "NTOUT": 10000,
    "VS": 3000,
    "SIGMA": 5e6,
    "L": 1,
    "W": 1,
    "FEAT_STRESS_COUPL": 0,
    "FEAT_TP": 0,
    "FEAT_LOCALISATION": 0,
    "D": 0,
    "HD": 0,
}

# RSF parameters
set_dict_RSF = {
    "RNS_LAW": 0,
    "THETA_LAW": 1,
    "A": 0.001,
    "B": 0.0015,
    "DC": 1e-5,
    "V1": 1e-2,
    "V2": 1e-7,
    "MU_SS": 0.6,
    "V_SS": 1e-6,
    "CO": 0,
}

# CNS parameters
set_dict_CNS = {
    "THICKNESS": 1e-5,
    "A_TILDE": 0.01,
    "MU_TILDE_STAR": 0.4,
    "Y_GR_STAR": 1e-6,
    "H": 0.577,
    "PHI_C": 0.32,
    "PHI0": 0.03,
    # Creep kinetics (note compatibility patch)
    "IPS_CONST_DIFF": 0.0,
    "IPS_CONST_DISS": 1e-10,
    "A": [1e-10],
    "N": [1],
    "M": [1],
}

# Add RSF/CNS dictionary to QDYN dictionary
set_dict["SET_DICT_RSF"] = set_dict_RSF
set_dict["SET_DICT_CNS"] = set_dict_CNS

qdyn_files = ["qdyn.in", "fort.18", "fort.19", "fort.22", "fort.121"]

# QDYN class object
p = qdyn()

# Test report header
print("".join(["="]*msg_width))
print("".join([" "]*((msg_width-25)//2)) + "Initiated QDYN Test Suite\n")

# Unit tests (hard coded in FORTRAN)
print(" Running unit tests...")
unit_failed = p.run(test=False, unit=True)

if unit_failed:
    print("Unit tests %s" % colored("FAILED", "red"))
    exit(1)

# Integration tests (executed from Python)
print(" Running integration tests...")

# Spring-block velocity-step simulation (RSF and CNS)
print(" - Testing spring-block (velocity-step)..")
p.settings(set_dict)
vstep = TestVstep(p)
vstep.import_results()
vstep.run_test("RSF")
vstep.run_test("CNS")

# Spring-block stick-slip simulation (RSF and CNS)
print(" - Testing spring-block (stick-slip)..")
p.settings(set_dict)
stickslip = TestStickSlip(p)
stickslip.import_results()
stickslip.run_test("RSF")
stickslip.run_test("CNS")

# 2D fault single asperity simulation (RSF)
print(" - Testing single asperity (will take a few minutes)...")
p.settings(set_dict)
single_asperity = TestSingleAsperity(p)
single_asperity.import_results()
single_asperity.run_test("RSF")
single_asperity.run_test("CNS")
# Plot results
# single_asperity.plot_results("RSF")
# single_asperity.plot_results("CNS")

# Tse & Rice (1986) example test (RSF)
# see https://doi.org/10.1029/JB091iB09p09452
print(" - Testing Tse & Rice (1986) example (will take a few minutes)...")
p.settings(set_dict)
tse_rice = TestTseRice(p)
tse_rice.import_results()
tse_rice.run_test()
# Plot results
# tse_rice.plot_results("RSF")

# Print out integration test report
print("".join(["-"]*msg_width))
print("Integration test results:\n")
print(" - Spring-block (velocity-step)")
print("     %s" % vstep.test_results["RSF"]["success_msg"])
print("     %s" % vstep.test_results["CNS"]["success_msg"])
print(" - Spring-block (stick-slip)")
print("     %s" % stickslip.test_results["RSF"]["success_msg"])
print("     %s" % stickslip.test_results["CNS"]["success_msg"])
print(" - Single asperity")
print("     %s" % single_asperity.test_results["RSF"]["success_msg"])
print("     %s" % single_asperity.test_results["CNS"]["success_msg"])
print(" - Tse & Rice (1986) example")
print("     %s" % tse_rice.test_results["RSF"]["success_msg"])
print("".join(["="]*msg_width))

# Clean-up
for file in qdyn_files:
    os.remove(file)
