# Importing some required modules
import os

import numpy as np
from aux import AuxiliaryFunctions


class TestTseRice(AuxiliaryFunctions):

    test_results = {}
    dirname = "tserice"
    # The git hash below is an identifier for the version of QDYN that was
    # used to generate the benchmark results. This hash is included in the
    # frozen benchmark results, and it is checked against each time the results
    # are imported. When a new benchmark is generated, this hash should be
    # updated.
    frozen_hash = "9c0281dda0e913633fa7f3c1c26422c0bec38db0"
    frozen_loaded = False

    def __init__(self, p):
        AuxiliaryFunctions.__init__(self)
        self.p = p
        pickle_name = "tserice_results_frozen_%s.tar.gz" % self.frozen_hash[:10]
        self.pickle_file = os.path.join(self.dirname, pickle_name)
        pass

    def run_test(self):

        p = self.p
        set_dict = p.set_dict
        t_yr = 3600 * 24 * 365.0
        N = np.power(2, 11)

        # Compute depth-dependent parameters

        # Polynomial coefficients for temperature data of Lachenbruch & Sass (1973)
        T_poly_coeffs = [3.96094744e-06, -6.94566205e-04, 4.37987070e-02,
                         -1.38543340e+00, 3.69151304e+01, 8.90082578e-01]
        T_int = np.poly1d(T_poly_coeffs)

        # Depth vector
        z = np.linspace(0, 30e3, N)

        # Temperature vector
        T = T_int(z * 1e-3)

        # RSF parameters
        a = 3.28e-5 * T - 9.288e-3
        a_min_b = a
        a = np.clip(a, a_min=0.004, a_max=10)
        a_min_b = np.clip(a_min_b, a_min=-0.0029, a_max=10)
        b = -(a_min_b - a)

        # Effective normal stress
        sigma = 18e3 * z + 1e7

        set_dict["N"] = N
        set_dict["L"] = 30e3
        set_dict["MESHDIM"] = 1
        set_dict["FINITE"] = 3
        set_dict["TMAX"] = 500*t_yr
        set_dict["NTOUT"] = 1000
        set_dict["NXOUT"] = np.power(2, 3)
        set_dict["V_PL"] = 35e-3 / t_yr
        set_dict["MU"] = 3e10
        set_dict["SIGMA"] = 1e8
        set_dict["ACC"] = 1e-7
        set_dict["DTTRY"] = 100
        set_dict["V_TH"] = 1e-2
        set_dict["SOLVER"] = 1
        set_dict["FRICTION_MODEL"] = "RSF"

        # Setting some RSF parameters
        set_dict["SET_DICT_RSF"]["RNS_LAW"] = 0
        set_dict["SET_DICT_RSF"]["THETA_LAW"] = 1
        set_dict["SET_DICT_RSF"]["DC"] = 40e-3
        set_dict["SET_DICT_RSF"]["V_SS"] = set_dict["V_PL"]
        set_dict["SET_DICT_RSF"]["V_0"] = set_dict["V_PL"]
        set_dict["SET_DICT_RSF"]["TH_0"] = \
            0.99*set_dict["SET_DICT_RSF"]["DC"] / set_dict["V_PL"]

        p.settings(set_dict)
        p.render_mesh()

        p.mesh_dict["SIGMA"][:] = sigma[:]
        p.mesh_dict["A"][:] = a[:]
        p.mesh_dict["B"][:] = b[:]

        # Write input file
        p.write_input()

        # Run simulation
        p.run(test=True)

        # Get our results
        p.read_output()

        # Store results
        result_t = (p.ot["t"] - p.ot["t"].iloc[0]).values
        result_var1 = p.ot["theta"].values
        result_var2 = np.log10(p.ot["v_max"].values)

        self.test_results["RSF"] = {
            "t": result_t,
            "var1": result_var1,
            "var2": result_var2,
        }

        self.compare_results("RSF")
        pass


