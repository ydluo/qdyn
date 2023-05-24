# Importing some required modules
import os

import numpy as np
from aux import AuxiliaryFunctions


class TestSingleAsperity(AuxiliaryFunctions):

    test_results = {}
    dirname = "singleasperity"
    # The git hash below is an identifier for the version of QDYN that was
    # used to generate the benchmark results. This hash is included in the
    # frozen benchmark results, and it is checked against each time the results
    # are imported. When a new benchmark is generated, this hash should be
    # updated.
    frozen_hash = "ed8645b9584ab24be87f908143d13dce0a336a71"
    frozen_loaded = False

    def __init__(self, p):
        AuxiliaryFunctions.__init__(self)
        self.p = p
        pickle_name = "singleasperity_results_frozen_%s.tar.gz" % self.frozen_hash[:10]
        self.pickle_file = os.path.join(self.dirname, pickle_name)
        pass

    @staticmethod
    def render_RSF_mesh(set_dict):

        b = set_dict["SET_DICT_RSF"]["B"]
        Dc = set_dict["SET_DICT_RSF"]["DC"]
        mu = set_dict["MU"]
        sigma = set_dict["SIGMA"]
        V_pl = set_dict["V_PL"]

        Lasp = 7
        L = 5
        ab_ratio = 0.9
        cab_ratio = 1 - ab_ratio
        resolution = 7

        Lb = mu*Dc/(b*sigma)
        Lc = Lb / cab_ratio
        Lasp *= Lc
        L *= Lasp

        N = int(np.power(2, np.ceil(np.log2(resolution*L/Lb))))
        x = np.linspace(-L/2, L/2, N, dtype=float)

        a = b*(1 + cab_ratio*(1 - 2*np.exp(-(2*x/Lasp)**6)))

        V0 = V_pl*(1 + 0.01*np.exp(-(2*x/Lasp)**6))
        V0 = V0 * V_pl / np.mean(V0)

        result = {
            "N": N,
            "L": L,
            "a": a,
            "V0": V0,
        }

        return result

    def run_test(self, mode):

        p = self.p
        set_dict = p.set_dict
        t_yr = 3600 * 24 * 365.0

        set_dict["MESHDIM"] = 1
        set_dict["FINITE"] = 0
        set_dict["TMAX"] = 10*t_yr
        set_dict["NTOUT"] = 1000
        set_dict["NXOUT"] = 1
        set_dict["V_PL"] = 1e-9
        set_dict["MU"] = 3e10
        set_dict["W"] = 50e3
        set_dict["SIGMA"] = 1e8
        set_dict["ACC"] = 1e-10
        set_dict["DTTRY"] = 100
        set_dict["V_TH"] = 1e-2
        set_dict["SOLVER"] = 1

        # Setting some RSF parameters
        set_dict["SET_DICT_RSF"]["A"] = 0.9e-2
        set_dict["SET_DICT_RSF"]["B"] = 1e-2
        set_dict["SET_DICT_RSF"]["DC"] = 4e-4
        set_dict["SET_DICT_RSF"]["V_SS"] = set_dict["V_PL"]
        set_dict["SET_DICT_RSF"]["TH_0"] = \
            set_dict["SET_DICT_RSF"]["DC"] / set_dict["V_PL"]

        # Setting some CNS parameters
        set_dict["SET_DICT_CNS"]["IPS_CONST_DISS"] = 1e-15
        set_dict["SET_DICT_CNS"]["A"] = [1e-15]

        result = self.render_RSF_mesh(set_dict)
        set_dict["N"] = result["N"]
        set_dict["L"] = result["L"]

        # Setting some CNS parameters
        set_dict["SET_DICT_CNS"]["PHI_INI"] = 0.1
        set_dict["SET_DICT_CNS"]["TAU"] = 0.2*set_dict["SIGMA"]

        if mode == "RSF":
            set_dict["FRICTION_MODEL"] = "RSF"

            p.settings(set_dict)
            p.render_mesh()

            p.mesh_dict["A"] = result["a"]
            p.mesh_dict["V_0"] = result["V0"]

        if mode == "CNS":
            set_dict["FRICTION_MODEL"] = "CNS"
            set_dict["SOLVER"] = 2

            p.settings(set_dict)
            p.render_mesh()

            if "IPS_CONST_DISS" in p.mesh_dict:
                p.mesh_dict["IPS_CONST_DISS"][1500:2500] = 1e-12
            elif "A0" in p.mesh_dict:
                p.mesh_dict["A0"][1500:2500] = 1e-12
            else:
                print("Unable to set CNS creep kinetics")
                exit()


        # Write input file
        p.write_input()

        # Run simulation
        p.run(test=True)

        # Get our results
        p.read_output()

        # Store results
        result_t = (p.ot[0]["t"] - p.ot[0]["t"].iloc[0]).values
        result_var1 = p.ot[0]["theta"].values
        result_var2 = np.log10(p.ot_vmax["v"].values)

        self.test_results[mode] = {
            "t": result_t,
            "var1": result_var1,
            "var2": result_var2,
        }

        self.compare_results(mode, cc=True)
        pass


