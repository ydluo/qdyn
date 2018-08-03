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
    frozen_hash = "f46e45e06a647519845c6f8e22bb5cdfd8c2a2f7"
    frozen_loaded = False

    def __init__(self, p):
        AuxiliaryFunctions.__init__(self)
        self.p = p
        pickle_name = "singleasperity_results_frozen_%s.tar.gz" % self.frozen_hash[:10]
        self.pickle_file = os.path.join(self.dirname, pickle_name)
        pass

    def run_test(self, mode):

        p = self.p
        set_dict = p.set_dict

        N = np.power(2, 10)
        t_yr = 3600 * 24 * 365.0

        set_dict["N"] = N
        set_dict["MESHDIM"] = 1
        set_dict["FINITE"] = 0
        set_dict["TMAX"] = 100*t_yr
        set_dict["L"] = 500.0
        set_dict["V_PL"] = 1e-9

        if mode == "CNS":
            set_dict["FRICTION_MODEL"] = "CNS"
            dict_name = "SET_DICT_CNS"
            set_dict[dict_name]["PHI_INI"] = 0.1
            set_dict[dict_name]["TAU"] = 0.2*set_dict[dict_name]["SIGMA"]

            # Feed our settings to QDYN and render mesh
            p.settings(set_dict)
            p.render_mesh()

        elif mode == "RSF":
            set_dict["FRICTION_MODEL"] = "RSF"
            dict_name = "SET_DICT_RSF"
            set_dict[dict_name]["A"] = 0.01
            set_dict[dict_name]["B"] = 0.009
            set_dict[dict_name]["DC"] = 1e-4
            set_dict[dict_name]["TH_0"] = 0.5*set_dict[dict_name]["DC"]/set_dict["V_PL"]

            # Feed our settings to QDYN and render mesh
            p.settings(set_dict)
            p.render_mesh()

            p.mesh_dict["B"][2*N//7:5*N//7] = 0.015

        # Write input file
        p.write_input()

        # Run simulation
        p.run(test=False)

        # Get our results
        p.read_output()

        self.test_results[mode] = {
            "t": result_t,
            "var1": result_var1,
            "var2": result_var2,
        }

        self.compare_results(mode)
        pass


