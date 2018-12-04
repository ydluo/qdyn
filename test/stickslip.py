# Importing some required modules
import os

from aux import AuxiliaryFunctions


class TestStickSlip(AuxiliaryFunctions):

    test_results = {}
    dirname = "stickslip"
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
        pickle_name = "stickslip_results_frozen_%s.tar.gz" % self.frozen_hash[:10]
        self.pickle_file = os.path.join(self.dirname, pickle_name)
        pass

    def run_test(self, mode):

        p = self.p
        set_dict = p.set_dict

        set_dict["V_PL"] = 1e-6
        set_dict["TMAX"] = 1e4
        set_dict["MU"] = 2e10
        set_dict["L"] = 1e1

        if mode == "CNS":
            set_dict["FRICTION_MODEL"] = "CNS"
        elif mode == "RSF":
            set_dict["FRICTION_MODEL"] = "RSF"
            dict_name = "SET_DICT_RSF"
            set_dict[dict_name]["A"] = 0.01
            set_dict[dict_name]["B"] = 0.025

        p.settings(set_dict)
        p.render_mesh()

        # Write input file
        p.write_input()

        # Run simulation
        p.run(test=True)

        # Get our results
        p.read_output()

        # Store results
        result_t = (p.ot["t"]-p.ot["t"].iloc[0]).values
        result_var1 = p.ot["theta"].values
        result_var2 = p.ot["tau"].values/set_dict["SIGMA"]

        self.test_results[mode] = {
            "t": result_t,
            "var1": result_var1,
            "var2": result_var2,
        }

        self.compare_results(mode)
        pass
