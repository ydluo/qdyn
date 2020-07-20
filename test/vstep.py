# Importing some required modules
import os

import numpy as np
from aux import AuxiliaryFunctions


class TestVstep(AuxiliaryFunctions):

    test_results = {}
    dirname = "vstep"
    # The git hash below is an identifier for the version of QDYN that was
    # used to generate the benchmark results. This hash is included in the
    # frozen benchmark results, and it is checked against each time the results
    # are imported. When a new benchmark is generated, this hash should be
    # updated.
    frozen_hash = "672933dd47a6c072c8971d1fd963ced595999875"
    frozen_loaded = False

    def __init__(self, p):
        AuxiliaryFunctions.__init__(self)
        self.p = p
        pickle_name = "vstep_results_frozen_%s.tar.gz" % self.frozen_hash[:10]
        self.pickle_file = os.path.join(self.dirname, pickle_name)
        pass

    def run_test(self, mode):

        p = self.p
        set_dict = p.set_dict

        # Velocity step sequence 10 > 15 > 10 micron/s
        Vs = [10e-6, 20e-6, 10e-6]

        # Define some starting points for the first simulation
        t_final = 0

        # Total slip distance per simulation
        x_ss = 100e-6

        if mode == "CNS":
            set_dict["FRICTION_MODEL"] = "CNS"
            dict_name = "SET_DICT_CNS"
            state_ini = "PHI_INI"
            state_final = 0.28
            var2 = "tau"
            var2_ini = "TAU"
            var2_final = 2.5e6
        elif mode == "RSF":
            set_dict["FRICTION_MODEL"] = "RSF"
            dict_name = "SET_DICT_RSF"
            state_ini = "TH_0"
            state_final = 0.5*set_dict[dict_name]["DC"]/Vs[0]
            var2 = "v"
            var2_ini = "V_0"
            var2_final = Vs[0]

        result_t = np.array([], dtype=float)
        result_var1 = np.array([], dtype=float)
        result_var2 = np.array([], dtype=float)

        # Loop over all velocity steps
        for i, V in enumerate(Vs):
            # Set load-point velocity
            set_dict["V_PL"] = V
            # Set initial values from previous step
            set_dict[dict_name][state_ini] = state_final
            set_dict[dict_name][var2_ini] = var2_final
            # Set simulated time
            set_dict["TMAX"] = x_ss/V

            # Feed our settings to QDYN and render mesh
            p.settings(set_dict)
            p.render_mesh()

            # Write input file
            p.write_input()

            # Run simulation
            p.run(test=True)

            # Get our results
            p.read_output()

            # Store results
            result_t = np.hstack([result_t, p.ot[0]["t"]-p.ot[0]["t"].iloc[0]+t_final])
            result_var1 = np.hstack([result_var1, p.ot[0]["theta"]])
            result_var2 = np.hstack([result_var2, p.ot[0]["tau"]/set_dict["SIGMA"]])

            # Set starting point for next simulation
            var2_final = p.ot[0][var2].values[-1]
            state_final = p.ot[0]["theta"].values[-1]
            t_final += p.ot[0]["t"].values[-1]

        self.test_results[mode] = {
            "t": result_t,
            "var1": result_var1,
            "var2": result_var2,
        }

        self.compare_results(mode)
        pass


