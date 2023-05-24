# Importing some required modules
import gzip
import pickle

import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import pearsonr
from scipy.interpolate import interp1d
from termcolor import colored


class AuxiliaryFunctions:

    def __init__(self):
        pass

    def compare_results(self, mode, cc=False):

        if not self.frozen_loaded:
            return

        benchmark = self.frozen_results[mode]
        results = self.test_results[mode]

        t_b = benchmark["t"]
        t_r = results["t"]
        var1_int = interp1d(t_r, results["var1"], bounds_error=False, fill_value="extrapolate")(t_b)
        var2_int = interp1d(t_r, results["var2"], bounds_error=False, fill_value="extrapolate")(t_b)
        
        b1_int = interp1d(t_b, benchmark["var1"], bounds_error=False, fill_value="extrapolate")(t_b)
        b2_int = interp1d(t_b, benchmark["var2"], bounds_error=False, fill_value="extrapolate")(t_b)

        try:
            if not cc:
                assert_allclose(b1_int, var1_int, rtol=1e-4)
                assert_allclose(b2_int, var2_int, rtol=1e-4)
            else:
                cc1 = pearsonr(b1_int, var1_int, alternative="greater").statistic
                cc2 = pearsonr(b2_int, var2_int, alternative="greater").statistic
                assert cc1 > 0.85
                assert cc2 > 0.85
            self.test_results[mode]["success"] = True
            self.test_results[mode]["success_msg"] = "%s benchmark comparison... [%s]" % (mode.upper(), colored("OK", "green"))
        except AssertionError:
            self.test_results[mode]["success"] = False
            self.test_results[mode]["success_msg"] = "%s benchmark comparison... [%s]" % (mode.upper(), colored("FAILED", "red"))

        pass

    def export_results(self):

        self.test_results["hash"] = self.frozen_hash
        with gzip.GzipFile(self.pickle_file, "w") as f:
            pickle.dump(self.test_results, f)
        pass

    def import_results(self):

        with gzip.GzipFile(self.pickle_file, "r") as f:
            self.frozen_results = pickle.load(f)
        if self.frozen_results["hash"] != self.frozen_hash:
            print("\nWARNING: frozen hash does not match git hash!\n")
        self.frozen_loaded = True
        pass

    def plot_results(self, mode):

        results = self.test_results[mode]
        if self.frozen_loaded:
            frozen_results = self.frozen_results[mode]

        plt.subplot(211)
        plt.plot(results["t"], results["var2"], label="Test result")
        if self.frozen_loaded:
            plt.plot(frozen_results["t"], frozen_results["var2"],
                     "k--", label="Benchmark result")
        plt.ylabel("friction [-]")
        plt.legend(loc="upper center", ncol=2)
        plt.subplot(212)
        plt.plot(results["t"], results["var1"])
        if self.frozen_loaded:
            plt.plot(frozen_results["t"], frozen_results["var1"], "k--")
        plt.ylabel("state")
        plt.xlabel("time [s]")
        plt.tight_layout()
        plt.show()

        pass
