# Importing some required modules
import gzip
import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from numpy.testing import assert_allclose
from termcolor import colored


class AuxiliaryFunctions:

    def __init__(self):
        pass

    def compare_results(self, mode):

        if not self.frozen_loaded:
            return

        benchmark = self.frozen_results[mode]
        results = self.test_results[mode]

        t_b = benchmark["t"]
        t_r = results["t"]
        dt = np.diff(t_r)
        inds = np.arange(len(dt))[dt == 0]
        # var1_int = np.interp(t_b, t_r, results["var1"])
        # var2_int = np.interp(t_b, t_r, results["var2"])
        var1_int = interp1d(t_r, results["var1"])(t_b)
        var2_int = interp1d(t_r, results["var2"])(t_b)

        # print((benchmark["var2"]-var2_int).sum())
        # ax = plt.subplot(211)
        # plt.plot(benchmark["var2"], ".-")
        # plt.plot(var2_int, ".-")
        # plt.subplot(212, sharex=ax)
        # # plt.plot(results["var2"], ".-")
        # plt.plot(t_r, ".-")
        # plt.plot(inds, t_r[inds], ".-")
        # # plt.plot(benchmark["var2"], "k.-")
        # plt.show()
        # exit()

        try:
            assert_allclose(benchmark["var1"], var1_int, rtol=1e-4)
            assert_allclose(benchmark["var2"], var2_int, rtol=1e-4)
            self.test_results[mode]["success"] = True
            self.test_results[mode]["success_msg"] = "%s benchmark comparison... [%s]" % (mode.upper(), colored("OK", "green"))
        except AssertionError:
            self.test_results[mode]["success"] = False
            self.test_results[mode]["success_msg"] = "%s benchmark comparison... [%s]" % (mode.upper(), colored("FAILED", "red"))

        pass

    def export_results(self):

        self.test_results["hash"] = self.frozen_hash
        with gzip.GzipFile(self.pickle_file, "w") as f:
            pickle.dump(self.test_results, f, pickle.HIGHEST_PROTOCOL)
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
