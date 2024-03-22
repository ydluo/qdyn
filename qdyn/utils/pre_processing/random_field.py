#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Random functions for QDYN

Generate 2-D random field using a power spectral density P(k), k is 
wavenumber.

Created on Fri Jul  2 11:26:01 2021
@author: yifan
"""

#from FyeldGenerator import generate_field
from scipy.special import kv
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
import six


def psd_gen(coef):
    def psd_vonKarman(k):
        """
        Parameters
        ----------
        coef : list-like
            coef contains [ax, az, H]
        H : float
            The Hurst exponent, 0 <= H <= 1.
        """
        ax, az, H = coef
        k2 = ax**2*k**2 + az**2*k**2
        return (4 * np.pi * H * ax * az / kv(H, 0.1)) / (1 + k2)**(H+1)
    return psd_vonKarman


def phase_distrib(shape, rseed):
    rng = np.random.default_rng(seed=rseed)
    a = rng.uniform(-1, 1, shape)
    b = rng.uniform(-1, 1, shape)
    return a + 1j * b


# FyeldGenerator
def generate_field(power_spectrum, statistic, rseed, shape, unit_length=1,
                   fft=np.fft, fft_args=dict(), stat_real=False):
    """
    Generates a field given a stastitic and a power_spectrum.
    Parameters.

    Author : FyeldGenerator
    
    ----------
    statistic: callable
        A function that takes returns a random array of a given signature,
        with signature (s) -> (B) with B.shape == s. Please note that the
        distribution is in *Fourier space* not in real space, unless you set
        stat_real=True. See the note below for more details.
    power_spectrum: callable
        A function that returns the power contained in a given mode,
        with signature (k) -> P(k) with k.shape == (ndim, n)
    shape: tuple
        The shape of the output field
    unit_length: float, optional
        How much physical length represent 1pixel. For example a value of 10
        mean that each pixel stands for 10 physical units. It has the
        dimension of a physical_unit/pixel.
    fft: a numpy-like fft API, optional
    fft_args: array, optional
        a dictionary of kwargs to pass to the FFT calls
    stat_real: boolean, optional
        Set to true if you want the distribution to be drawn in real space and
        then transformed into Fourier space.
    Returns
    -------
    field: a real array of shape `shape` following the statistic
        with the given power_spectrum
    Note
    ----
    When generation the distribution in Fourier mode, the result
    should be complex and unitary. Only the phase is random.
    """

    if not six.callable(statistic):
        raise Exception('`statistic` should be callable')
    if not six.callable(power_spectrum):
        raise Exception('`power_spectrum` should be callable')

    try:
        fftfreq = fft.fftfreq
        rfftfreq = fft.rfftfreq
    except NameError:
        # Fallback on numpy for the frequencies
        fftfreq = np.fft.fftfreq
        rfftfreq = np.fft.rfftfreq

    # Compute the k grid
    all_k = [fftfreq(s, d=unit_length) for s in shape[:-1]] + \
            [rfftfreq(shape[-1], d=unit_length)]
#    print(all_k)
    kgrid = np.meshgrid(*all_k, indexing='ij')
    knorm = np.sqrt(np.sum(np.power(kgrid, 2), axis=0))
#    print(knorm)
    fourier_shape = knorm.shape

    if stat_real:
        field = statistic(shape, rseed)
        # Compute the FFT of the field
        fftfield = fft.rfftn(field, **fft_args)
    else:
        # Draw a random sample in Fourier space
        fftfield = statistic(fourier_shape, rseed)

    power_k = np.where(knorm == 0, 0, np.sqrt(power_spectrum(knorm)))
    fftfield *= power_k

    return np.real(fft.irfftn(fftfield, **fft_args))


def generate_Dc(dc_min, dc_max, shape, coef, rseed, plot=False):
    """Generate Dc field using the min and max inputs. The min and max are the
    natural log of the Dc length. The field is log-normal when generated, then
    transformed to log-uniform. The returned field is transfered back in
    meters.
    """

    field = generate_field(psd_gen(coef), phase_distrib, rseed, shape)
    floc, fscale = norm.fit(field)
    field_cdf = norm.cdf(field, floc, fscale)
#    refield = np.exp(field_cdf * (dc_max - dc_min) + dc_min)
    refield = 10**(field_cdf * (dc_max - dc_min) + dc_min)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(refield, cmap='Spectral')
        fig.colorbar(im, ax=ax)
        plt.show()
    return refield
    
# def psd_gaussian(ax, az, kr):
#     return 0.5*ax*az*np.exp(-0.25*kr**2)

# %%
if __name__ == '__main__':

    # def Pkgen(n):
    #     def Pk(k):
    #         return np.power(k, -n)
    
    #     return Pk
    
    # # Draw samples from a normal distribution
    # def distrib(shape):
    #     a = np.random.normal(loc=0, scale=1, size=shape)
    #     b = np.random.normal(loc=0, scale=1, size=shape)
    #     return a + 1j * b

    # Define the two functions

    shape = (120, 160)
    coef = [10, 10, 0.8]
    rseed = 42
    Dc_max = 0.1
    Dc_min = 0.015


#    field = generate_field(distrib, Pkgen(2), shape)

#   Rescale the Dc distribution from normal to uniform

