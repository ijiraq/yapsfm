#! /usr/bin/env python
"""
Contains functions to compute Zernike modes.

Inputs:
    n: mode's radial
    m: mode's asimuthal
    r: a meshgrid object of radial positions
    theta: a meshgrid object of asimuthal positions

Outputs:
    an array giving the shape of the wavefront
"""

import numpy as np
import pyfits
import glob
import logging


def even_zernike(n, m, r, theta):
    """
    Computes the even Zernike mode following Noll (1976)

    :param n: radial order of the mode
    :type n: int
    :param m: azimuthal order of the mode
    :type m: int
    :param r: radial position
    :type r: float
    :param theta: angle
    :type theta: float
    :return: the value of the Zernike mode corresponding to n and m, at position (r,theta).
    :rtype: float
    """

    zernike_value = np.sqrt(n+1)*r_nm(n, m, r)*np.sqrt(2)*np.cos(m*theta)
    return zernike_value


def odd_zernike(n, m, r, theta):
    """
    Computes the odd Zernike mode following Noll (1976)

    :param n: radial order of the mode
    :type n: int
    :param m: azimuthal order of the mode
    :type m: int
    :param r: radial position
    :type r: float
    :param theta: angle
    :type theta: float
    :return: the value of the Zernike mode corresponding to n and m, at position (r,theta).
    :rtype: float
    """

    zernike_value = np.sqrt(n+1)*r_nm(n, m, r)*np.sqrt(2)*np.sin(m*theta)
    return zernike_value


def r_nm(n, m, r):
    """
    Computes the radial part R_n^m(r), with r a meshgrid object

    :param n: radial order of the polynomial
    :type n: int
    :param m: azimuthal order of the polynomial
    :type m: int
    :param r: radial position
    :type r: float
    :return: radial part of the Zernike mode at position r
    :rtype: float
    """

    r_nm_value = 0
    for s in range((n-m)/2+1):
        r_nm_value += (((-1)**s * np.math.factorial(n-s)) / (np.math.factorial(s) * np.math.factorial((n+m)/2-s) *
                                                             np.math.factorial((n-m)/2-s))) * r**(n-2*s)
    return r_nm_value


def select_coeff(chip, wavelength, position):
    """
    Get the Zernike coefficients values from file and create a sub table with the values measured on the detector at
    positions 1 to 5 (1=center, 2 to 5 are the corners starting from top-left and going counter-clockwise).

    :param chip: which one of the detector
    :type chip: int
    :param wavelength:
    :type wavelength: float in microns
    :param position: where on the detector, for interpolation of the coefficient
    :type position: tuple of (x,y) position on the detector
    :return: list of interpolated coefficients, from Z1 to Z11
    :rtype: list of float
    """
    if glob.glob('zernike_coefficients.fits'):
        logging.debug("opening 'zernike_coefficients.fits' and reading data...")
        tbl = pyfits.open('zernike_coefficients.fits')
        tbl = tbl[1]  # table extensions can NOT be PrimaryHDU. tbl has a header and a data method.
        logging.debug("... done")
    else:
        logging.info("'zernike_coefficients.fits' not found. Skipping ahead.")
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # return no distortion

    rows = list()
    for i in range(len(tbl.data)):
        if tbl.data['Chip'][i] == str(chip):
            if tbl.data['Wave'][i] == str(wavelength):
                rows.append(i)
    new_tbl = tbl.data[rows]  # sub table containing all the coefficients for *chip* and *wavelength*

    ##########
    # interpolation part:
    # for each zernike mode, interpolate the 5 coefficients on a 4088 x 4088 grid times the (over)sampling factor,
    # using spline (order 3 ?)
    # (over)sampling factor = 0.11 / scale_factor (both in ''/px)
    ##########

    dist = list()
    return dist