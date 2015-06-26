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
from astropy.io import fits
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


def get_dist(chip, wavelength, position):
    """
    Get the Zernike coefficients values from file and create a sub table with the values measured on the detector at
    positions 1 to 5 (1=center, 2 to 5 are the corners starting from top-left and going counter-clockwise).
    The code uses the closest measured wavelength.

    :param chip: which one of the detector
    :type chip: int
    :param wavelength:
    :type wavelength: float in microns
    :param position: position of the PSF on the detector, for the coefficients. 0 < position < 4088
    :type position: [x,y] position on the detector
    :return: dist = list of interpolated coefficients, from Z1 to Z11
    :rtype: list of float
    """
    if glob.glob('zernike_coefficients.fits'):
        logging.debug("opening 'zernike_coefficients.fits' and reading data...")
        tbl = fits.open('zernike_coefficients.fits')
        tbl = tbl[1]  # table extensions can NOT be PrimaryHDU. tbl has a header and a data method.
        logging.debug("... done")
    else:
        logging.info("'zernike_coefficients.fits' not found. Skipping ahead.")
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # return no distortion

    wv = list()  # the measured wavelengths
    for i in range(18):
        wv.append(float(tbl.data['Wave'][i*5]))  # 5 measurements @ every wavelength
    wv.sort()
    # the actual wavelength to use (ie: the closest to measured):
    to_use = wv[min(range(len(wv)), key=lambda k: abs(wv[k]-wavelength))]
    if to_use == 2.0:
        to_use = 2  # workaround for format
    logging.debug("%sum will be used instead of %sum, to find the zernike coefficients." % (to_use, wavelength))

    rows = list()
    for i in range(len(tbl.data)):
        if tbl.data['Chip'][i] == str(chip):
            if tbl.data['Wave'][i] == str(to_use):
                rows.append(i)
    tmp_tbl = tbl.data[rows]  # sub table containing all the coefficients for *chip* and *wavelength*

    pos = [float(position[i])/100.-20.44 for i in range(2)]  # new position on a 40.88 pixel grid
    logging.debug("the reduced position of the PSF is: %s" % pos)
    local_x = tmp_tbl['Local_x']
    minx = get_min_diff(pos[0], local_x)
    local_y = tmp_tbl[tmp_tbl['Local_x'] == minx]['Local_y']
    miny = get_min_diff(pos[1], local_y)
    logging.debug("closest measured x,y: %s, %s" % (minx, miny))

    new_tbl = tmp_tbl[(tmp_tbl['Local_x'] == minx) & (tmp_tbl['Local_y'] == miny)]  # the corresponding parameters
    logging.debug("measured parameters closest to requested position: %s" % new_tbl)
    dist = list()
    for i in range(1, 12):  # from Z1 to Z11
        dist.append(float(new_tbl['Z%s' % i][0]))
    logging.debug("corresponding distortion coefficients: %s" % dist)
    return dist


def get_min_diff(pos, local):
    """
    Gets the index of the minimum difference between pos and local
    :param pos: the single position to be compared with local
    :type pos: float
    :param local: a set of position in a list format
    :type local: list of floats
    :return: the index of the minimum difference between pos and local, in local
    :rtype: int
    """

    diff = list()
    for i in range(len(local)):
        diff.append(abs(float(pos) - float(local[i])))
    pos_min = diff.index(min(diff))  # get the index of the minimum in diff, hence the minimum in local
    min_diff = local[pos_min]
    return min_diff