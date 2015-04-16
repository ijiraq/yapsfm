#! /usr/bin/env python
"""
Test class for pupil computation
"""

import numpy as np
import glob
import class_aperture as ca
import scipy.ndimage.interpolation


class Pupil(object):
    def __init__(self, size, wavelength, dist):
        self.size = size
        self.center = [size/2, size/2]
        self.wavelength = wavelength
        self.a = np.zeros((size, size))
        self.dist = dist

    def compute_pupil(self):
        print 'Computing pupil...'

        aperture = ca.Aperture(self.size)
        if glob.glob('aperture.fits'):
            print "... using 'aperture.fits'..."
            aperture.fits2aperture()
        else:
            aperture.make_hst_ap()

        opd = path_diff(self.size, self.wavelength, self.dist)  # the optical path difference (OPD)
        self.a = np.multiply(self.a, np.exp(np.divide(2j*np.pi*opd, self.wavelength)))
        print '... done'


def path_diff(size, wavelength, dist):
    zernike_modes = [(1, 1), (1, -1), (2, 0), (2, -2), (2, 2), (3, -1), (3, 1), (3, -3), (3, 3), (4, 0)]  # z2..z11

    rho_range = np.linspace(0, 1, size)
    theta_range = np.linspace(0, 2*np.pi, size)
    rho, theta = np.meshgrid(rho_range, theta_range)

    zernike_total = np.zeros((size, size))
    for i in range(len(dist)):
        aj = dist[i]*.547/wavelength  # Zernike coefficient in microns, .547um is the reference wavelength
        print 'Computing Z%s with aj=%s' % (4+i, aj)
        n, m = zernike_modes[i][0], zernike_modes[i][1]
        if m < 0.:
            zernike_value = odd_zernike(n, -m, rho, theta)
        else:
            zernike_value = even_zernike(n, m, rho, theta)
        zernike_total += np.multiply(zernike_value, aj)  # OPD = Sum aj Zj
    print 'OPD computed...'
    cartesian_z = scipy.ndimage.interpolation.geometric_transform(zernike_total, polar2cart, extra_arguments=(size,))
    print '... and converted back to cartesian space.'

    return cartesian_z


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


def polar2cart(coords, size=101):
    """
    Change between polar and cartesian coordinates system

    :param coords: list of array position
    :type coords: np.array
    :param size: total size of the image
    :type size: int
    :return: the (theta,r) values in (x,y)
    """

    # origin back at the center of the image
    x = (coords[1]-size//2.)/(size//2.)
    y = (coords[0]-size//2.)/(size//2.)

    # compute -1<r<1 and 0<theta<2pi
    r = np.sqrt(x*x+y*y)
    theta = np.arctan2(y, x)
    theta = theta < 0 and theta+2*np.pi or theta

    # bin r,theta back to pixel space (size,size)
    r *= size-1
    theta *= (size-1)/(2*np.pi)
    return theta, r