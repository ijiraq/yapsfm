#! /usr/bin/env python
"""
Pupil inherits from Distorted and has properties *opd* and *name*. Its methods are *_compute_pupil*, *_path_diff*
"""

import glob
import numpy as np

from distorted import Distorted
from aperture import Aperture
from optiarray import OpticalArray
import zernike as ze


class Pupil(Distorted):
    """
    Pupil inherits from Distorted, it has properties *opd* (the optical path differences) and *name*.
    It can be computed using *_compute_pupil* which will compute the *_path_diff* of the wavefront.
    """
    def __init__(self, wavelength, array_size=101):
        super(Pupil, self).__init__(array_size, wavelength)
        self.opd = self._path_diff()
        self._compute_pupil()
        self.name = 'Pupil'

    def _compute_pupil(self):
        print 'Computing pupil...'

        aperture = Aperture(self.array_size)
        if glob.glob('aperture.fits'):
            print "... using 'aperture.fits'..."
            aperture.fits2aperture()
        else:
            aperture.make_hst_ap()

        pass
        self.a = np.multiply(aperture.a, np.exp(np.divide(2j*np.pi*self.opd, self.wavelength)))
        print '... done'

    def _path_diff(self):
        zernike_modes = [(1, 1), (1, -1), (2, 0), (2, -2), (2, 2), (3, -1), (3, 1), (3, -3), (3, 3), (4, 0)]  # z2..z11

        rho_range = np.linspace(0, 1, self.array_size)
        theta_range = np.linspace(0, 2*np.pi, self.array_size)
        rho, theta = np.meshgrid(rho_range, theta_range)

        zernike_total = OpticalArray(polar=True, size=self.array_size)
        for i in range(len(self.dist)):
            aj = self.dist[i]*.547/self.wavelength  # Zernike coefficient in microns, .547um is the reference wavelength
            print 'Computing Z%s with aj=%s' % (2+i, aj)
            n, m = zernike_modes[i][0], zernike_modes[i][1]
            if m < 0.:
                zernike_value = ze.odd_zernike(n, -m, rho, theta)
            else:
                zernike_value = ze.even_zernike(n, m, rho, theta)
            zernike_total.a += np.multiply(zernike_value, aj)  # OPD = Sum aj Zj
        print 'OPD computed...'
        zernike_total.polar2cart()
        print '... and converted back to cartesian space.'

        return zernike_total.a