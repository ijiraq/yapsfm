#! /usr/bin/env python
"""
OpticalArray is the common ancestor class between aperture, pupil and the PSF object: it's an object with
properties *array_size*, *center*, *name*, *_a* and *_polar*.
One can save the OpticalArray instance in a fits format using *save* or change it from polar coordinate system to
cartesian using *polar2cart*.
"""

import pyfits
import numpy as np
import scipy.ndimage.interpolation
from datetime import datetime as dt


class OpticalArray(object):
    """
    The ancestor class, it has the basic properties *array_size*, *center*, *name*, the data array *_a*, and a polar
    coordinates trigger *_polar*, allowing to switch between polar and cartesian coordinate systems.
    """
    def __init__(self, size, polar=False):
        self.array_size = size
        self.center = [self.array_size//2., self.array_size//2.]
        self.name = ''
        self._a = None
        self._polar = polar

    def save(self):
        hdu = pyfits.PrimaryHDU(self.a)
        header = hdu.header
        now = dt.now()
        created = '%a %b %d %X %Y'
        header['CREATED'] = (now.strftime(created), 'Time and date file was created')

        """
        if dist:
            now = dt.now()
            created = '%a %b %d %X %Y'
            header['CREATED'] = (now.strftime(created), 'Time and date file was created')
            header['INSTRUME'] = ('WFI', 'Simulated instrument')
            header['FOCUS'] = (dist[0], 'PSF RMS focus (waves @ 547 nm)')
            header['X_ASTIG'] = (dist[1], 'PSF RMS 0d astig (waves @ 547 nm)')
            header['Y_ASTIG'] = (dist[2], 'PSF RMS 45d astig (waves @ 547 nm)')
            header['X_COMA'] = (dist[3], 'PSF RMS X-coma (waves @ 547 nm)')
            header['Y_COMA'] = (dist[4], 'PSF RMS Y-coma (waves @ 547 nm)')
            header['X_CLOVER'] = (dist[5], 'PSF RMS X-clover (waves @ 547 nm)')
            header['Y_CLOVER'] = (dist[6], 'PSF RMS Y-clover (waves @ 547 nm)')
            header['SPHEICL'] = (dist[7], 'PSF RMS spherical (waves @ 547 nm)')
            header['PIXSCALE'] = ('work in progress', 'Pixel scale in arcseconds')
            header['WAVELNTH'] = (wavelength, 'PSF wavelength in microns')
        """
        hdu.writeto('%s.fits' % self.name, clobber=True)
        print 'saving %s' % self.name

    @property
    def a(self):
        if self._a is not None:
            return self._a
        self._a = np.zeros((self.array_size, self.array_size))
        return self._a

    @a.setter
    def a(self, a):
        self._a = a

    def polar2cart(self):
        """
        Change between polar and cartesian coordinates system

        :return: the (theta,r) values in (x,y)
        """
        if not self._polar:
            return
        size = self.array_size

        def func(coords):

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
        self.a = scipy.ndimage.interpolation.geometric_transform(self.a, func)
        self._polar = False