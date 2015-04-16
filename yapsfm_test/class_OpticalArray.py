#! /usr/bin/env python
"""
The common ancestor class between aperture, pupil and the PSF object is OpticalArray: a ImageHDU object with
properties *array_size*, *center*, *a* and *name*.
One can save the OpticalArray instance in a fits format using *save*.

Aperture inherits from OpticalArray and has methods *circular* to make a circular aperture or *make_hst_ap* to make an
HST-like aperture. Furthermore, one can *reset* the aperture to circular shape, or *fits2aperture* to load an aperture
from a .fits image.
"""

import numpy as np
import pyfits


class OpticalArray(object):
    def __init__(self, size):
        self.array_size = size
        self.center = [size/2, size/2]
        self.a = np.zeros((size, size))
        self.name = ''

    def save(self):
        hdu = pyfits.PrimaryHDU(self.a)
        hdu.writeto('%s.fits' % self.name, clobber=True)
        print 'saving %s' % self.name


class Aperture(OpticalArray):
    """
    An aperture class, handling simple circular, HST-like or input file (.fits)
    """
    def __init__(self, size):
        super(Aperture, self).__init__(size)
        self._circular()

    def circular(self):
        """
        Makes a simple circular aperture of size *self.array_size*
        """
        for y in range(self.array_size):
            for x in range(self.array_size):
                radial_position = np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
                if radial_position <= self.array_size/2:
                    self.a[y][x] = 1.
        self.name = 'circular_aperture'
        return

    def make_hst_ap(self):
        """
        Makes an HST-like aperture of size *self.array_size*
        """
        sec_mir = 0.330*self.array_size/2  # secondary mirror
        sp_width = 0.022*self.array_size/2  # spider width
        pad1 = [0.8921*self.array_size/2, 0.0000, 0.065*self.array_size/2]  # x, y, radius
        pad2 = [-0.4615*self.array_size/2, 0.7555*self.array_size/2, 0.065*self.array_size/2]
        pad3 = [-0.4564*self.array_size/2, -0.7606*self.array_size/2, 0.065*self.array_size/2]
        for y in range(self.array_size):
            for x in range(self.array_size):
                # main aperture (including secondary obstruction)
                radial_position = np.sqrt((x - self.center[0])**2+(y - self.center[1])**2)
                if self.array_size/2 >= radial_position > sec_mir:
                    self.a[y][x] = 1.
                # Obstructions:
                # spiders
                if self.center[0] - sp_width/2. <= x <= self.center[0]+sp_width/2:
                    self.a[y][x] = 0.
                if self.center[0]-sp_width/2 <= y <= self.center[1]+sp_width/2:
                    self.a[y][x] = 0.
                # pads
                if np.sqrt((x-self.center[0]-pad1[0])**2+(y-self.center[1]-pad1[1])**2) <= pad1[2]:
                    self.a[y][x] = 0.
                if np.sqrt((x-self.center[0]-pad2[0])**2+(y-self.center[1]-pad2[1])**2) <= pad2[2]:
                    self.a[y][x] = 0.
                if np.sqrt((x-self.center[0]-pad3[0])**2+(y-self.center[1]-pad3[1])**2) <= pad3[2]:
                    self.a[y][x] = 0.
        self.name = 'hst_aperture'
        return

    def reset(self):
        """
        Reset the aperture to simple circular.
        """
        for y in range(self.array_size):
            for x in range(self.array_size):
                radial_position = np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
                if radial_position <= self.array_size/2:
                    self.a[y][x] = 1.
        self.name = 'circular_aperture'
        return

    def fits2aperture(self):
        """
        Opens an aperture fits file and reads it.
        """
        hdu = pyfits.open('aperture.fits')
        if len(hdu) < 2:  # check if aperture.fits has a header (it should not)
            self.a = hdu[0].data
        else:
            self.a = hdu[1].data

        self.array_size = np.shape(self.a)[0]
        self.name = 'aperture'
        return


class Distorted(OpticalArray):
    """
    The ancestor class of both *Pupil* and *PSF* objects.
    Based on *OpticalArray*, but have extra attributes *wavelength* and optical distortions *dist*
    """
    def __init__(self, size, wavelength, dist):
        super(OpticalArray, self).__init__(self, size)
        self.wavelength = wavelength
        self.dist = dist


class Pupil(Distorted):
    def __init__(self):
        self.wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns') or .76)
        self.dist = self._read_distortions()
        self.array_size = 101
        super(Distorted, self).__init__(self)
        self._compute_pupil(101)
        self.name = 'Pupil'

    def _read_distortions(self, file_path='distortions.par'):
        data = open(file_path).readlines()

        dist = []
        for i, line in enumerate(data):
            dist.append(float(data[i].partition('#')[0].strip()))
        return dist

    def _compute_pupil(self):
        print 'Computing pupil...'

        aperture = Aperture(self.array_size)
        if glob.glob('aperture.fits'):
            print "... using 'aperture.fits'..."
            aperture.fits2aperture()
        else:
            aperture.make_hst_ap()

        pass
        opd = path_diff(self.array_size, self.wavelength, self.dist)  # the optical path difference (OPD)
        self.a = np.multiply(self.a, np.exp(np.divide(2j*np.pi*opd, self.wavelength)))
        print '... done'