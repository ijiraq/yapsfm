#! /usr/bin/env python
"""
OpticalArray is the common ancestor class between aperture, pupil and the PSF object: it's an object with
properties *array_size*, *center*, *name*, *_a*, *_polar*, *wavelength* and *_dist*.
One can save the OpticalArray instance in a fits format using *save* or change it from polar coordinate system to
cartesian using *polar2cart*.

Aperture inherits from OpticalArray and has methods *circular* to make a circular aperture or *make_hst_ap* to make an
HST-like aperture. Furthermore, one can *reset* the aperture to circular shape, or *fits2aperture* to load an aperture
from a .fits image.

Pupil inherits from Distorted and has properties *opd* and *name*. Its methods are *_compute_pupil*, *_path_diff*

PSF inherits from Distorted and has properties *array_size*, *pupil*, *_a* and *name*. Its method is *rezie_psf* to
a desired resolution.
"""

import glob
import pyfits
import numpy as np
import scipy.ndimage.interpolation
from datetime import datetime as dt
import zernike as ze


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
        self.wavelength = None
        self._dist = None

    @property
    def dist(self, file_path='distortions.par'):
        if self._dist is not None:
            return self._dist
        data = open(file_path).readlines()

        self._dist = []
        for i in range(len(data)):
            self._dist.append(float(data[i].partition('#')[0].strip()))
        return self._dist

    def save(self):
        hdu = pyfits.PrimaryHDU(self.a)
        header = hdu.header
        now = dt.utcnow()
        created = "%Y-%m-%dT%H:%M:%S"
        header['DATE'] = (now.strftime(created), 'UTC time and date file was created')

        """
        # Work in progress
        if dist:
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


class Aperture(OpticalArray):
    """
    An aperture class, handling simple circular, HST-like or input file (.fits)
    """
    def __init__(self, size):
        super(Aperture, self).__init__(size)
        self.circular()

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


class Pupil(OpticalArray):
    """
    Pupil inherits from Distorted, it has properties *opd* (the optical path differences) and *name*.
    It can be computed using *_compute_pupil* which will compute the *_path_diff* of the wavefront.
    """
    def __init__(self, wavelength, array_size=101):
        super(Pupil, self).__init__(array_size)
        self.wavelength = wavelength
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


class PSF(OpticalArray):
    """
    PSF inherits from Distorted and has properties *array_size*, *pupil*, *_a* the array containing the data, and
    *name*. Its resolution can be modified using *resize_psf*, which will change the pixel resolution to *scale*.
    """
    def __init__(self, pupil):
        self.array_size = 505
        self.pupil = pupil
        self._a = None
        super(PSF, self).__init__(self.array_size, pupil.wavelength)
        self.name = 'Psf'

    @property
    def a(self):
        if self._a is not None:
            return self._a
        # print 'Starting FFT with zero-padding factor of %s...' % (self.array_size/np.shape(pupil.a)[0])
        tmp = np.fft.fft2(self.pupil.a, s=[self.array_size, self.array_size])  # padding with 0s
        tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array
        print "... done"
        self._a = np.real(np.multiply(tmp, np.conjugate(tmp)))
        return self._a

    def resize_psf(self, wavelength=.76, size=505, scale=0.110):
        print "resizing PSF to match pixel resolution of %s''/px..." % scale
        new_psf = scipy.ndimage.interpolation.geometric_transform(self._a, rebin, extra_arguments=(
            wavelength, size, scale))
        # new_psf = new_psf[size//2.-32:size//2.+32, size//2.-32:size//2.+32]
        print '... done'
        self._a = new_psf


def rebin(coords, wavelength, size, detector_scale):
    """
    Scale the PSF to the desired size (detector_scale, in ''/pixel)
    This function is called by *resize_psf* in geometric_transform.

    :param coords: list of array coordinates
    :param wavelength: the wavelength at which the PSF is computed at
    :param size: size of the array
    :param detector_scale: the pixel resolution at which to rebin the PSF
    :return: new coordinates (python starts with y)
    """
    scale_factor = 0.0797/6.*(wavelength/0.76)  # at 5x zero-padding

    x = (coords[1]-size//2.)*detector_scale/scale_factor+(size//2.)
    y = (coords[0]-size//2.)*detector_scale/scale_factor+(size//2.)
    return y, x