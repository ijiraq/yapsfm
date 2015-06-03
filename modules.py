#! /usr/bin/env python
"""
OpticalArray is the common ancestor class between aperture, pupil and the PSF object: it's an object with
properties *array_size*, *center*, *name*, *_a*, *_polar*, *wavelength* and *_dist*.
One can save the OpticalArray instance in a fits format using *save* or change it from polar coordinate system to
cartesian using *polar2cart*.

Aperture inherits from OpticalArray and has methods *circular* to make a circular aperture or *make_hst_ap* to make an
HST-like aperture. Furthermore, one can *reset* the aperture to circular shape, or *fits2aperture* to load an aperture
from a .fits image.

Pupil inherits from OpticalArray and has properties *opd* and *name*. Its methods are *_compute_pupil*, *_path_diff*

PSF inherits from OpticalArray and has properties *array_size*, *pupil*, *_a* and *name*. Its method is *rezie_psf* to
a desired resolution.

PolyPSF inherits from OpticalArray and has properties *band*, *spectral_type*, *size*, *scale*. Its methods are
*get_sed*, *wavelength_contribution*, *create_polychrome* and *check_sed*.
"""

import glob
import sys
import pyfits
import logging
import numpy as np
import scipy.ndimage.interpolation
from scipy.interpolate import interp1d
from datetime import datetime as dt
import zernike as ze


class OpticalArray(object):
    """
    The ancestor class, it has the basic properties *array_size*, *center*, *name*, the data array *_a*, and a polar
    coordinates trigger *_polar*, allowing to switch between polar and cartesian coordinate systems.

    init takes: size, polar=False, scale=None, poly=False
    """
    def __init__(self, size, polar=False, scale=None, poly=False):
        self.array_size = size
        self.center = [self.array_size//2., self.array_size//2.]
        self.name = ''
        self._a = None
        self._polar = polar
        self.wavelength = None
        self._dist = None
        self._poly = poly
        self.band = None
        self.spectral_type = None
        self._wavelength_contributions = None
        self.scale = scale
        self.b = None

    @property
    def dist(self, file_path='distortions.par'):
        if self._dist is not None:
            return self._dist
        data = open(file_path).readlines()

        self._dist = []
        for i in range(len(data)):
            self._dist.append(float(data[i].partition('#')[0].strip()))
        logging.debug("distortions (@ 0.547um): %s" % self._dist)
        logging.debug("distortions (@ %sum): %s" % (self.wavelength, [self._dist[i]*.547/self.wavelength for i in
                                                    range(len(self._dist))]))
        return self._dist

    def save(self, name=None):
        """
        Will save the data array "a" in a fits file with header information as follow:
        :param name: name of the .fits file to be written. If None, will call it self.name
        :return: None
        """

        if self.b is None:
            hdu = pyfits.PrimaryHDU(self.a)
        else:
            hdu = pyfits.PrimaryHDU()
            list_arrays = list()
            list_arrays.append(self.a)
            for i in self.b:
                list_arrays.append(i)
            hdu.data = np.array(list_arrays)
        header = hdu.header
        now = dt.utcnow()
        created = "%Y-%m-%dT%H:%M:%S"
        header['DATE'] = (now.strftime(created), 'Date and UTC time file was created')

        if self._dist:
            header['INSTRUME'] = ('work in progress', 'Simulated instrument')
            header['FOCUS'] = (self._dist[2], 'PSF RMS focus (waves @ 547 nm)')
            header['X_ASTIG'] = (self._dist[3], 'PSF RMS 0d astig (waves @ 547 nm)')
            header['Y_ASTIG'] = (self._dist[4], 'PSF RMS 45d astig (waves @ 547 nm)')
            header['X_COMA'] = (self._dist[5], 'PSF RMS X-coma (waves @ 547 nm)')
            header['Y_COMA'] = (self._dist[6], 'PSF RMS Y-coma (waves @ 547 nm)')
            header['X_CLOVER'] = (self._dist[7], 'PSF RMS X-clover (waves @ 547 nm)')
            header['Y_CLOVER'] = (self._dist[8], 'PSF RMS Y-clover (waves @ 547 nm)')
            header['SPHERICL'] = (self._dist[9], 'PSF RMS spherical (waves @ 547 nm)')
            header['PIXSCALE'] = (self.scale, 'Pixel scale in arcseconds')
            header['WAVELNTH'] = (self.wavelength, 'PSF wavelength in microns')

        if self._poly:
            header['SPECTYPE'] = (self.spectral_type.upper(), 'Spectral type of target star')
            header['BAND'] = (self.band.upper(), 'filter band used for polychromatic PSF')
            header['WAVELNTH'] = ('polychromatic', 'PSF wavelength in microns')
            for i, wavel in enumerate(self._wavelength_contributions[0]):
                header['WAVEL%s' % i] = (float('%.3f' % wavel), 'wavelength in microns')

        name = "{}.fits".format(name is not None and name or self.name)  # if name != None: name = name, else: self.name
        hdu.writeto(name, clobber=True)
        print >> sys.stdout, 'Saving %s' % name

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

    init takes: size
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

    def fits2aperture(self, input_file='aperture.fits'):
        """
        Opens an aperture fits file and reads it.
        """
        hdu = pyfits.open(input_file)
        if len(hdu) < 2:  # check if aperture.fits has a header (it should not)
            self.a = hdu[0].data
        else:
            self.a = hdu[1].data

        self.array_size = np.shape(self.a)[0]
        self.name = input_file.split('.')[0]
        logging.info("'%s' loaded, image size=%s" % (input_file, self.array_size))
        return


class Pupil(OpticalArray):
    """
    Pupil inherits from OpticalArray, it has properties *opd* (the optical path differences) and *name*.
    It can be computed using *_compute_pupil* which uses the opd, or *_path_diff* of the wavefront.

    init takes: wavelength, array_size
    """
    def __init__(self, wavelength, array_size):
        super(Pupil, self).__init__(array_size)
        self.wavelength = wavelength
        self.opd = self._path_diff()
        self._compute_pupil()
        self.name = 'Pupil'

    def _compute_pupil(self):
        logging.info("Computing pupil...")

        aperture = Aperture(self.array_size)
        if glob.glob('aperture.fits'):
            logging.info("... using 'aperture.fits'...")
            aperture.fits2aperture()
        else:
            aperture.make_hst_ap()

        self.a = np.multiply(aperture.a, np.exp(np.divide(2j*np.pi*self.opd, self.wavelength)))
        logging.info("... done")

    def _path_diff(self):
        zernike_modes = [(1, 1), (1, -1), (2, 0), (2, -2), (2, 2), (3, -1), (3, 1), (3, -3), (3, 3), (4, 0)]  # z2..z11

        rho_range = np.linspace(0, 1, self.array_size)
        theta_range = np.linspace(0, 2*np.pi, self.array_size)
        rho, theta = np.meshgrid(rho_range, theta_range)

        zernike_total = OpticalArray(polar=True, size=self.array_size)
        for i in range(len(self.dist)):
            aj = self.dist[i]*.547/self.wavelength  # Zernike coefficient in microns, .547um is the reference wavelength
            logging.info("Computing Z%s with aj=%s" % (2+i, aj))
            n, m = zernike_modes[i][0], zernike_modes[i][1]
            if m < 0.:
                zernike_value = ze.odd_zernike(n, -m, rho, theta)
            else:
                zernike_value = ze.even_zernike(n, m, rho, theta)
            zernike_total.a += np.multiply(zernike_value, aj)  # OPD = Sum aj Zj
        logging.info("OPD computed...")
        zernike_total.polar2cart()
        logging.info("... and converted back to cartesian space.")

        return zernike_total.a


class PSF(OpticalArray):
    """
    PSF inherits from OpticalArray and has properties *array_size*, *pupil*, *_a* the array containing the data, and
    *name*. Its resolution can be modified using *resize_psf*, which will change the pixel resolution to *scale*.

    init takes: pupil, scale, array_size
    """
    def __init__(self, pupil, scale, array_size):
        self.array_size = array_size*5  # for 0 padding
        self.pupil = pupil
        self._a = None
        super(PSF, self).__init__(self.array_size, pupil.wavelength, scale)
        self.name = 'Psf'
        self._dist = pupil.dist
        self.wavelength = pupil.wavelength
        self.scale = scale

    @property
    def a(self):
        if self._a is not None:
            return self._a
        tmp = np.fft.fft2(self.pupil.a, s=[self.array_size, self.array_size])  # padding with 0s
        tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array
        logging.info("Updating PSF ... done")
        self._a = np.real(np.multiply(tmp, np.conjugate(tmp)))
        return self._a

    @a.setter
    def a(self, a):
        self._a = a

    def resize_psf(self):
        """
        Re-scale the PSF to the desired pixel resolution (*self.scale*). When calling geometric_transform, the
        prefilter option has to be turned off (False) because the data is already filtered (no noise).
        """
        logging.debug("resize_psf parameters: wavelength=%s, array_size=%s, scale=%s" % (self.wavelength,
                                                                                         self.array_size, self.scale))
        logging.info("Resizing PSF to match pixel resolution of %s''/px..." % self.scale)
        logging.debug("before interpolation: max(psf.a)=%s, min(psf.a)=%s" % (np.max(self.a), np.min(self.a)))
        new_psf = scipy.ndimage.interpolation.geometric_transform(self.a, rebin, order=3, prefilter=False,
                                                                  extra_arguments=(self.wavelength, self.array_size,
                                                                                   self.scale))
        logging.debug("after interpolation: max(psf.a)=%s, min(psf.a)=%s" % (np.max(new_psf), np.min(new_psf)))
        logging.info("... done")
        self.a = new_psf


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


class PolyPSF(OpticalArray):
    """
    Polychromatic PSF class. Inherits from OpticalArray.

    init takes: band, spectral_type, size, scale=0.01
    """
    def __init__(self, band, spectral_type, size, scale=0.01):
        super(PolyPSF, self).__init__(size=size, poly=True, scale=scale)
        self.band = band.title()
        self.spectral_type = spectral_type.title()
        self._wavelength_contributions = None
        self.sed_file = 'sed_%s0_fe0.txt' % spectral_type.lower()
        self._x = None
        self._flux = None
        self.name = 'polychromatic_psf'
        self._dist = None
        self.b = []

    def get_sed(self):
        """
        Reads the SED profile from files and get the spectral energy distribution.
        :return: self._x, self._flux : two lists, wavelength and its associated flux.
        """
        if self._x is not None and self._flux is not None:
            return self._x, self._flux

        data = np.genfromtxt('SED/%s' % self.sed_file)
        x = data.transpose()[0]
        flux = data.transpose()[1]
        self._x = x
        self._flux = flux

    def wavelength_contributions(self):
        """
        Computes the contributions of all 10 evenly-spaced wavelength in the filter band.
        :return: wavelength_contributions. tuple of lists, containing the 10 wavelength and their relative contribution
        """
        if self._wavelength_contributions is not None:
            return self._wavelength_contributions

        bands = dict(Z=(0.76, 0.977), Y=(0.927, 1.192), J=(1.131, 1.454), H=(1.380, 1.774), F=(1.683, 2.000),
                     Wide=(0.927, 2.000))

        spectral_interpolation = interp1d(self._x, self._flux, kind='cubic')

        waves = np.linspace(bands[self.band][0], bands[self.band][1], 10)  # Takes 10 wavelengths in band
        self._wavelength_contributions = [waves, spectral_interpolation(waves)]

    def create_polychrome(self):
        """
        Creates a polychromatic PSF by adding 10 PSFs computed at 10 wavelengths from self._wavelength_contributions
        and add them to the list self.b
        """
        logging.debug("polychrome array_size=%s" % self.array_size)
        tmp = np.zeros((self.array_size*5, self.array_size*5))  # after FFT, the array will be 5*bigger -> 0-padding
        for i, wavel in enumerate(self._wavelength_contributions[0]):
            pupil = Pupil(wavel, self.array_size)
            logging.debug("pupil array_size=%s" % pupil.array_size)
            psf = PSF(pupil, self.scale, self.array_size)
            psf.resize_psf()  # scale is supposed to be 0.01
            if logging.debug:
                psf.save('polyPSF_%s' % wavel)  # will save all the 10 PSFs in separate files if debug mode is on
            logging.debug("psf.a size: %s" % psf.array_size)
            psf.a *= self._wavelength_contributions[1][i]
            self.b.append(psf.a)  # add the array to the list of arrays for data cube creation
            tmp += psf.a
        self._a = tmp

    def check_sed(self):
        """
        plot the SED of the star and the 10 wavelengths used to create the PSF, at the filter position.
        """
        import matplotlib.pyplot as plt
        plt.plot(self._x, self._flux, 'b')
        plt.plot(self._wavelength_contributions[0], self._wavelength_contributions[1], 'ro')
        plt.show()