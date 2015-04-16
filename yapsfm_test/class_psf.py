#! /usr/bin/env python
"""
Test class for PSF generation and handling
"""

import numpy as np
import scipy.ndimage.interpolation
import pyfits  # fits file handling and creation
from datetime import datetime as dt


class PSF(pyfits.ImageHDU):
    def __init__(self, size, wavelength, padding, pupil, aperture):
        super(PSF, self).__init__()
        self.a = np.zeros((size, size))  # actual array
        self.array_size = size  # size of array, the class already has a size method, giving size in bytes of the data
        self.center = [size/2, size/2]  # center of array
        self.wavelength = wavelength  # wavelength at which to compute PSF
        self.padding = padding  # zero-padding factor for FFT
        self.dist = dist  # Zernike coefficients
        self.aperture = aperture
        self.pupil = pupil

    def compute(self):
        padded = self.array_size*padding

        print 'Starting FFT with zero-padding factor of %s...' % scale_factor
        tmp = np.fft.fft2(self.pupil, s=[padded, padded])  # padding with 0s
        tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array
        print '... done'

        self.a = np.real(np.multiply(tmp, np.conjugate(tmp)))
        print '----------\nPSF image size: (%s,%s)' % (np.shape(psf)[0], np.shape(psf)[1])
        print 'lambda = %s' % wavelength
        print "pixel size @ 0.76um: 0.0133''/px"
        print '----------\nPSF computed'