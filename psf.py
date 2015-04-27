#! /usr/bin/env python
"""
PSF inherits from Distorted and has properties *array_size*, *pupil*, *_a* and *name*. Its method is *rezie_psf* to
a desired resolution.
"""

import numpy as np
import scipy.ndimage.interpolation

from distorted import Distorted


class PSF(Distorted):
    """
    PSF inherits from Distorted and has properties *array_size*, *pupil*, *_a* the array containing the data, and
    *name*. Its resolution can be modified using *resize_psf*, which will change the pixel resolution to *scale*.
    """
    def __init__(self, pupil):
        self.array_size = 505
        self.pupil = pupil
        self._a = None
        super(Distorted, self).__init__(self.array_size, pupil.wavelength)
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