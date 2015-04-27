#! /usr/bin/env python
"""
Distorted inherits from OpticalArray and has properties *wavelength* and *_dist*, corresponding to the optical
distortion affecting the psf.
"""

from optiarray import OpticalArray


class Distorted(OpticalArray):
    """
    The ancestor class of both *Pupil* and *PSF* objects.
    Inherits from *OpticalArray*, but have extra attributes *wavelength* and optical distortions *dist*
    """
    def __init__(self, size, wavelength):
        super(Distorted, self).__init__(size)
        self.wavelength = wavelength
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