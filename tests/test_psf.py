import unittest

import glob
import numpy as np
from PSFModel import aperture
from PSFModel import odd_zernike
from PSFModel import fits2aperture

class PSFModelTest(unittest.TestCase):

    def test_aperture(self):
        if glob.glob('aperture.fits'):
            result = fits2aperture('aperture.fits')
        else:
            result = aperture(101, 'HST')
        self.assertEqual(0, result[0, 0])
        self.assertEqual(0, result[50, 0])
        self.assertEqual(0, result[50, 50])
        self.assertEqual(0, result[100, 100])

    def test_zodd(self):
        """
        Test the Zodd function as it's the only one I can reasonably
        understand in less that a few seconds of looking at it.

        No idea what it actually does, of course. Let's not let that get in the
        way.
        """
        result = odd_zernike(1, 1, 1, np.pi)
        self.assertEqual(1.6829419696157935, result)

    # You can add further tests by simply adding methods starting with test_*
    # and importing your test target at the top of this file.
