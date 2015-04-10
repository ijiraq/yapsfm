import unittest

from PSFModel import odd_zernike


class PSFModelTest(unittest.TestCase):

    def test_zodd(self):
        """
        Test the Zodd function as it's the only one I can reasonably
        understand in less that a few seconds of looking at it.

        No idea what it actually does, of course. Let's not let that get in the
        way.
        """
        result = odd_zernike(1, 1, 1, 1)
        self.assertEqual(1.6829419696157935, result)  # Woo I'm doing astronomy

    # You can add further tests by simply adding methods starting with test_*
    # and importing your test target at the top of this file.
