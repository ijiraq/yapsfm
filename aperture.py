#! /usr/bin/env python
"""
Aperture handling:
aperture() creates a HST-like aperture and saves it as a .fits file
fits2aperture() opens an already existing .fits "aperture.fits" file and reads it
"""

import numpy as np
import pyfits


def aperture(size=101, ap=None):
    """
    Defines the aperture as an array containing 0s (obstructed) and 1s (unobstructed).

    This method is deprecated and will be replaced with a .fits image aperture input.

    :param size: size of the array defining the aperture.
    :type size: int, taken as odd to have a middle.
    :param ap: name of the telescope to model.
    :type ap: str, currently only "HST" or None is supported. If None, a simple circular aperture is used.
    :return a: array containing 0s and 1s representing aperture (light obstruction)
    :rtype: np.array
    """
    """
        for HST-like aperture:
    secondary mirror obstruction: 0.330
    spider width: 0.022
    pads:(v3(x), v2(y), radius) x, y are arbitrarily chosen here
    1: 0.8921  0.0000  0.065
    2: -0.4615 0.7555  0.065
    3: -0.4564 -0.7606 0.065
    ----------
    if nothing specified: circular aperture is used.
    """

    a = np.zeros((size, size))

    center = [size/2, size/2]  # center of the image
    sec_mir = 0.330*size/2  # secondary mirror
    sp_width = 0.022*size/2  # spider width
    pad1 = [0.8921*size/2, 0.0000, 0.065*size/2]  # x, y, radius
    pad2 = [-0.4615*size/2, 0.7555*size/2, 0.065*size/2]
    pad3 = [-0.4564*size/2, -0.7606*size/2, 0.065*size/2]

    for y in range(size):
        for x in range(size):
            # main aperture (including secondary obstruction)
            radial_position = np.sqrt((x-center[0])**2+(y-center[1])**2)
            if ap == 'HST':
                if size/2 >= radial_position > sec_mir:
                    a[y][x] = 1.
                    # Obstructions:
                    # spiders
                    if center[0]-sp_width/2. <= x <= center[0]+sp_width/2:
                        a[y][x] = 0.
                    if center[0]-sp_width/2 <= y <= center[1]+sp_width/2:
                        a[y][x] = 0.
                    # pads
                    if np.sqrt((x-center[0]-pad1[0])**2+(y-center[1]-pad1[1])**2) <= pad1[2]:
                        a[y][x] = 0.
                    if np.sqrt((x-center[0]-pad2[0])**2+(y-center[1]-pad2[1])**2) <= pad2[2]:
                        a[y][x] = 0.
                    if np.sqrt((x-center[0]-pad3[0])**2+(y-center[1]-pad3[1])**2) <= pad3[2]:
                        a[y][x] = 0.
            else:
                if radial_position <= size/2:
                    a[y][x] = 1.

    print 'Aperture image size: (%s,%s)' % (len(a), len(a[0]))
    hdu = pyfits.PrimaryHDU(a)
    hdu.writeto('aperture.fits', clobber=True)
    print 'Aperture created'

    return a

# --------------------------------------------------


def fits2aperture(input_file):
    """
    Opens an aperture fits file and reads it.

    :param input_file: the aperture .fits file path
    :type input_file: str
    :return: an aperture array
    :rtype: np.array
    """
    hdu = pyfits.open(input_file)
    if len(hdu) < 2:  # check if aperture.fits has a header (it should not)
        a = hdu[0].data
    else:
        a = hdu[1].data

    return a

# --------------------------------------------------