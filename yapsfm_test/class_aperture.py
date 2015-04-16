#! /usr/bin/env python
"""
Test classes for aperture generation using an Aperture class

When initialized, the object is a simple circular aperture of size *size* and center *center*.
It can then be changed into an HST-like aperture using .make_hst_ap()
or a file 'aperture.fits' can be loaded using .fits2aperture()
finally, the aperture can be saved into a .fits file using .save()
"""

import numpy as np
import pyfits


class Aperture(object):
    def __init__(self, size):
        self.size = size
        self.center = [size/2, size/2]

        self.a = np.zeros((self.size, self.size))
        for y in range(self.size):
            for x in range(self.size):
                radial_position = np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
                if radial_position <= self.size/2:
                    self.a[y][x] = 1.

    def make_hst_ap(self):
        """
        Makes an HST-like aperture of size *self.size*
        """
        sec_mir = 0.330*self.size/2  # secondary mirror
        sp_width = 0.022*self.size/2  # spider width
        pad1 = [0.8921*self.size/2, 0.0000, 0.065*self.size/2]  # x, y, radius
        pad2 = [-0.4615*self.size/2, 0.7555*self.size/2, 0.065*self.size/2]
        pad3 = [-0.4564*self.size/2, -0.7606*self.size/2, 0.065*self.size/2]
        for y in range(self.size):
            for x in range(self.size):
                # main aperture (including secondary obstruction)
                radial_position = np.sqrt((x - self.center[0])**2+(y - self.center[1])**2)
                if self.size/2 >= radial_position > sec_mir:
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

    def reset(self):
        """
        Reset the aperture to simple circular.
        """
        for y in range(self.size):
            for x in range(self.size):
                radial_position = np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
                if radial_position <= self.size/2:
                    self.a[y][x] = 1.

    def fits2aperture(self):
        """
        Opens an aperture fits file and reads it.
        """
        hdu = pyfits.open('aperture.fits')
        if len(hdu) < 2:  # check if aperture.fits has a header (it should not)
            self.a = hdu[0].data
        else:
            self.a = hdu[1].data

        self.size = np.shape(self.a)[0]

    def save(self):
        hdu = pyfits.PrimaryHDU(self.a)
        hdu.writeto('aperture.fits', clobber=True)
        print "aperture saved to 'aperture.fits'"