#! /usr/bin/env python
"""
The common ancestor class between aperture, pupil and the PSF object is OpticalArray: a ImageHDU object with
properties *array_size*, *center*, *a* and *name*.
One can save the OpticalArray instance in a fits format using *save*.

Aperture inherits from OpticalArray and has methods *circular* to make a circular aperture or *make_hst_ap* to make an
HST-like aperture. Furthermore, one can *reset* the aperture to circular shape, or *fits2aperture* to load an aperture
from a .fits image.
"""

import glob
import numpy as np
import pyfits
import scipy.ndimage.interpolation


class OpticalArray(object):
    def __init__(self, size):
        self.array_size = size
        self.center = [self.array_size//2., self.array_size//2.]
        self.a = np.zeros((self.array_size, self.array_size))
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


class Distorted(OpticalArray):
    """
    The ancestor class of both *Pupil* and *PSF* objects.
    Based on *OpticalArray*, but have extra attributes *wavelength* and optical distortions *dist*
    """
    def __init__(self, wavelength, dist):
        super(OpticalArray, self).__init__()
        self.wavelength = wavelength
        self.dist = dist


class Pupil(Distorted):
    def __init__(self, wavelength):
        # self.wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns) ') or .76)
        self.wavelength = wavelength
        self._read_distortions()
        self.array_size = 101
        self.opd = self._path_diff()
        super(Distorted, self).__init__(self.array_size)
        self._compute_pupil()
        self.name = 'Pupil'

    def _read_distortions(self, file_path='distortions.par'):
        data = open(file_path).readlines()

        dist = []
        for i, line in enumerate(data):
            dist.append(float(data[i].partition('#')[0].strip()))
        self.dist = dist

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

        zernike_total = np.zeros((self.array_size, self.array_size))
        for i in range(len(self.dist)):
            aj = self.dist[i]*.547/self.wavelength  # Zernike coefficient in microns, .547um is the reference wavelength
            print 'Computing Z%s with aj=%s' % (4+i, aj)
            n, m = zernike_modes[i][0], zernike_modes[i][1]
            if m < 0.:
                zernike_value = odd_zernike(n, -m, rho, theta)
            else:
                zernike_value = even_zernike(n, m, rho, theta)
            zernike_total += np.multiply(zernike_value, aj)  # OPD = Sum aj Zj
        print 'OPD computed...'
        opd = scipy.ndimage.interpolation.geometric_transform(zernike_total, polar2cart, extra_arguments=(
            self.array_size,))
        print '... and converted back to cartesian space.'
        return opd


def r_nm(n, m, r):
    """
    Computes the radial part R_n^m(r), with r a meshgrid object

    :param n: radial order of the polynomial
    :type n: int
    :param m: azimuthal order of the polynomial
    :type m: int
    :param r: radial position
    :type r: float
    :return: radial part of the Zernike mode at position r
    :rtype: float
    """

    r_nm_value = 0
    for s in range((n-m)/2+1):
        r_nm_value += (((-1)**s * np.math.factorial(n-s)) / (np.math.factorial(s) * np.math.factorial((n+m)/2-s) *
                                                             np.math.factorial((n-m)/2-s))) * r**(n-2*s)
    return r_nm_value


def even_zernike(n, m, r, theta):
    """
    Computes the even Zernike mode following Noll (1976)

    :param n: radial order of the mode
    :type n: int
    :param m: azimuthal order of the mode
    :type m: int
    :param r: radial position
    :type r: float
    :param theta: angle
    :type theta: float
    :return: the value of the Zernike mode corresponding to n and m, at position (r,theta).
    :rtype: float
    """

    zernike_value = np.sqrt(n+1)*r_nm(n, m, r)*np.sqrt(2)*np.cos(m*theta)
    return zernike_value


def odd_zernike(n, m, r, theta):
    """
    Computes the odd Zernike mode following Noll (1976)

    :param n: radial order of the mode
    :type n: int
    :param m: azimuthal order of the mode
    :type m: int
    :param r: radial position
    :type r: float
    :param theta: angle
    :type theta: float
    :return: the value of the Zernike mode corresponding to n and m, at position (r,theta).
    :rtype: float
    """

    zernike_value = np.sqrt(n+1)*r_nm(n, m, r)*np.sqrt(2)*np.sin(m*theta)
    return zernike_value


def polar2cart(coords, size=101):
    """
    Change between polar and cartesian coordinates system

    :param coords: list of array position
    :type coords: np.array
    :param size: total size of the image
    :type size: int
    :return: the (theta,r) values in (x,y)
    """

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


class PSF(Distorted):
    def __init__(self, pupil):
        self.array_size = 505
        super(Distorted, self).__init__(self.array_size)
        self._compute_psf(pupil)
        self.name = 'Psf'

    def _compute_psf(self, pupil):
        print 'Starting FFT with zero-padding factor of %s...' % (self.array_size/np.shape(pupil.a)[0])
        tmp = np.fft.fft2(pupil.a, s=[self.array_size, self.array_size])  # padding with 0s
        tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array
        print "... done"
        self.a = np.real(np.multiply(tmp, np.conjugate(tmp)))


def main():
    wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns) ') or .76)
    pupil = Pupil(wavelength)

    psf = PSF(pupil)
    psf.save()  # does not fill the header of the .fits image

if __name__ == "__main__":
    main()