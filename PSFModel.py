#! /usr/bin/env python
"""
The script can load a gray scale fits image with fits2aperture(.fits)
or create a circular aperture using aperture(size)

Feb 25th: PSF computed as perfect theoretical Airy Disk: OPD=0.

Mar 3rd: computing the OPD in cartesian space is problematic because of the polar nature of Zernike modes.
Compute it in polar coordinates, interpolate to a continuous function and map into cartesian space.

Mar 6th: OPD computation is working. Linear combination of Zernike mode done.

Mar 9th: PSF scaling to detector space (0.005''/pixel).

Mar 10th: Code cleaned and polar2cart: theta fixed.

Mar 11th: .fits creation with header information corresponding to Tiny Tim's. Possibility to change pixel scale,
defaulted as constant with wavelength.

Apr 9th: PEP8 code clean and docstrings creation

Apr 10th: Aperture creation in .fits file and .fits aperture input handling added.
"""

import glob
import numpy as np
import scipy.ndimage.interpolation  # image handling
import pyfits  # fits file handling and creation
from datetime import datetime as dt

# --------------------------------------------------
# Methods definition
# --------------------------------------------------


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
    # hdulist = pyfits.HDUList([hdu])
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


def psf_compute(a, wavelength=.76, scale_factor=5, dist=[0, 0, 0, 0, 0, 0, 0, 0]):
    """
    Compute the PSF at wavelength L, using aperture A and optical distortions dist
    corresponding to Zernike modes 4 to 11 (defocus to spherical)

    :param a: aperture array of 0s and 1s
    :type a: np.array
    :param wavelength: wavelength at which to compute the PSF, in microns
    :type wavelength: float
    :param scale_factor: zero-padding factor for FFT. It adds 2x the size of the array in 0s on all sides of the
    initial array
    :type scale_factor: int
    :param dist: Zernike coefficient values for Z4 (defocus) to Z11 (spherical). Given in values of rms waves at .547um.
    :type dist: np.array
    :return: PSF image
    :rtype: np.array
    """

    """fft=complex numbers: amplitude AND phase
    We take the modulus square --> distribution of light
    L is the wavelength, same units as OPD (microns)
    np.fft.fft2 manages zero-padding on each axis with s[0],s[1] corresponding to both x and y axis.
    """

    pup = pupil(a, wavelength, dist)
    size = np.shape(pup)
    scaled = [size[i]*scale_factor for i in range(len(size))]
    print 'Starting FFT with zero-padding factor of %s...' % scale_factor
    tmp = np.fft.fft2(pup, s=[scaled[0], scaled[1]])  # padding with 0s

    tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array

    print '... done'
    PSF = np.real(np.multiply(tmp, np.conjugate(tmp)))
    print '----------\nPSF image size: (%s,%s)' % (np.shape(PSF)[0], np.shape(PSF)[1])
    print 'lambda = %s' % wavelength
    print "pixel size @ 0.76um: 0.0133''/px"

    print '----------\nPSF computed'
    return PSF

# --------------------------------------------------


def pupil(a, wavelength=.76, dist=[0, 0, 0, 0, 0, 0, 0, 0]):
    """
    Compute the pupil function for aperture A, wavelength L and optical distortions dist

    :param a: aperture function
    :type a: np.array
    :param wavelength: wavelength in microns
    :type wavelength: float
    :param dist: Zernike coefficient values for Z4 (defocus) to Z11 (spherical)
    :type dist: np.array
    :return: the pupil function evaluated at each point of an array
    :rtype: np.array
    """

    """P = A exp(2pi i OPD / L), L=lambda"""

    print 'Computing pupil...'
    size = np.shape(a)[0]
    opd = path_diff(size, wavelength, dist)  # the optical path difference (OPD)
    pup = np.multiply(a, np.exp(np.divide(2j*np.pi*opd, wavelength)))
    print '... done'
    return pup

# --------------------------------------------------


def path_diff(size=101, wavelength=.76, dist=[0., 0., 0., 0., 0., 0., 0., 0.]):
    """
    Computes the optical path differences at wavelength L for optical distortions dist.

    :param size: the array size on which to compute the OPD
    :type size: int
    :param wavelength: wavelength in microns
    :type wavelength: float
    :param dist: Zernike coefficient values for Z4 (defocus) to Z11 (spherical)
    :type dist: np.array
    :return: an array of OPD
    :rtype: np.array
    """

    """
    ==================================================
    Optical Path Differences for pupil characterization
    ==================================================
    from Noll (1976):
    If phi(r,theta) is an arbitrary function, its expansion over a circle of radius R is:
    phi(R*rho, theta) = Sum_j a_j Z_j(rho,theta), rho=r/R
    --------------------
    OPD(R*rho, theta) is therefore equals to phi(R*rho, theta)
    which is, in cartesian space: OPD(sqrt(x**2+y**2), atan(y/x))
    the wavelength dependence is hidden in a_j. -> each Zernike mode
    has a different coefficient, depending on the wavelength.
    ==================================================
    Method:
    Compute the Zernike mode(s) and multiply each of them by its coefficient.
    The matrix of the Zernike values is in polar coordinates. Has to be
    transformed back to cartesian.
    ==================================================
    1 micron of despace induces 6.1nm of rms wave
    Zernike coefficient values are given in microns RMS of wave at 547 microns
    """

    zernike_modes = [(2, 0), (2, -2), (2, 2), (3, -1), (3, 1), (3, -3), (3, 3), (4, 0)]  # defocus,z5,z6,z7,z8,z9,z10,z11

    rho_range = np.linspace(0, 1, size)
    theta_range = np.linspace(0, 2*np.pi, size)
    rho, theta = np.meshgrid(rho_range, theta_range)
    
    zernike_total = np.zeros((size, size))
    for i in range(len(dist)):
        aj = dist[i]*.547/wavelength  # Zernike coefficient at wavelength in microns
        print 'Computing Z%s with aj=%s' % (4+i, aj)
        n, m = zernike_modes[i][0], zernike_modes[i][1]
        if m < 0.:
            zernike_value = odd_zernike(n, -m, rho, theta)
        else:
            zernike_value = even_zernike(n, m, rho, theta)
        zernike_total += np.multiply(zernike_value, aj)  # OPD = Sum aj Zj
    print 'OPD computed...'
    # print 'Saving zernike_total'
    # plt.imshow(zernike_total)
    # plt.savefig('zernike_total.png')

    cartesian_z = scipy.ndimage.interpolation.geometric_transform(zernike_total, polar2cart, extra_arguments=(size,))
    print '... and converted back to cartesian space.'
    # print 'Saving cartesian_Z'
    # plt.imshow(cartesian_z)
    # plt.savefig('cartesian_z.png')

    return cartesian_z

# --------------------------------------------------


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
        r_nm_value += (((-1)**s*np.math.factorial(n-s))/(np.math.factorial(s)*np.math.factorial((n+m)/2-s)*np.math.factorial((n-m)/2-s)))*r**(n-2*s)
    return r_nm_value

# --------------------------------------------------


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

# --------------------------------------------------


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

# --------------------------------------------------


def polar2cart(coords, size=101):
    """
    Change between polar and cartesian coordinates system

    :param coords: list of array position
    :type coords: array
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

# --------------------------------------------------


def jitter(psf, jitter_size):
    """
    Computes the optical transfer function (OTF) and multiply it by the gaussian jitter function

    Work in progress: currently just a skeleton

    :param psf:
    :param jitter_size:
    :return:
    """

    """WORK IN PROGRESS"""
    jitter = 0.
    OTF = np.fft.fft2(psf)*jitter
    OTF = np.fft.ifft2(OTF)
    return

# --------------------------------------------------


def bin2detector(coords, wavelength, size, detector_scale):
    """
    re-bin the Optical Transfer Function to the detector's scale.

    Currently not used, since the desired PSF has to be oversampled.

    :param coords: list of (x,y) positions
    :type coords: np.array
    :param wavelength: wavelength in microns
    :type wavelength: float
    :param size: size of the image
    :type size: int
    :param detector_scale:
    :type detector_scale: int
    :return: the positions re-binned to detector space
    :rtype:
    """

    """
    Used in resizePSF()
    The OTF can be the PSF, or if jitter is specified, its convolution with a gaussian.
    """

    # scale the PSF to the desired size (0.106 arcsec)
    scale_factor = 0.0797/6.*(wavelength/0.76)  # at 5x zero-padding

    x = (coords[1]-size//2.)*detector_scale/scale_factor+(size//2.)
    y = (coords[0]-size//2.)*detector_scale/scale_factor+(size//2.)
    return y, x

# --------------------------------------------------


def resize_psf(psf, wavelength=.76, size=505, scale=0.110):
    """
    Resize the PSF to match pixel size and resolution of instrument (0.12'' at .76um)

    :param psf: the PSF array
    :type psf: np.array
    :param wavelength: wavelength in microns
    :type wavelength: float
    :param size: size of the image
    :type size: int
    :param scale: detector pixel size in arcsec
    :type scale: float
    :return: re-binned PSF
    :rtype: np.array
    """

    print "resizing PSF to match detector pixel size of %s''/px..." % scale
    new_psf = scipy.ndimage.interpolation.geometric_transform(
        psf, bin2detector, extra_arguments=(wavelength, size, scale))
    new_psf = new_psf[size//2.-32:size//2.+32, size//2.-32:size//2.+32]
    print '... done'
    return new_psf

# --------------------------------------------------


def create_fits(psf, dist=[0, 0, 0, 0, 0, 0, 0, 0], wavelength=0.76):
    """
    Creates a .fits file containing the PSF image with header information: Created, Instrument,
    Focus, Astigmatism (0,45), Coma (x,y), Trefoil (x,y), Spherical, Pixel scale and Wavelength

    :param psf: PSF array to convert to .fits
    :type psf: np.array
    :param dist: the optical distortion coefficients
    :type dist: array
    :param wavelength: wavelength in microns
    :type wavelength: float
    :return: None
    """

    name = 'psf_%s.fits' % wavelength
    print 'Writing psf to file %s...' % name

    hdu = pyfits.PrimaryHDU(psf)
    header = hdu.header

    now = dt.now()
    header['CREATED'] = ('%s %s %s %s %s' % (dt.strftime(now, '%a'), dt.strftime(now, '%b'), dt.strftime(now, '%d'),
                                             dt.strftime(now, '%X'), dt.strftime(now, '%Y')),
                         'Time and date file was created')
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

    hdu.writeto('%s' % name, clobber=True)
    print '... done'
    return

# ==================================================


def main():
    wavelength = float(raw_input('Lambda? (0.76-2.00 microns) ') or .76)

    if glob.glob('aperture.fits'):
        ap = fits2aperture('aperture.fits')
    else:
        ap = aperture(101, 'HST')

    # Zernike coefficients for distortions Z4 to Z11 (defocus to spherical)
    dist = [0.0026, 0.0089, 0.0222, -0.0018, 0.0113, 0., 0., 0.]
    
    psf = psf_compute(ap, wavelength, 5, dist)

    pixel_scale = 0.110  # of the instrument, in arcseconds per pixel

    # newPSF = resizePSF(psf, wavelength, size=np.shape(psf)[0], pixel_scale)  # re-bin the psf to detector size.

    create_fits(psf, wavelength=wavelength, dist=dist)
    return

# ==================================================

if __name__ == '__main__':
    main()
