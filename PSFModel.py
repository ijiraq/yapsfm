#! /usr/bin/env python
"""
PSF modeling script:
the pupil function P(x,y) = A(x,y) exp(2 \pi i OPD(x,y) / \lambda)
the PSF = | FFT(P(x,y)) |**2
+ convolution with a gaussian for jitter and imprecision
--------------------
The script can load a gray scale image with imageToAperture(image)
or create a circular aperture using aperture(size)

Feb 25th: PSF computed as perfect theoretical Airy Disk: OPD=0.

Mar 3rd: computing the OPD in cartesian space is problematic because of the polar nature of Zernike modes.
Can't get a good representation in cartesian space right away. Decided to compute it in polar coordinates,
interpolate it to a continuous function and map it into cartesian space.

Mar 6th: OPD computation is working. Linear combination of Zernike mode done.

Mar 9th: PSF scaling to detector space (0.005''/pixel).

Mar 10th: Code cleaned and polar2cart: theta fixed.

Mar 11th: .fits creation with header information corresponding to Tiny Tim's. Possibility to change pixel scale,
defaulted as constant with wavelength.

Apr 9th: PEP8 code clean and docstrings creation
"""

import numpy as np
import png
import matplotlib.pyplot as plt
from scipy.misc import imread
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
    :return A: array containing 0s and 1s representing obstruction of light
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


    A = np.zeros((size, size))

    center = [size/2, size/2]  # center of the image
    secMir = 0.330*size/2  # secondary mirror
    spWidth = 0.022*size/2  # spider width
    pad1 = [0.8921*size/2, 0.0000, 0.065*size/2]  # x, y, radius
    pad2 = [-0.4615*size/2, 0.7555*size/2, 0.065*size/2]
    pad3 = [-0.4564*size/2, -0.7606*size/2, 0.065*size/2]
    
    for y in range(size):
        for x in range(size):
            # main aperture (including secondary obstruction)
            radPos = np.sqrt((x-center[0])**2+(y-center[1])**2)
            if ap == 'HST':
                if size/2 >= radPos > secMir:
                    A[y][x] = 1.
                    # Obstructions:
                    # spiders                
                    if center[0]-spWidth/2. <= x <= center[0]+spWidth/2:
                        A[y][x] = 0.
                    if center[0]-spWidth/2 <= y <= center[1]+spWidth/2:
                        A[y][x] = 0.
                    # pads
                    if np.sqrt((x-center[0]-pad1[0])**2+(y-center[1]-pad1[1])**2) <= pad1[2]:
                        A[y][x] = 0.
                    if np.sqrt((x-center[0]-pad2[0])**2+(y-center[1]-pad2[1])**2) <= pad2[2]:
                        A[y][x] = 0.
                    if np.sqrt((x-center[0]-pad3[0])**2+(y-center[1]-pad3[1])**2) <= pad3[2]:
                        A[y][x] = 0.
            else:
                if radPos <= size/2:
                    A[y][x] = 1.

    print 'Aperture image size: (%s,%s)' % (len(A), len(A[0]))
    png.from_array(A, mode='L;1').save('analyticAp.png')
    print 'Aperture created'

    return A

# --------------------------------------------------


def psf(A, L=.76, scaleFactor=5, dist=[0, 0, 0, 0, 0, 0, 0, 0]):
    """
    Compute the PSF at wavelength L, using aperture A and optical distortions dist
    corresponding to Zernike modes 4 to 11 (defocus to spherical)

    :param A: aperture array of 0s and 1s
    :type A: np.array
    :param L: wavelength at which to compute the PSF, in microns
    :type L: float
    :param scaleFactor: zero-padding factor for FFT. It adds 2x the size of the array in 0s on all sides of the
    initial array
    :type scaleFactor: int
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

    P = pupil(A, L, dist)
    size = np.shape(P)
    scaled = [size[i]*scaleFactor for i in range(len(size))]
    print 'Starting FFT with zero-padding factor of %s...' % (scaleFactor)
    tmp = np.fft.fft2(P, s=[scaled[0], scaled[1]])  # padding with 0s

    tmp = np.fft.fftshift(tmp)  # switch quadrant to place the origin in the middle of the array

    print '... done'
    PSF = np.real(np.multiply(tmp, np.conjugate(tmp)))
    print '----------\nPSF image size: (%s,%s)' % (np.shape(PSF)[0], np.shape(PSF)[1])
    print 'lambda = %s' % (L)
    print "pixel size @ 0.76um: 0.0133''/px"

    print '----------\nPSF computed'
    return PSF

# --------------------------------------------------


def pupil(A, L=.76, dist=[0, 0, 0, 0, 0, 0, 0, 0]):
    """
    Compute the pupil function for aperture A, wavelength L and optical distortions dist

    :param A: aperture function
    :type A: np.array
    :param L: wavelength in microns
    :type L: float
    :param dist: Zernike coefficient values for Z4 (defocus) to Z11 (spherical)
    :type param: np.array
    :return: the pupil function evaluated at each point of an array
    :rtype: np.array
    """

    """P = A exp(2pi i OPD / L), L=lambda"""

    print 'Computing pupil...'
    size = np.shape(A)[0]
    OPD = pathDiff(size, L, dist)  # the optical path difference (OPD)
    P = np.multiply(A, np.exp(np.divide(2j*np.pi*OPD, L)))
    print '... done'
    return P

# --------------------------------------------------


def pathDiff(size=101, L=.76, dist=[0., 0., 0., 0., 0., 0., 0., 0.]):
    """
    Computes the optical path differences at wavelength L for optical distortions dist.

    :param size: the array size on which to compute the OPD
    :type size: int
    :param L: wvelength in microns
    :type L: float
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

    Znm = [(2, 0), (2, -2), (2, 2), (3, -1), (3, 1), (3, -3), (3, 3), (4, 0)]  # defocus,z5,z6,z7,z8,z9,z10,z11

    rhorange = np.linspace(0, 1, size)
    thetarange = np.linspace(0, 2*np.pi, size)
    rho, theta = np.meshgrid(rhorange, thetarange)
    
    Ztot = np.zeros((size, size))
    for i in range(len(dist)):
        aj = dist[i]*.547/L  # Zernike coefficient at wavelength L. L and .547 in microns
        print 'Computing Z%s with aj=%s' % (4+i, aj)
        n, m = Znm[i][0], Znm[i][1]
        if m < 0.:
            Z = Zodd(n, -m, rho, theta)
        else:
            Z = Zeven(n, m, rho, theta)
        Ztot += np.multiply(Z, aj)  # OPD = Sum aj Zj
    print 'OPD computed...'
    # print 'Saving Ztot'
    # plt.imshow(Ztot)
    # plt.savefig('Ztot.png')

    cartesian_Z = scipy.ndimage.interpolation.geometric_transform(Ztot, polar2cart, extra_arguments=(size,))
    print '... and converted back to cartesian space.'
    # print 'Saving cartesian_Z'
    # plt.imshow(cartesian_Z)
    # plt.savefig('cartesian_Z.png')

    return cartesian_Z

# --------------------------------------------------


def Rnm(n, m, r):
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

    R = 0
    for s in range((n-m)/2+1):
        R += (((-1)**s*np.math.factorial(n-s))/(np.math.factorial(s)*np.math.factorial((n+m)/2-s)*np.math.factorial((n-m)/2-s)))*r**(n-2*s)
    return R 

# --------------------------------------------------


def Zeven(n, m, r, theta):
    """
    Computes the even Zernike mode

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

    Z = np.sqrt(n+1)*Rnm(n, m, r)*np.sqrt(2)*np.cos(m*theta)
    return Z

# --------------------------------------------------


def Zodd(n, m, r, theta):
    """
    Computes the odd Zernike mode

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

    Z = np.sqrt(n+1)*Rnm(n, m, r)*np.sqrt(2)*np.sin(m*theta)
    return Z

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


def jitter(PSF, jitterSize):
    """
    Computes the optical transfer function (OTF) and multiply it by the gaussian jitter function

    Work in progress: currently just a skeleton

    :param PSF:
    :param jitterSize:
    :return:
    """

    """WORK IN PROGRESS"""
    jitter = 0.
    OTF = np.fft.fft2(PSF)*jitter
    OTF = np.fft.ifft2(OTF)
    return

# --------------------------------------------------


def bin2detector(coords, L, size, detectorScale):
    """
    re-bin the Optical Transfer Function to the detector's scale.

    Currently not used, since the desired PSF has to be oversampled.

    :param coords: list of (x,y) positions
    :type coords: np.array
    :param L: wavelength in microns
    :type L: float
    :param size: size of the image
    :type size: int
    :param detectorScale:
    :type detectorScale: int
    :return: the positions re-binned to detector space
    :rtype:
    """

    """
    Used in resizePSF()
    The OTF can be the PSF, or if jitter is specified, its convolution with a gaussian.
    """

    # scale the PSF to the desired size (0.106 arcsec)
    scaleFactor = 0.0797/6.*(L/0.76)  # at 5x zero-padding

    x = (coords[1]-size//2.)*detectorScale/scaleFactor+(size//2.)
    y = (coords[0]-size//2.)*detectorScale/scaleFactor+(size//2.)
    return y, x

# --------------------------------------------------


def resizePSF(PSF, L=.76, size=505, scale=0.110):
    """
    Resize the PSF to match pixel size and resolution of instrument (0.12'' at .76um)

    :param PSF: the PSF array
    :type PSF: np.array
    :param L: wavelength in microns
    :type L: float
    :param size: size of the image
    :type size: int
    :param scale: detector pixel size in arcsec
    :type scale: float
    :return: re-binned PSF
    :rtype: np.array
    """

    print "resizing PSF to match detector pixel size of %s''/px..." % scale
    newPSF = scipy.ndimage.interpolation.geometric_transform(PSF, bin2detector, extra_arguments=(L, size, scale))
    newPSF = newPSF[size//2.-32:size//2.+32, size//2.-32:size//2.+32]
    print '... done'
    return newPSF

# --------------------------------------------------


def createFits(PSF, disto=[0, 0, 0, 0, 0, 0, 0, 0], wavelength=0.76):
    """
    Creates a .fits file containing the PSF image with header information: Created, Instrument,
    Focus, Astigmatism (0,45), Coma (x,y), Trefoil (x,y), Spherical, Pixel scale and Wavelength

    :param PSF: PSF array to convert to .fits
    :type PSF: np.array
    :param disto: the optical distortion coefficients
    :param disto: array
    :param wavelength: wavelength in microns
    :type wavelength: float
    :return: None
    """

    name = 'psf_%s.fits' % wavelength
    print 'Writing psf to file %s...' % name

    hdu = pyfits.PrimaryHDU(PSF)
    header = hdu.header

    now = dt.now()
    header['CREATED'] = ('%s %s %s %s %s' % (dt.strftime(now, '%a'), dt.strftime(now, '%b'), dt.strftime(now, '%d'),
                                             dt.strftime(now, '%X'), dt.strftime(now, '%Y')),
                         'Time and date file was created')
    header['INSTRUME'] = ('WFIRST_WFI', 'Simulated instrument')
    header['FOCUS'] = (disto[0], 'PSF RMS focus (waves @ 547 nm)')
    header['X_ASTIG'] = (disto[1], 'PSF RMS 0d astig (waves @ 547 nm)')
    header['Y_ASTIG'] = (disto[2], 'PSF RMS 45d astig (waves @ 547 nm)')
    header['X_COMA'] = (disto[3], 'PSF RMS X-coma (waves @ 547 nm)')
    header['Y_COMA'] = (disto[4], 'PSF RMS Y-coma (waves @ 547 nm)')
    header['X_CLOVER'] = (disto[5], 'PSF RMS X-clover (waves @ 547 nm)')
    header['Y_CLOVER'] = (disto[6], 'PSF RMS Y-clover (waves @ 547 nm)')
    header['SPHEICL'] = (disto[7], 'PSF RMS spherical (waves @ 547 nm)')
    header['PIXSCALE'] = ('work in progress', 'Pixel scale in arcseconds')
    header['WAVELNTH'] = (wavelength, 'PSF wavelength in microns')

    hdu.writeto('%s' % name, clobber=True)
    print '... done'
    return

# ==================================================


def main():
    L = float(raw_input('Lambda? (0.76-2.00 microns) ') or .76)
    A = aperture(101, 'HST')

    # Zernike coefficients for distortions Z4 to Z11 (defocus to spherical)
    dist = [0.0026, 0.0089, 0.0222, -0.0018, 0.0113, 0., 0., 0.]
    
    PSF = psf(A, L, 5, dist)

    pixelScale = 0.110  # of the instrument

    # newPSF = resizePSF(PSF, L, size=np.shape(PSF)[0], pixelScale)  # if one wants to re-bin the PSF to detector size.

    createFits(PSF, wavelength=L, disto=dist)
    return

# ==================================================

if __name__ == '__main__':
    main()
