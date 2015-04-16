#! /usr/bin/env python
"""
The script can load a gray scale fits image with fits2aperture(.fits)
or create a circular aperture using aperture(size).
It will then produce a PSF .fits file, optically distorted with defocus, astigmatism, coma, trefoil and spherical
(Zernike modes 4 to 11, from Noll (1976))


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
import aperture
import psf as point_spread_function


def main():
    wavelength = float(raw_input('Wavelength? (0.76-2.00 microns) ') or .76)

    if glob.glob('aperture.fits'):
        ap = aperture.fits2aperture('aperture.fits')
    else:
        ap = aperture.aperture(101, 'HST')

    if glob.glob('distortions.par'):
        dist = point_spread_function.read_distortions('distortions.par')  # par file containing list of coefficients
    else:
        dist = [0., 0., 0.0026, 0.0089, 0.0222, -0.0018, 0.0113, 0., 0., 0.]

    psf = point_spread_function.compute_psf(ap, wavelength, 5, dist)

    # pixel_scale = 0.110  # of the instrument, in arcseconds per pixel

    # newPSF = resizePSF(psf, wavelength, size=np.shape(psf)[0], pixel_scale)  # re-bin the psf to detector size.

    point_spread_function.create_fits(psf, wavelength=wavelength, dist=dist)
    return

# ==================================================

if __name__ == '__main__':
    main()