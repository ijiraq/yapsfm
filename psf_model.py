#! /usr/bin/env python
"""
Main function to create PSF model.
Asks for user input on *wavelength* and pixel *scale* (or resolution) in ''/pixel.
It then computes an on-axis PSF model using the file distortions.par as source of optical distortions and aperture.fits
as aperture source file. Aperture.fits must be located in the same folder as this script.
"""

import numpy as np
from argparse import ArgumentParser
import glob
import logging
from astropy.io import fits
from modules import Pupil, PSF, PolyPSF


def main():
    usage = "Creates a PSF model for wavelength at input scale.\ncommand line options:\n" \
            "for monochrome: psf_model.py -s -w -p -c\n" \
            "for polychrome: psf_model.py -s -b -t -p -c"
    parser = ArgumentParser(usage=usage)
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-s', '--scale', type=float, dest='scale', default=None, action='store',
                        help="final pixel scale in ''/pixel. Required.")
    group.add_argument('-w', '--wavelength', type=float, dest='wavelength', default=None, action='store',
                       help='wavelength at which to compute the PSF. Between 0.76 and 2.00 microns. '
                            'Required if monochromatic PSF.')
    group.add_argument('-b', '--band', type=str, dest='band', default=None, action='store',
                       help="filter band to use for polychromatic PSF. Z, Y, J, H, F, Wide. "
                            "Required if polychromatic PSF.")
    parser.add_argument('-t', '--type', type=str, dest='spectral_type', default=None, action='store',
                        help="spectral type of target star. B, A, F, G, K, M handled. "
                             "Required if polychromatic PSF.")
    parser.add_argument('-v', '--verbose', type=str, dest='verbose', default='info', action='store',
                        help="verbose level, from most to least exhaustive: debug, info, warning, error, critical. "
                             "Default: info")
    parser.add_argument('-a', '--aperture', type=str, dest='aperture', default='aperture.fits', action='store',
                        help="the name of the aperture .fits file to use during PSF creation. Default=%(default)s")
    parser.add_argument('-i', '--individual', dest='switch', action='store_true',
                        help="switch to save all the individual PSFs used to create the polychrome.")
    parser.add_argument('-p', '--position', nargs='+', dest='position', default=[0, 0], action='store',
                        help="position of the computed PSF on the detector. syntax: -p x y with 0 < x y < 4088. Center"
                             " of detector is at 2044 2044."
                             "Required.")
    parser.add_argument('-c', '--chip', type=int, dest='chip', default=None, action='store',
                        help="chip number (1 to 18). Required")
    parser.add_argument('-j', '--jitter', type=float, dest='jitter', default=0.01, action='store',
                        help="jitter (in arcseconds) affecting the PSF. Default=%(default)s")
    args = parser.parse_args()

    # defining verbose level and output options
    levels = {'debug': logging.DEBUG,
              'info': logging.INFO,
              'warning': logging.WARNING,
              'error': logging.ERROR,
              'critical': logging.CRITICAL}
    level = levels.get(args.verbose, logging.NOTSET)
    if args.verbose == 'info':
        logging.basicConfig(format='%(message)s', level=level)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=level)

    # required args
    if args.band and args.spectral_type is None:
        parser.error("--band (-b) requires --type (-t).")
    if args.position is None or args.chip is None:
        parser.error("--position (-p) and --chip (-c) are required.")

    scale = float(raw_input("final pixel scale in ''/pixel? (0.01) ") or 0.01) if args.scale is None else args.scale

    if glob.glob(args.aperture):  # if there's an aperture.fits file, open it, otherwise, create an HST-like
        with fits.open(args.aperture) as hdu:
            array = hdu[0].data
            size = np.shape(array)[0]  # aperture.fits is a square array
            aperture_name = args.aperture
    else:
        logging.info("'%s' not found, creating new HST-like aperture: size = 101 pixels" % args.aperture)
        size = 101
        aperture_name = None
    logging.debug("array size: %s" % size)

    # mono- or polychromatic
    if args.wavelength:  # if wavelength : monochromatic
        wavelength = args.wavelength
        pupil = Pupil(wavelength, size, aperture_name, args.position, args.chip)
        psf = PSF(pupil, scale, size, args.position, args.chip, args.jitter)
        psf.resize_psf()
        psf.save(name="psf_%s" % wavelength)
    elif args.band:  # if band : polychromatic
        poly = PolyPSF(band=args.band, spectral_type=args.spectral_type, scale=scale, size=size, aperture=aperture_name,
                       position=args.position, chip=args.chip)
        poly.get_sed()
        poly.wavelength_contributions()
        poly.create_polychrome(args.switch, args.jitter)
        poly.save()

if __name__ == "__main__":
    main()