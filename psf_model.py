#! /usr/bin/env python
"""
Main function to create PSF model.
Asks for user input on *wavelength* and pixel *scale* (or resolution) in ''/pixel.
It then computes an on-axis PSF model using the file distortions.par as source of optical distortions.
"""

import numpy as np
from argparse import ArgumentParser
from modules import Pupil, PSF, PolyPSF


def main():
    usage = "Creates a PSF model for wavelength at input scale.\ncommand line options:\n" \
            "for monochrome: psf_model.py -s -w\n" \
            "for polychrome: psf_model.py -s -b -t"
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
    args = parser.parse_args()

    if args.band and args.spectral_type is None:
        parser.error("--band (-b) requires --type (-t).")

    scale = float(raw_input("final pixel scale in ''/pixel? (0.01) ") or 0.01) if args.scale is None else args.scale

    if args.wavelength:
        wavelength = args.wavelength
        pupil = Pupil(wavelength)
        psf = PSF(pupil, scale)
        psf.resize_psf(wavelength=wavelength, size=np.shape(psf.a)[0], scale=scale)
        psf.save(name="psf_%s" % wavelength)
    elif args.band:
        poly = PolyPSF(band=args.band, spectral_type=args.spectral_type, scale=scale)
        poly.get_sed()
        poly.wavelength_contributions()
        poly.create_polychrome()
        poly.save()

if __name__ == "__main__":
    main()