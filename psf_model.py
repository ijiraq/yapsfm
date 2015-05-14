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
    usage = "Creates a PSF model for wavelength at input scale. If no option is specified, " \
            "user will be asked for input."
    parser = ArgumentParser(usage=usage)
    parser.add_argument('-s', type=float, dest='scale', default=None, action='store',
                        help="final pixel scale in ''/pixel. Required argument.")
    parser.add_argument('-p', type=bool, dest='polychromatic', default=0, action='store',
                        help='Boolean switch to output a polychromatic PSF. Default=0')
    parser.add_argument('-b', type=str, dest='band', default=None, action='store',
                        help="filter band to use for polychromatic PSF. 'Z', 'Y', 'J', 'H', 'F', 'Wide'.")
    parser.add_argument('-t', type=str, dest='spectral_type', default=None, action='store',
                        help="spectral type of target star. b, a, f, g, k, m handled. "
                             "Required if polychromatic PSF is computed")
    parser.add_argument('-w', type=float, dest='wavelength', default=None, action='store',
                        help='wavelength at which to compute the PSF. Between 0.76 and 2.00 microns. '
                             'Required if monochromatic PSF is computed.')
    args = parser.parse_args()

    scale = float(raw_input("final pixel scale in ''/pixel? (0.01) ") or 0.01) if args.scale is None else args.scale

    if args.polychromatic == 0:
        wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns) ') or .76) if args.wavelength is None else \
            args.wavelength

        pupil = Pupil(wavelength)
        psf = PSF(pupil, scale)
        psf.resize_psf(wavelength=wavelength, size=np.shape(psf.a)[0], scale=scale)
        psf.save()
    else:
        poly = PolyPSF(band=args.band, spectral_type=args.spectral_type, scale=args.scale)
        poly.get_sed()
        poly.wavelength_contributions()
        poly.create_polychrome()
        poly.save()

if __name__ == "__main__":
    main()