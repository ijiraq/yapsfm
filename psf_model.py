#! /usr/bin/env python
"""
Main function to create PSF model.
Asks for user input on *wavelength* and pixel *scale* (or resolution) in ''/pixel.
It then computes an on-axis PSF model using the file distortions.par as source of optical distortions.
"""

import numpy as np
from argparse import ArgumentParser
from modules import Pupil, PSF


def main():
    usage = "Creates a PSF model for wavelength at input scale. If no option is specified, " \
            "user will be asked for input."
    parser = ArgumentParser(usage=usage)
    parser.add_argument('-w', type=float, dest='wavelength', default=None, action='store',
                        help='wavelength at which to compute the PSF. Between 0.76 and 2.00 microns.')
    parser.add_argument('-s', type=float, dest='scale', default=None, action='store',
                        help='final pixel scale in ''/pixel.')
    args = parser.parse_args()

    wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns) ') or .76) if args.wavelength is None else \
        args.wavelength
    scale = float(raw_input("final pixel scale in ''/pixel? (0.01) ") or 0.01) if args.scale is None else args.scale

    pupil = Pupil(wavelength)

    psf = PSF(pupil, scale)
    psf.resize_psf(wavelength=wavelength, size=np.shape(psf.a)[0], scale=scale)
    psf.save()  # does not fill the header of the .fits image (yet)

if __name__ == "__main__":
    main()