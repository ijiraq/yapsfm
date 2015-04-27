#! /usr/bin/env python
"""
Main function to create PSF model.
Asks for user input on *wavelength* and pixel *scale* (or resolution) in ''/pixel.
It then computes an on-axis PSF model using the file distortions.par as source of optical distortions.
"""

import numpy as np

from modules import Pupil, PSF


def main():
    wavelength = float(raw_input('wavelength? (0.76 to 2.00 microns) ') or .76)
    scale = float(raw_input("final pixel scale in ''/pixel? (0.01) ") or 0.01)
    pupil = Pupil(wavelength)

    psf = PSF(pupil)
    psf.resize_psf(wavelength=wavelength, size=np.shape(psf.a)[0], scale=scale)
    psf.save()  # does not fill the header of the .fits image (yet)

if __name__ == "__main__":
    main()