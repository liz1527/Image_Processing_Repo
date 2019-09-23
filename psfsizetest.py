#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 12:01:23 2018

Code to test appropriate VIGNET sizes

@author: ppxee
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from photutils import CircularAperture
from photutils import aperture_photometry
import vari_funcs
plt.close('all')

#sems = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
sem='10B' # test on worst and it should then fint the rest
psf = fits.getdata('PSFs/limited_'+sem+'_K_PSF.fits')


mask = np.zeros(np.shape(psf))
mask[33:91,33:91] = 1
mask = mask.astype(bool)
psfmasked = psf.copy()
psfmasked[~mask] = 1

#plot before and after
plt.imshow(np.log(psf))
plt.figure()
plt.imshow(np.log(psfmasked))