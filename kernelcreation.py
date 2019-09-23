#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 11:36:12 2018

Code to create a kernel that will match the PSFs of two stacks

@author: ppxee
"""

### Import Modules ###
import numpy as np
from photutils import create_matching_kernel, TopHatWindow, CosineBellWindow
import matplotlib.pyplot as plt
from astropy.io import fits
plt.close('all')

psf11 = fits.open('11B_K_PSF.fits')[0].data
psf12 = fits.open('12B_K_PSF.fits')[0].data

#plt.figure()
#plt.imshow(psf12-psf11)
#plt.figure()
#plt.imshow(psf11-psf12)


kernel = create_matching_kernel(psf11, psf12, window=CosineBellWindow(0.35)) #need to check which is larger

plt.figure()
plt.imshow(kernel)


