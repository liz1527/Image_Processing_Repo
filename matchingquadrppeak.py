#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 12:24:26 2018

@author: ppxee
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import CircularAperture, aperture_photometry
import vari_funcs
plt.close('all')

def radial_profile(data, center):
    y, x = np.indices((data.shape)) #create coordinate grid
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2) #get radius values for grid
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel()) # counts number of times value
                                                # of radius occurs in the psf
                                                # weighted by the data
    nr = np.bincount(r.ravel()) # counts number of radii values in psf
    radialprofile = tbin / nr # as weighted is r*data then get profile by 
                              # dividing by unweighted counts of r values.
    return radialprofile 

def FWHM2sigma(FWHM, const):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    FWHM /= const
    return FWHM/np.sqrt(8*np.log(2))

def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))

def convolve_psf(sem, psf, newpsf, sigmakernel):
    ## Open image ###
    im05B = psf[sem]
    
    ## Convolve Image ###
    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf[sem] = convolve(im05B, kernel, normalize_kernel=True) 

def convolve_one_psf(psf, sigmakernel):
#    print(np.sum(psf))
    ## Convolve Image ###
#    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf = convolve(psf, kernel, normalize_kernel=True) 
    newpsf /= np.nansum(newpsf)
#    print(np.sum(newpsf))
    return newpsf

def quadrants(initdata,sem):
    
    ira = initdata['X_IMAGE_'+sem]
    idec = initdata['Y_IMAGE_'+sem]

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec < middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec >= middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad4data = initdata[mask1*mask2]
    
    return quad1data, quad2data, quad3data, quad4data

#sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
quads = np.arange(1,5,1)#np.array([1,3])
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
centre = [29,29]

for n, sem in enumerate(semesters):
    # limit data
    mag = oldsdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 19 #removes very faint stars
    mask = mask1 * mask2
    if sem == '05B':
        ids = oldsdata['NUMBER_05B'][mask]
    else:
        ids = np.intersect1d(ids, oldsdata['NUMBER_05B'][mask])

mask = np.isin(oldsdata['NUMBER_05B'], ids)
oldsdata = oldsdata[mask] 

for k, sem in enumerate(semesters):
    ### get quaddata
    quaddata = quadrants(oldsdata,sem)
    
    avgFWHM = np.zeros(len(quads))
    peakrp = np.zeros(len(quads))
    oldavgFWHM = np.zeros(len(quads))
    psf = {}
    oldpsf = {}
    oldrp = {}
    for n, quad in enumerate(quads):
        # for new
        colnames = colname+sem
        # for old
        oldmag = quaddata[n]['MAG_APER_'+sem][:,4]
        oldavgFWHM[n] = np.median(quaddata[n][colnames]) * 3600
        oldpsf[quad] = fits.open('PSFs/quadPSFs/1519_'+sem+'_'+str(quad)+'_K_PSF.fits')[0].data
    #    if sem == '10B':
    #        oldpsf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
    #    else:
    #        oldpsf[sem] = fits.open('PSFs/extra_'+sem+'_K_PSF.fits')[0].data
        oldrp[quad] = radial_profile(oldpsf[quad], centre)
        peakrp[n] = np.nanmax(oldrp[quad])

    ### Find minimum rp peak value as this is what all the others willl become ###
    aimind = np.argmin(peakrp)
    aimquad = quads[aimind]
    aimpsf = oldpsf[aimquad]
    aimrp = oldrp[aimquad]
    sqrtaimrp = np.sqrt(aimrp)
    
    ### Convert FWHM into a sigma ###
    sigmaold = np.array([fluxrad2sigma(fwhm) for fwhm in oldavgFWHM])
    sigmabroad = sigmaold[aimind]
    
    phot = {}
    flux = np.zeros(len(quads))
    oldphot = {}
    oldflux = np.zeros(len(quads))
    
    ### testing the extra factor method ###
    tests =np.linspace(0,1.1,1000)
    r = np.arange(0,42,1) * const * 3600 # define radius values   
    
    newpsf = {}
    newphot = {}
    newflux = np.zeros(len(quads))
    extras = np.zeros(len(quads))
    sigmafinal = np.zeros(len(quads))
    aperflux = np.empty(len(quads))
    oldaperflux = np.empty(len(quads))
    pixelr = (0.5/3600) / const
    aperture = CircularAperture(centre, pixelr)
    t = np.linspace(1, 4, num=4)
    
    plt.figure(1)
    plt.subplot(4,2,k+1)
    plt.plot(r, sqrtaimrp,label=aimquad)   
    plt.figure(2)
    plt.subplot(4,2,k+1)
    plt.hlines(0,0,1.5,label=str(aimquad))
    for n, quad in enumerate(quads):
        if quad == aimquad:
#            phot = aperture_photometry(aimpsf, aperture)
#            aperflux[n] = phot['aperture_sum'][0]
#            oldaperflux[n] = phot['aperture_sum'][0]
    #        plt.figure()
    #        plt.imshow(aimpsf)
    #        vari_funcs.no_ticks()
            continue
        
        sigma = sigmaold[n]
        singlephot = {}
        singleflux = np.zeros(len(tests))
        sumdiffold = 10
        for m, extra in enumerate(tests):
    #        print(extra)
            basesig = sigmabroad**2 - sigma**2
            
            # make sure extra is large enough to make the sqrt +ve
            if (basesig < 0 and extra < np.abs(basesig)) or (extra <0 and basesig < np.abs(extra)):
#                print('yep')
                continue
            else:
                sigmakernel = np.sqrt(basesig + extra)
                
            if sigmakernel <= 0:
#                print('yep')
                continue
#            print('no')
            new = convolve_one_psf(oldpsf[quad], sigmakernel)    
            
            ### Get radial profile ###
            radialprofile = radial_profile(new, centre)
            sqrtrp = np.sqrt(radialprofile)
            
            diff = aimrp[:12] - radialprofile[:12]
            sumdiff = np.nansum(diff)
    #        print(sumdiff)
            
            ### Plot the profile and difference so can see pattern ###
#            plt.figure(k+1)
#            plt.plot(r,sqrtrp, '--', label=str(quad)+' '+str(extra))
#            plt.ylabel('sqrt(Flux)')
#            plt.xlabel('Radius (arcsec)')
##            plt.xlim(xmax=1.5)
#            plt.legend()
#            plt.figure(k+9)
#            plt.plot(r[:12], diff, label=str(quad))
#            plt.xlim(xmin=0, xmax=1.5)
#            plt.legend()
#            plt.ylabel('Difference from aim rp')
#            plt.xlabel('Radius (arcsec)')
#            plt.tight_layout()
            if sumdiff > 0:
                if sumdiffold < sumdiff:
                    extras[n] = tests[m-1]
                    # make sure correct things plotted
                    sigmafinal[n] = np.sqrt(sigmabroad**2 - sigma**2+ extras[n]) 
                    radialprofile = radial_profile(new, centre)
                    sqrtrp = np.sqrt(radialprofile)
                    diff = aimrp[:12] - radialprofile[:12]
                else:
#                    print('extra ='+str(extra))
                    extras[n] = extra    
                    sigmafinal[n] = sigmakernel
#                
                plt.figure(1)
                plt.subplot(4,2,k+1)
    #            plt.plot(r, sqrtaimrp,label='10B')
                plt.plot(r,sqrtrp, '--', label=str(quad)+' '+str(extras[n]))
                plt.ylabel('sqrt(Flux)')
                plt.xlabel('Radius (arcsec)')
    #            plt.xlim(xmax=1.5)
                plt.legend()
                plt.figure(2)
                plt.subplot(4,2,k+1)
                plt.plot(r[:12], diff, label=str(quad))
                plt.xlim(xmin=0, xmax=1.5)
                plt.legend()
                plt.ylabel('Difference from aim rp')
                plt.xlabel('Radius (arcsec)')
#                plt.tight_layout()
                break
            else:
                sumdiffold = sumdiff    
                
    print(sem + str(extras))
    np.save('quadpeakextras'+sem, extras)