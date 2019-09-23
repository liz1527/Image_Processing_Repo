from __future__ import print_function, division
import astropy.io.fits as pyfits
import numpy as np
import sys
import subprocess
from scipy import signal
import weightedstats as ws
from skimage.filters.rank import median as skmed

"""
__author__ = 'Will Hartley'

Code to clean up the VISTA VIDEO single-chip coadds' backgrounds.

Steps:
- run sextractor, with some default set of params
- expand segmentation map by 4 px (choice of 6px for 0.1342" / px is optimal for UDS, VISTA pixels are 0.2636"). Actually using a convolution with a square array
- use expanded segmap as mask in a median filtering of the background
- Filter is a simple square (129x129 px, i.e., 34x34 arcsec) 
- to make the algorithm run easier, we should embed the image in a larger array.
"""

# Image Object
# properties: file_name, weightmap, size, weight_name, mask, cleaned_im, cleaned_name, conf (configuration constants - see __main__)
# methods: run sextractor, embed image, make mask (inc. read segmap, expand segmap), clean image, save cleaned

class video_image:

    def __init__(self, fname, conf):

        self.fname = fname
        self.conf = conf

        # read the image, set-up weightmap filename
        self.read_image()


    def pad_array(self, im_array):
        # pad the array on all sides with the number of px in conf
        # and return the padded array
        pad_array = np.zeros((im_array.shape[0]+self.conf.pad_px*2,im_array.shape[1]+self.conf.pad_px*2))
        pad_array[self.conf.pad_px-1:pad_array.shape[0]-self.conf.pad_px-1, self.conf.pad_px-1:pad_array.shape[1]-self.conf.pad_px-1] = im_array
        return pad_array
    
        
    def read_image(self):
        # read header keywords and image data, pad image data, construct weightmap filename 
        with pyfits.open(self.fname) as f:
            self.im_size = f[0].header['NAXIS2'], f[0].header['NAXIS1'] # axes switched in python w.r.t. fits expectation.
            self.im_data = self.pad_array(f[0].data)
            self.im_head = f[0].header
        self.weight_name = '.'.join(fname.split('.')[:-1])+'.weight.fits' # assumes weight extension, could allow this to change via conf.
        self.clean_name = '.'.join(fname.split('.')[:-1])+'.cleaned.fits'
        
            
    def run_sex(self):
        # run sextractor
        subprocess.call('sex -c VIDEO.sex {0} -WEIGHT_IMAGE {1}'.format(self.fname, self.weight_name), shell=True)


    def expand_seg(self, seg):
        seg = signal.convolve2d(seg, np.ones((self.conf.exp_seg,self.conf.exp_seg)), mode='same')
        seg = np.floor(seg/self.conf.exp_seg**2)
        return seg

    
    def save_mask(self):
        # save the mask
        hdu = pyfits.PrimaryHDU(self.mask)
        hdu.writeto('mask.fits', clobber=True)
        

    def read_mask(self):
        # read in segmap
        segmap = pyfits.open('seg.fits')[0].data
        
        # change values such that source pixels are 0 and background are 1
        segmap[segmap==1] = 2
        segmap[segmap==0] = 1
        segmap[segmap>1] = 0

        # expand sources in segmap by pre-defined amount
        segmap = self.expand_seg(segmap)

        # embed in larger array
        self.mask = self.pad_array(segmap)

        # mask border region
        self.mask[0:self.conf.pad_px+self.conf.border,:] = 0
        self.mask[self.mask.shape[0]-(1+self.conf.pad_px+self.conf.border):,:] = 0
        self.mask[:, 0:self.conf.pad_px+self.conf.border] = 0
        self.mask[:, self.mask.shape[1]-(1+self.conf.pad_px+self.conf.border):] = 0

        # For testing purposes, we might want to save the mask
        if self.conf.bug_check == True:
            self.save_mask()


    def save_bugcheck_im(self, im, name):
        # save the bugcheck image
        hdu = pyfits.PrimaryHDU(im)
        hdu.writeto(name, clobber=True)
        

    def clean_im_BAD(self):
        # THIS DOESN'T WORK AT ALL - BIT DEPTH IS INSUFFICIENT
        
        # using masked median filter from scikit
        # filter the image, rejecting source pixels and masked regions
        tmp_im = self.im_data
        tmp_im[self.mask==0] = 0.
        tmp_im[tmp_im > self.conf.clip] = self.conf.clip
        tmp_im[tmp_im < -1.* self.conf.clip] = -1. * self.conf.clip

        # if testing, save the temporary image
        if self.conf.bug_check == True:
            self.save_bugcheck_im(tmp_im, 'tmp.fits')

        # invert mask for use with the scikit function - don't need
        #self.mask = -1 * self.mask + 1 
        
        # scale to range (-1, 1)
        im_max = np.max(np.abs(tmp_im))
        tmp_im /= im_max
        self.cleaned_im = skmed(tmp_im, selem=np.ones((self.conf.exp_seg*2+1,self.conf.exp_seg*2+1)), mask=self.mask)
        self.cleaned_im = self.cleaned_im.astype(np.float) * im_max # this is a lost cause

        # invert back - don't need
        #self.mask = -1 * (self.mask - 1.) 

        # if testing, save the intermediate background image
        if self.conf.bug_check == True:
            self.save_bugcheck_im(self.cleaned_im, 'bkgnd.fits')
        
        # subtract from original data
        self.cleaned_im = self.im_data - self.cleaned_im

        # de-pad the array (because that is what I assume when saving)
        self.cleaned_im = self.depad_array(self.cleaned_im)
        

    def fix_px(self, px_val, data, weights):
        # weighted median. Masked pixels have weight zero.
        try:
            return px_val - np.average(data, weights=weights)
        except:
            return px_val
        
        
    def fix_px_med(self, px_val, data, weights):
        # weighted median. Masked pixels have weight zero.
        try:
            return px_val - ws.numpy_weighted_median(data, weights=weights)
        except:
            return px_val
        
        
    def clean_im(self):
        # correct the background via a square median filter, ignoring masked pixels
        self.cleaned_im = np.zeros((self.im_size[0], self.im_size[1]))
        for i in range(self.im_size[0]):
            for j in range(self.im_size[1]):
                # the slices are of padded arrays, and extend the number of pixels given
                # by the half-size of the filter, in self.conf.filt_size, either side of
                # the target pixel. Slice dimension is always odd. 
                self.cleaned_im[i,j] = self.fix_px(self.im_data[i+self.conf.pad_px-1,j+self.conf.pad_px-1], self.im_data[i+self.conf.pad_px-1-self.conf.filt_size:i+self.conf.pad_px-1+self.conf.filt_size+1,j+self.conf.pad_px-1-self.conf.filt_size:j+self.conf.pad_px-1+self.conf.filt_size+1], self.mask[i+self.conf.pad_px-1-self.conf.filt_size:i+self.conf.pad_px-1+self.conf.filt_size+1,j+self.conf.pad_px-1-self.conf.filt_size:j+self.conf.pad_px-1+self.conf.filt_size+1])

    
    def depad_array(self, im_array):
        # remove the padded region
        depad_array = np.zeros((im_array.shape[0]-self.conf.pad_px*2,im_array.shape[1]-self.conf.pad_px*2))
        depad_array = im_array[self.conf.pad_px-1:im_array.shape[0]-self.conf.pad_px, self.conf.pad_px-1:im_array.shape[1]-self.conf.pad_px]
        return depad_array

    
    def save_clean_im(self):
        # save the cleaned image - unpad the array first. Copy original header
        hdu = pyfits.PrimaryHDU(self.cleaned_im) # don't need to depad the array, done during cleaning
        hdu.header = self.im_head
        hdu.writeto(self.clean_name, clobber=True)
    

# main prog.

if __name__ == '__main__':

    # get filename arg.
    if len(sys.argv) != 2:
        print('Error: No filename supplied.')
        print('Call as: python clean_im.py <filename>')

    else:
        fname = sys.argv[1]

    # here are some constants that can be tuned (Using a lambda function as an empty object)
    conf = lambda: None
    conf.exp_seg = 8 # number of pixels by which to expand the seg map in mask making (conv. w. sqyare array, so need twice the required n-px expansion).
    conf.pad_px = 65 # number of pixel to pad the image with - to simplify coding.
    conf.filt_size = 64 # half-dimension of the square filter
    conf.border = 175 # depth in pixels that forms mask region around the border
    conf.clip = 500. # 
    conf.bug_check = False # do we want to save bug checking images?
    print('Using config: exp_seg = {0}, pad_px = {1}, (half-)filt_size = {2}, border = {3}, clip_value = {4}, bug_check = {5}'.format(conf.exp_seg, conf.pad_px, conf.filt_size, conf.border, conf.clip, conf.bug_check))

    # create the VISTA image object (will read the data and set up weightmap filename
    im = video_image(fname, conf)

    # run sextractor
    im.run_sex()

    # read in the segmap and prepare it as a mask
    im.read_mask()

    # clean the image
    im.clean_im()

    # save cleaned image
    im.save_clean_im()
    
