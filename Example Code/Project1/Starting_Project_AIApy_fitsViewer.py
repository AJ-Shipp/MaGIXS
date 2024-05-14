### Website version: https://sunpy.org/posts/2022/2022-01-06-aiapy-demo/#Wavelength-Response-Functions ###

## aiapy: A SunPy affiliated package for analyzing data from the Atmospheric Imaging Assembly ##

###################################
## Importing all needed packages &
## Outputting the current versions 
## of astropy, sunpy, and aiapy
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.time
from astropy.visualization import time_support, ImageNormalize, LogStretch
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy
import sunpy.map
import sunpy.io._fits
from sunpy.map import Map
from sunpy.net import Fido, attrs as a
from sunpy.time import parse_time
import cupy as cp

import fits_Files

import aiapy
from aiapy.psf import psf, deconvolve
from aiapy.calibrate import (register,update_pointing,correct_degradation, estimate_error,
                             degradation,normalize_exposure, respike, fetch_spikes)
from aiapy.calibrate.util import get_correction_table
from aiapy.response import Channel

import matplotlib as mpl

print(f'astropy v{astropy.__version__}')
print(f'sunpy v{sunpy.__version__}')
print(f'aiapy v{aiapy.__version__}')

# Counter for the for-loop
num = 1

# Date
time = '2024-04-14T22:45:09'
# Wavelengths
lambda_a = 171
lambda_b = 335

dates = ['2024-04-14T22:45:09','2024-04-15T22:45:45','2024-04-16T10:28:57','2024-04-17T22:45:33','2024-04-18T22:45:09','2024-04-19T22:45:09','2024-04-20T22:45:33']
l1pt5_171 = []
l1pt5_335 = []

#: l1pt5_171.append(data)

print(l1pt5_171)
#: print(data)

###################################
## Prints the 171 and 335 lists
## This should be empty... should...
print(l1pt5_171)
print(l1pt5_335)

###################################
## Iterates through the list of 
## necessary dates and runs the 
## function to prepare the lvl 1.5
## images 
for i in dates:
    img = sunpy.map.Map('Example Code/Project1/fits_Files/aia171_'+str(num)+'.fits')
    l1pt5_171.append(img)
    img = sunpy.map.Map('Example Code/Project1/fits_Files/aia335_'+str(num)+'.fits')
    l1pt5_335.append(img)
    #("aia335_"+str(num)+".fits")
    print("Iteration #",num,"finished")
    num = num + 1

print("Finished Now")

###################################
## Creates a sequencing variable, 
## and sets which images are going 
## to be sequenced
## Also saves the sequence to a file
map_seq = sunpy.map.Map(l1pt5_171, sequence=True)  
ani = map_seq.plot(interval=1000)
plt.show()

###################################
## Creates a sequencing variable,
## and sets which images are going
## to be sequenced
## Also saves the sequence to a file
map_seq = sunpy.map.Map(l1pt5_335, sequence=True)
ani = map_seq.plot(interval=1000)
plt.show()