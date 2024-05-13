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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sunpy
import sunpy.map
from sunpy.map import Map
from sunpy.net import Fido, attrs as a
from sunpy.time import parse_time
import cupy as cp

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

# Date
time = '2024-04-14T22:45:09'
# Wavelengths
lambda_a = 171
lambda_b = 335

dates = ['2024-04-14T22:45:09','2024-04-15T22:45:45','2024-04-16T10:28:57','2024-04-17T22:45:33','2024-04-18T22:45:09','2024-04-19T22:45:09','2024-04-20T22:45:33']
l1pt5_171 = []
l1pt5_335 = []

def prep(smap):
        return register(update_pointing(smap))

def lvl_1pt5(day, l1, l2):

    ###################################
    ## Obtaining AIA Data ##
    t_start = parse_time(day)
    search_results = Fido.search(
        a.Time(t_start, t_start+11*u.s),
        a.Instrument.aia,
        a.Wavelength(l1*u.angstrom) | a.Wavelength(l2*u.angstrom),
    )
    search_results

    files = Fido.fetch(search_results, max_conn=1)

    m_a, m_b = sunpy.map.Map(sorted(files))

    ###################################
    ## PSF Deconvolution ##
    psf_a = psf(m_a.wavelength)
    psf_b = psf(m_b.wavelength)

    m_a_deconvolved = deconvolve(m_a, psf=psf_a)
    m_b_deconvolved = deconvolve(m_b, psf=psf_b)

    ## Respiking LVL 1 Images ##
    m_a_respiked = respike(m_a)
    m_b_respiked = respike(m_b)

    ###################################
    ## Transforming LVL 1 Images to LVL
    ## 1.5
    ## prep() Updates the pointing and 
    ## registers the images in one step

    m_a_L15 = prep(m_a)
    m_b_L15 = prep(m_b)

    ###################################
    ## Degrading Correction ##
    ## Extra Corrections Begin
    t_begin = parse_time('2010-03-25T00:00:00')
    now = astropy.time.Time.now()
    time_window = t_begin + np.arange(0, (now - t_begin).to(u.day).value, 7) * u.day

    correction_table = get_correction_table()

    d_a = degradation(m_a.wavelength, time_window, correction_table=correction_table)
    d_b = degradation(m_b.wavelength, time_window, correction_table=correction_table)

    d_a_map = degradation(m_a.wavelength, m_a.date, correction_table=correction_table)
    d_b_map = degradation(m_b.wavelength, m_b.date, correction_table=correction_table)

    m_a_corrected = correct_degradation(m_a, correction_table=correction_table)
    m_b_corrected = correct_degradation(m_b, correction_table=correction_table)

    l1pt5_171.append(m_a_corrected)
    print("\n171 Image Appended\n")

    l1pt5_335.append(m_b_corrected)
    print("\n335 Image Appended\n")

    ##Checking Functionality: print("\n\n Here is this iteration's 171: ")
    ##Checking Functionality: print(m_a_corrected)
    ##Checking Functionality: print("\n")
    ##Checking Functionality: 
    ##Checking Functionality: print("Here is this iteration's 335: ")
    ##Checking Functionality: print(l1pt5_335)
    ##Checking Functionality: print("\n\n")

    # fig = plt.figure()
    # ax = fig.add_subplot(projection=m_a_corrected)
    # m_a_corrected.plot(axes=ax)
    # plt.colorbar()
    # plt.show()

    return

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
    lvl_1pt5(i, lambda_a, lambda_b)

print("Finished Now")


###################################
# Creates a sequencing variable, 
# and sets which images are going to be sequenced
map_seq = sunpy.map.Map(l1pt5_171, sequence=True)  
ani = map_seq.plot(interval=1000)   
plt.show()

###################################
# Creates a sequencing variable, 
# and sets which images are going to be sequenced
map_seq = sunpy.map.Map(l1pt5_335, sequence=True)  
ani = map_seq.plot(interval=1000)   
plt.show()




##Check Extra: print("\nDone!\n ")
##Check Extra: 
##Check Extra: print("Here is 171: ")
##Check Extra: print(l1pt5_171)
##Check Extra: print("\n")
##Check Extra: 
##Check Extra: print("Here is 335: ")
##Check Extra: print(l1pt5_335)
##Check Extra: print("\n")


