### Website version: https://sunpy.org/posts/2022/2022-01-06-aiapy-demo/#Wavelength-Response-Functions ###

## aiapy: A SunPy affiliated package for analyzing data from the Atmospheric Imaging Assembly ##
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

## Obtaining AIA Data ##
t_start = parse_time('2024-04-14T22:45:09')
search_results = Fido.search(
    a.Time(t_start, t_start+11*u.s),
    a.Instrument.aia,
    a.Wavelength(171*u.angstrom) | a.Wavelength(335*u.angstrom),
)
search_results

files = Fido.fetch(search_results, max_conn=1)

m_171_a14, m_335_a14 = sunpy.map.Map(sorted(files))

## PSF Deconvolution ##
psf_171_a14 = psf(m_171_a14.wavelength)
psf_335_a14 = psf(m_335_a14.wavelength)

m_171_a14_deconvolved = deconvolve(m_171_a14, psf=psf_171_a14)
m_335_a14_deconvolved = deconvolve(m_335_a14, psf=psf_335_a14)

## Respiking LVL 1 Images ##
m_171_a14_respiked = respike(m_171_a14)
m_335_a14_respiked = respike(m_335_a14)

## Transforming LVL 1 Images to LVL 1.5
m_171_a14_up = update_pointing(m_171_a14)
m_335_a14_up = update_pointing(m_335_a14)

m_171_a14_L15 = register(m_171_a14_up)
m_335_a14_L15 = register(m_335_a14_up)

def prep(smap):
    return register(update_pointing(smap))

m_335_a14_L15 = prep(m_335_a14)

## Degrading Correction ##
t_begin = parse_time('2010-03-25T00:00:00')
now = astropy.time.Time.now()
time_window = t_begin + np.arange(0, (now - t_begin).to(u.day).value, 7) * u.day

correction_table = get_correction_table()

d_171_a14 = degradation(m_171_a14.wavelength, time_window, correction_table=correction_table)
d_335_a14 = degradation(m_335_a14.wavelength, time_window, correction_table=correction_table)

d_171_a14_map = degradation(m_171_a14.wavelength, m_171_a14.date, correction_table=correction_table)
d_335_a14_map = degradation(m_335_a14.wavelength, m_335_a14.date, correction_table=correction_table)

m_171_a14_corrected = correct_degradation(m_171_a14, correction_table=correction_table)
m_335_a14_corrected = correct_degradation(m_335_a14, correction_table=correction_table)

fig = plt.figure()
ax = fig.add_subplot(projection=m_171_a14_corrected)
m_171_a14_corrected.plot(axes=ax)
plt.colorbar()
plt.show()

map_seq = sunpy.map.Map([m_171_a14_corrected, m_335_a14_corrected], sequence=True)  
ani = map_seq.plot()   
plt.show()
