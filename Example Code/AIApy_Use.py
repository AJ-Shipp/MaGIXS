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
# Increases the figure size in this notebook.
mpl.rcParams["savefig.dpi"] = 150
mpl.rcParams["figure.dpi"] = 150

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

#= m_171_a14.peek(vmin=0)
#= m_335_a14.peek(vmin=0)

## PSF Deconvolution ##
psf_171_a14 = psf(m_171_a14.wavelength)
psf_335_a14 = psf(m_335_a14.wavelength)

#Testing: plt.imshow(psf_171_a14, origin='lower', norm=ImageNormalize(vmax=1e-6, stretch=LogStretch()))
#Testing: plt.colorbar()

m_171_a14_deconvolved = deconvolve(m_171_a14, psf=psf_171_a14)
m_335_a14_deconvolved = deconvolve(m_335_a14, psf=psf_335_a14)

#Testing: blc = SkyCoord(750,-375,unit='arcsec',frame=m_171_a14.coordinate_frame)
#Testing: fov = {'width': 400*u.arcsec, 'height': 400*u.arcsec}
#Testing: m_171_a14_cutout = m_171_a14.submap(blc, **fov)
#Testing: m_171_a14_deconvolved_cutout = m_171_a14_deconvolved.submap(blc, **fov)

#Testing: fig = plt.figure(figsize=(7,3))
#Testing: ax = fig.add_subplot(121, projection=m_171_a14_cutout)
#Testing: m_171_a14_cutout.plot(axes=ax, title='Before Deconvolution')
#Testing: ax = fig.add_subplot(122, projection=m_171_a14_deconvolved_cutout)
#Testing: m_171_a14_deconvolved_cutout.plot(axes=ax, title='After Deconvolution')
#Testing: ax.coords[1].set_axislabel(' ')

#Testing: x = np.linspace(m_171_a14_deconvolved_cutout.dimensions.x.value*0.55, m_171_a14_deconvolved_cutout.dimensions.x.value*0.7)
#Testing: y = 0.59 * m_171_a14_deconvolved_cutout.dimensions.y.value * np.ones(x.shape)
#Testing: sl = np.s_[np.round(y[0]).astype(int), np.round(x[0]).astype(int):np.round(x[-1]).astype(int)]

#Testing: fig = plt.figure(figsize=(7,3))
#Testing: ax = fig.add_subplot(121, projection=m_171_a14_deconvolved_cutout)
#Testing: m_171_a14_deconvolved_cutout.plot(axes=ax)
#Testing: ax.plot(x, y, lw=1)
#Testing: ax = fig.add_subplot(122)
#Testing: Tx = sunpy.map.all_coordinates_from_map(m_171_a14_cutout)[sl].Tx
#Testing: ax.plot(Tx,m_171_a14_cutout.data[sl], label='Original')
#Testing: ax.plot(Tx,m_171_a14_deconvolved_cutout.data[sl], label='Deconvolved')
#Testing: ax.set_ylabel(f'Intensity [{m_171_a14_cutout.unit}]')
#Testing: ax.set_xlabel(r'Helioprojective Longitude [arcsec]')
#Testing: ax.legend(loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5,1.15))
#Testing: plt.tight_layout()

## Respiking LVL 1 Images ##
m_171_a14_respiked = respike(m_171_a14)
m_335_a14_respiked = respike(m_335_a14)

#Testing: fig = plt.figure(figsize=(7, 3))
#Testing: ax = fig.add_subplot(121, projection=m_171)
#Testing: m_171_a14.plot(axes=ax)
#Testing: ax.set_title(f"Despiked (Level {m_171_a14.processing_level:.0f})")
#Testing: ax = fig.add_subplot(122, projection=m_171_a14_respiked)
#Testing: m_171_a14_respiked.plot(axes=ax)
#Testing: ax.set_title(f"Respiked (Level {m_171_a14_respiked.processing_level})")
#Testing: ax.coords[1].set_axislabel(' ')

#Testing: pix_coords, vals = fetch_spikes(m_171,)

#Testing: vals_despiked = m_171_a14.data[pix_coords.y.value.round().astype(int), pix_coords.x.value.round().astype(int)]

#Testing: plt.hist(vals, bins='scott', log=True, histtype='step', label='Respiked');
#Testing: plt.hist(vals_despiked, bins='scott', log=True, histtype='step', label='Despiked');
#Testing: plt.legend()
#Testing: plt.xlabel(f'Intensity [{m_171_a14.unit.to_string()}]')
#Testing: plt.ylabel('Number of Pixels')

#Testing: spike_coords = m_171_a14.pixel_to_world(pix_coords.x, pix_coords.y)

#Testing: fig = plt.figure()
#Testing: ax = fig.add_subplot(111, projection=m_171_a14_respiked)
#Testing: m_171_a14_respiked.plot(axes=ax)
#Testing: ax.plot_coord(spike_coords, marker='.', ls=' ', markersize=1)

## Transforming LVL 1 Images to LVL 1.5
m_171_a14_up = update_pointing(m_171_a14)
m_335_a14_up = update_pointing(m_335_a14)

#Testing: m_171_a14.reference_pixel

#Testing: m_171_a14_up.reference_pixel

m_171_a14_L15 = register(m_171_a14_up)
m_335_a14_L15 = register(m_335_a14_up)

#Testing: print(m_171_a14_up.scale)
#Testing: print(m_171_a14_up.rotation_matrix)

#Testing: print(m_171_a14_L15.scale)
#Testing: print(m_171_a14_L15.rotation_matrix)

def prep(smap):
    return register(update_pointing(smap))

m_335_a14_L15 = prep(m_335_a14)

#Testing: print(m_335_a14_L15.scale)
#Testing: print(m_335_a14_L15.rotation_matrix)

## Degrading Correction ##
t_begin = parse_time('2010-03-25T00:00:00')
now = astropy.time.Time.now()
time_window = t_begin + np.arange(0, (now - t_begin).to(u.day).value, 7) * u.day

correction_table = get_correction_table()

d_171_a14 = degradation(m_171_a14.wavelength, time_window, correction_table=correction_table)
d_335_a14 = degradation(m_335_a14.wavelength, time_window, correction_table=correction_table)

d_171_a14_map = degradation(m_171_a14.wavelength, m_171_a14.date, correction_table=correction_table)
d_335_a14_map = degradation(m_335_a14.wavelength, m_335_a14.date, correction_table=correction_table)

#Testing: d_335_a14_v9 = degradation(m_335_a14.wavelength, time_window, calibration_version=9, correction_table=correction_table)
#Testing: d_335_a14_v8 = degradation(m_335_a14.wavelength, time_window, calibration_version=8, correction_table=correction_table)

#Testing: with time_support(format='jyear'):
#Testing:     plt.plot(time_window, d_335, label='v10')
#Testing:     plt.plot(time_window, d_335_a14_v9, label='v9')
#Testing:     plt.plot(time_window, d_335_a14_v8, label='v8')
#Testing:     plt.plot(m_335_a14.date, d_335_a14_map, linestyle='', marker='o', color='C0', label=m_335_a14.date.iso)
#Testing: plt.ylabel('Degradation 335 Å')
#Testing: plt.legend()

m_171_a14_corrected = correct_degradation(m_171_a14, correction_table=correction_table)
m_335_a14_corrected = correct_degradation(m_335_a14, correction_table=correction_table)

#37: fig = plt.figure(figsize=(7, 3))
#37: ax = fig.add_subplot(121, projection=m_335)
#37: m_335_a14.plot(axes=ax, vmin=0, vmax=2.5e3,title='Uncorrected')
#37: ax = fig.add_subplot(122, projection=m_335_a14_corrected)
#37: m_335_a14_corrected.plot(axes=ax, vmin=0, vmax=2.5e3, title='Corrected')
#37: ax.coords[1].set_axislabel(' ')

fig = plt.figure()
ax = fig.add_subplot(projection=m_171_a14_corrected)
m_171_a14_corrected.plot(axes=ax)
plt.colorbar()
plt.show()

map_seq = sunpy.map.Map([m_171_a14_corrected, m_335_a14_corrected], sequence=True)  
ani = map_seq.plot()   
plt.show()

## Computing Uncertainties ##
#38: errors_171 = estimate_error(m_171_a14_L15.quantity/u.pix, m_171_a14_L15.wavelength)

#39: m_171_a14_errors = sunpy.map.Map(errors_171_a14.value, m_171_a14_L15.meta)

#40: m_171_a14_errors.peek(norm=ImageNormalize(vmax=50))

#41: errors_171_a14_chianti = estimate_error(m_171_a14_L15.quantity/u.pix, m_171_a14_L15.wavelength, include_chianti=True)

#42: errors_171_a14_eve = estimate_error(m_171_a14_L15.quantity/u.pix, m_171_a14_L15.wavelength, include_eve=True)

#43: hist_params = {'bins': np.logspace(0,4,50), 'histtype': 'step', 'log': True}
#43: plt.hist(errors_171_a14.value.flatten(), **hist_params, label='Nominal');
#43: plt.hist(errors_171_a14_chianti.value.flatten(), **hist_params, label='CHIANTI');
#43: plt.hist(errors_171_a14_eve.value.flatten(), **hist_params, label='Photometric (EVE)');
#43: plt.xlabel('Uncertainty [ct/pix]')
#43: plt.ylabel('Number of Pixels')
#43: plt.xscale('log')
#43: plt.legend(frameon=False)

## Exposure Time Normalization ##
#44: m_171_a14_norm = normalize_exposure(m_171_a14_L15)

#45: print(m_171_a14_L15.unit)
#45: print(m_171_a14_norm.unit)

#46: print(m_171_a14_L15.exposure_time)
#46: print(m_171_a14_norm.exposure_time)

#47: m_171_a14_norm = m_171_a14_L15 / m_171_a14_L15.exposure_time

#48: print(m_171_a14_norm.unit)

#49: print(m_171_a14_norm.exposure_time)

## Wavelength Response Functions ##
#50: c = Channel(m_335_a14.wavelength)

#51: print(c.channel)
#51: print(c.telescope_number)

#52: r = c.wavelength_response()

#53: print(c.geometrical_collecting_area)
#53: print(c.primary_reflectance)
#53: print(c.focal_plane_filter_efficiency)
#53: print(c.contamination)
#53: print(c.quantum_efficiency)
#53: print(c.gain)

#54: r_time = c.wavelength_response(obstime=m_335_a14.date)
#54: r_time_eve = c.wavelength_response(obstime=m_335_a14.date, include_eve_correction=True)

#55: plt.plot(c.wavelength,r,label='Uncorrected')
#55: plt.plot(c.wavelength,r_time,label='Degradation')
#55: plt.plot(c.wavelength,r_time_eve,label='Degradation + EVE')
#55: plt.xlim([315,355])
#55: plt.ylim([0,0.03])
#55: plt.xlabel('$\lambda$ [Å]')
#55: plt.ylabel(f'$R(\lambda)$ [{r.unit.to_string(format="latex_inline")}]')
#55: plt.legend(loc=2, frameon=False)
