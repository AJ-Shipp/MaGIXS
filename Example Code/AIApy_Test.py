## aiapy: A SunPy affiliated package for analyzing data from the Atmospheric Imaging Assembly ##
#1: import astropy
#1: import astropy.units as u
#1: from astropy.coordinates import SkyCoord
#1: import astropy.time
#1: from astropy.visualization import time_support, ImageNormalize, LogStretch
#1: import numpy as np
#1: import matplotlib.pyplot as plt
#1: import sunpy
#1: import sunpy.map
#1: from sunpy.net import Fido, attrs as a
#1: from sunpy.time import parse_time
#1: 
#1: import aiapy
#1: from aiapy.psf import psf, deconvolve
#1: from aiapy.calibrate import (register,update_pointing,correct_degradation, estimate_error,
#1:                              degradation,normalize_exposure, respike, fetch_spikes)
#1: from aiapy.calibrate.util import get_correction_table
#1: from aiapy.response import Channel
#1: 
#1: import matplotlib as mpl
#1: # Increases the figure size in this notebook.
#1: mpl.rcParams["savefig.dpi"] = 150
#1: mpl.rcParams["figure.dpi"] = 150

#2: print(f'astropy v{astropy.__version__}')
#2: print(f'sunpy v{sunpy.__version__}')
#2: print(f'aiapy v{aiapy.__version__}')

## Obtaining AIA Data ##
#3: t_start = parse_time('2017-09-10T20:00:00')
#3: search_results = Fido.search(
#3:     a.Time(t_start, t_start+11*u.s),
#3:     a.Instrument.aia,
#3:     a.Wavelength(171*u.angstrom) | a.Wavelength(335*u.angstrom),
#3: )
#3: search_results

#4: files = Fido.fetch(search_results, max_conn=1)

#5: m_171, m_335 = sunpy.map.Map(sorted(files))

#6: m_171.peek(vmin=0)
#6: m_335.peek(vmin=0)

## PSF Deconvolution ##
#7: psf_171 = psf(m_171.wavelength)

#8: plt.imshow(psf_171, origin='lower', norm=ImageNormalize(vmax=1e-6, stretch=LogStretch()))
#8: plt.colorbar()

#9: m_171_deconvolved = deconvolve(m_171, psf=psf_171)

#10: blc = SkyCoord(750,-375,unit='arcsec',frame=m_171.coordinate_frame)
#10: fov = {'width': 400*u.arcsec, 'height': 400*u.arcsec}
#10: m_171_cutout = m_171.submap(blc, **fov)
#10: m_171_deconvolved_cutout = m_171_deconvolved.submap(blc, **fov)

#11: fig = plt.figure(figsize=(7,3))
#11: ax = fig.add_subplot(121, projection=m_171_cutout)
#11: m_171_cutout.plot(axes=ax, title='Before Deconvolution')
#11: ax = fig.add_subplot(122, projection=m_171_deconvolved_cutout)
#11: m_171_deconvolved_cutout.plot(axes=ax, title='After Deconvolution')
#11: ax.coords[1].set_axislabel(' ')

#12: x = np.linspace(m_171_deconvolved_cutout.dimensions.x.value*0.55, m_171_deconvolved_cutout.dimensions.x.value*0.7)
#12: y = 0.59 * m_171_deconvolved_cutout.dimensions.y.value * np.ones(x.shape)
#12: sl = np.s_[np.round(y[0]).astype(int), np.round(x[0]).astype(int):np.round(x[-1]).astype(int)]

#13: fig = plt.figure(figsize=(7,3))
#13: ax = fig.add_subplot(121, projection=m_171_deconvolved_cutout)
#13: m_171_deconvolved_cutout.plot(axes=ax)
#13: ax.plot(x, y, lw=1)
#13: ax = fig.add_subplot(122)
#13: Tx = sunpy.map.all_coordinates_from_map(m_171_cutout)[sl].Tx
#13: ax.plot(Tx,m_171_cutout.data[sl], label='Original')
#13: ax.plot(Tx,m_171_deconvolved_cutout.data[sl], label='Deconvolved')
#13: ax.set_ylabel(f'Intensity [{m_171_cutout.unit}]')
#13: ax.set_xlabel(r'Helioprojective Longitude [arcsec]')
#13: ax.legend(loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5,1.15))
#13: plt.tight_layout()

## Respiking LVL 1 Images ##
#14: m_171_respiked = respike(m_171)

#15: fig = plt.figure(figsize=(7, 3))
#15: ax = fig.add_subplot(121, projection=m_171)
#15: m_171.plot(axes=ax)
#15: ax.set_title(f"Despiked (Level {m_171.processing_level:.0f})")
#15: ax = fig.add_subplot(122, projection=m_171_respiked)
#15: m_171_respiked.plot(axes=ax)
#15: ax.set_title(f"Respiked (Level {m_171_respiked.processing_level})")
#15: ax.coords[1].set_axislabel(' ')

#16: pix_coords, vals = fetch_spikes(m_171,)

#17: vals_despiked = m_171.data[pix_coords.y.value.round().astype(int), pix_coords.x.value.round().astype(int)]

#18: plt.hist(vals, bins='scott', log=True, histtype='step', label='Respiked');
#18: plt.hist(vals_despiked, bins='scott', log=True, histtype='step', label='Despiked');
#18: plt.legend()
#18: plt.xlabel(f'Intensity [{m_171.unit.to_string()}]')
#18: plt.ylabel('Number of Pixels')

#19: spike_coords = m_171.pixel_to_world(pix_coords.x, pix_coords.y)

#20: fig = plt.figure()
#20: ax = fig.add_subplot(111, projection=m_171_respiked)
#20: m_171_respiked.plot(axes=ax)
#20: ax.plot_coord(spike_coords, marker='.', ls=' ', markersize=1)

## Transforming LVL 1 Images to LVL 1.5
#21: m_171_up = update_pointing(m_171)

#22: m_171.reference_pixel

#23: m_171_up.reference_pixel

#24: m_171_L15 = register(m_171_up)

#25: print(m_171_up.scale)
#25: print(m_171_up.rotation_matrix)

#26: print(m_171_L15.scale)
#26: print(m_171_L15.rotation_matrix)

#27: def prep(smap):
#27:     return register(update_pointing(smap))

#28: m_335_L15 = prep(m_335)

#29: print(m_335_L15.scale)
#29: print(m_335_L15.rotation_matrix)

## Degrading Correction ##
#30: t_begin = parse_time('2010-03-25T00:00:00')
#30: now = astropy.time.Time.now()
#30: time_window = t_begin + np.arange(0, (now - t_begin).to(u.day).value, 7) * u.day

#31: correction_table = get_correction_table()

#32: d_335 = degradation(m_335.wavelength, time_window, correction_table=correction_table)

#33: d_335_map = degradation(m_335.wavelength, m_335.date, correction_table=correction_table)

#34: d_335_v9 = degradation(m_335.wavelength, time_window, calibration_version=9, correction_table=correction_table)
#34: d_335_v8 = degradation(m_335.wavelength, time_window, calibration_version=8, correction_table=correction_table)

#35: with time_support(format='jyear'):
#35:     plt.plot(time_window, d_335, label='v10')
#35:     plt.plot(time_window, d_335_v9, label='v9')
#35:     plt.plot(time_window, d_335_v8, label='v8')
#35:     plt.plot(m_335.date, d_335_map, linestyle='', marker='o', color='C0', label=m_335.date.iso)
#35: plt.ylabel('Degradation 335 Å')
#35: plt.legend()

#36: m_335_corrected = correct_degradation(m_335, correction_table=correction_table)

#37: fig = plt.figure(figsize=(7, 3))
#37: ax = fig.add_subplot(121, projection=m_335)
#37: m_335.plot(axes=ax, vmin=0, vmax=2.5e3,title='Uncorrected')
#37: ax = fig.add_subplot(122, projection=m_335_corrected)
#37: m_335_corrected.plot(axes=ax, vmin=0, vmax=2.5e3, title='Corrected')
#37: ax.coords[1].set_axislabel(' ')

## Computing Uncertainties ##
#38: errors_171 = estimate_error(m_171_L15.quantity/u.pix, m_171_L15.wavelength)

#39: m_171_errors = sunpy.map.Map(errors_171.value, m_171_L15.meta)

#40: m_171_errors.peek(norm=ImageNormalize(vmax=50))

#41: errors_171_chianti = estimate_error(m_171_L15.quantity/u.pix, m_171_L15.wavelength, include_chianti=True)

#42: errors_171_eve = estimate_error(m_171_L15.quantity/u.pix, m_171_L15.wavelength, include_eve=True)

#43: hist_params = {'bins': np.logspace(0,4,50), 'histtype': 'step', 'log': True}
#43: plt.hist(errors_171.value.flatten(), **hist_params, label='Nominal');
#43: plt.hist(errors_171_chianti.value.flatten(), **hist_params, label='CHIANTI');
#43: plt.hist(errors_171_eve.value.flatten(), **hist_params, label='Photometric (EVE)');
#43: plt.xlabel('Uncertainty [ct/pix]')
#43: plt.ylabel('Number of Pixels')
#43: plt.xscale('log')
#43: plt.legend(frameon=False)

## Exposure Time Normalization ##
#44: m_171_norm = normalize_exposure(m_171_L15)

#45: print(m_171_L15.unit)
#45: print(m_171_norm.unit)

#46: print(m_171_L15.exposure_time)
#46: print(m_171_norm.exposure_time)

#47: m_171_norm = m_171_L15 / m_171_L15.exposure_time

#48: print(m_171_norm.unit)

#49: print(m_171_norm.exposure_time)

## Wavelength Response Functions ##
#50: c = Channel(m_335.wavelength)

#51: print(c.channel)
#51: print(c.telescope_number)

#52: r = c.wavelength_response()

#53: print(c.geometrical_collecting_area)
#53: print(c.primary_reflectance)
#53: print(c.focal_plane_filter_efficiency)
#53: print(c.contamination)
#53: print(c.quantum_efficiency)
#53: print(c.gain)

#54: r_time = c.wavelength_response(obstime=m_335.date)
#54: r_time_eve = c.wavelength_response(obstime=m_335.date, include_eve_correction=True)

#55: plt.plot(c.wavelength,r,label='Uncorrected')
#55: plt.plot(c.wavelength,r_time,label='Degradation')
#55: plt.plot(c.wavelength,r_time_eve,label='Degradation + EVE')
#55: plt.xlim([315,355])
#55: plt.ylim([0,0.03])
#55: plt.xlabel('$\lambda$ [Å]')
#55: plt.ylabel(f'$R(\lambda)$ [{r.unit.to_string(format="latex_inline")}]')
#55: plt.legend(loc=2, frameon=False)
