import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a

import aiapy
from aiapy.calibrate import (
    correct_degradation,
    degradation,
    estimate_error,
    fetch_spikes,
    register,
    respike,
    update_pointing,
)
from aiapy.calibrate.util import get_correction_table, get_error_table, get_pointing_table
from aiapy.psf import deconvolve, psf
from aiapy.response import Channel


import cupy

query = Fido.search(a.Time('2024-07-16 4:00:00', '2024-07-16 7:00:00'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample((1/12)*u.h))

print(query)

files = Fido.fetch(query)

aia_sequence = sunpy.map.Map(files, sequence=True)



firstImg = aia_sequence[0]
psf_171 = psf(firstImg.wavelength)

for img in aia_sequence:
    img = deconvolve(img, psf=psf_171)
    img = respike(img)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence[0])
norm = norm=ImageNormalize(vmin=0, vmax=3e3, stretch=SqrtStretch())
ani = aia_sequence.plot(axes=ax, norm=norm)
plt.show()

corner = SkyCoord(Tx=-150*u.arcsec, Ty=-750*u.arcsec, frame=aia_sequence[6].coordinate_frame)
cutout_map = aia_sequence[6].submap(corner, width=750*u.arcsec, height=750*u.arcsec)

fig = plt.figure()
ax = fig.add_subplot(projection=cutout_map)
cutout_map.plot(axes=ax)
plt.show()

with propagate_with_solar_surface():
    aia_sequence_aligned = sunpy.map.Map([m.reproject_to(cutout_map.wcs) for m in aia_sequence], sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_sequence_aligned[0])
ani = aia_sequence_aligned.plot(axes=ax, cmap='sdoaia94', norm=norm, interval=1000)

plt.show()
