import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus-1_mask_5sec.fits"

with fits.open(fits_image) as hdul:
    hdul.info()
    print(hdul._data)
    data = hdul._data
data[30:40, 10:20]


