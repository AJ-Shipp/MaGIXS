import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

#"C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus-5_mask_5sec.fits"
myData = "myData+5.fits"
inputImg = myData

img = np.zeros((2048,2048))

with open(dataFile, 'r') as csvfile:
    csvReader = csv.reader(csvfile)

    dataList = []
    for row in csvReader:
        dataList.append(row)
    for row in dataList:
        r = row
        img[int(r[0])][int(r[1])] = int(r[2])        

print(dataList[0:5])
print("\n")

with fits.open(fits_image) as hdul:
    hdul.info()
    
fits.writeto(myData, data=img, overwrite = True)

image_data = fits.getdata(inputImg)
print("\n")
plt.imshow(image_data, cmap='gray', origin='lower')
plt.colorbar()
#: plt.xlim()
#: plt.ylim()
plt.contour(image_data)
plt.show()

print(img)