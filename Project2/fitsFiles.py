import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/throughFocusProgram/dataOutPut.csv"
fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus-1_mask_5sec.fits"
myData = "myData.fits"
inputImg = myData

img = np.zeros((1024,2048))

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
plt.imshow(image_data, cmap='gray')
plt.colorbar()
plt.show()

print(img)