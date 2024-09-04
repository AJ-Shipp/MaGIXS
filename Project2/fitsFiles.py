import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import matplotlib.cm as cm

#"C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus-5_mask_5sec.fits"
myData = "myData_M5_rot+shift.fits"
inputImg = myData
datRot = False
contour = True

img = np.zeros((2048,2048))
th = np.radians(25)
rotMtrx = np.array([[np.cos(th),np.sin(th)],[-np.sin(th),np.cos(th)]])

with open(dataFile, 'r') as csvfile:
    csvReader = csv.reader(csvfile)

    dataList = []
    for row in csvReader:
        dataList.append(row)

    for row in dataList:
        r = row
        xDiff = 0
        yDiff = 0
        if datRot == True:
            xDiff = 393
            yDiff = 549

            x = float(r[0])
            y = float(r[1])
            oldPair = np.array([[x],[y]])
            newPair = rotMtrx @ oldPair
        
            r[0] = int(np.floor(newPair[0][0]))
            r[1] = int(np.floor(newPair[1][0]))
        img[int(r[0])-xDiff][int(r[1])+yDiff] = float(r[2])

print(dataList[0:])
print("\n")

with fits.open(fits_image) as hdul:
    hdul.info()
    
fits.writeto(myData, data=img, overwrite = True)

image_data = fits.getdata(inputImg)
print("\n")
plt.imshow(image_data, cmap='jet', origin='lower')
plt.colorbar()
#: plt.xlim()
#: plt.ylim()

if contour == True:
    plt.contour(image_data, [1000])
    
plt.show()

print(img)