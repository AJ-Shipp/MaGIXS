import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import matplotlib.cm as cm
from matplotlib import colors

#"C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/simCsvs/arc_m4.csv"
myData = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/simFits/myData_M4_rot+shift.fits"                #0
fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus-4_mask_5sec.fits"    #1
datRot = True
contour = True
dataSelect = 1
simCountMax = 55

if dataSelect == 0:
    inputImg = myData
    datRot = True
elif dataSelect == 1:
    inputImg = fits_image

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
image_data2 = fits.getdata(myData)

print("\n")

plot = plt.imshow(image_data, cmap='jet', origin='lower', norm='linear', vmin=0)
plt.xlim(1000,1060)
plt.ylim(960,1020)

if contour == True:
    plot.axes.set_title('SLTF(m4) and csSim(m4)',fontsize=20)
    cs = plt.contour(image_data2, levels=[simCountMax*1/5, simCountMax*2/5, simCountMax*3/5, simCountMax*4/5])
    cs.set_linewidth(1)
    cbc = plt.colorbar(cs, shrink=1)
    cb = plt.colorbar(plot)
    
    print(cbc.get_ticks())
    cbc.set_ticks(ticks=cbc.get_ticks(),labels=['20%','40%','60%','80%'])
    cbc.set_label('Percent of Max Intensity')
else:
    plot.axes.set_title('SLTF(m4)',fontsize=20)
    plot.axes.set_title('csSim(m4)',fontsize=20)
    cb = plt.colorbar(plot)

plt.show()

# print(image_data)