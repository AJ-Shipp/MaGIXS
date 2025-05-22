import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import matplotlib.cm as cm
from matplotlib import colors

#"C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project3/simCsvs/Det20x10k_Sm5upTest3.csv"
myData = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project3/simFits/myData_Det20x10k_Sm5upTest3.csv.fits"                #0
# fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus+6_mask_5sec.fits"    #1
datRot = False
contour = True
dataSelect = 0
simCountMax = 50

if dataSelect == 0:
    inputImg = myData
    datRot = False
# elif dataSelect == 1:
    # inputImg = fits_image

img = np.zeros((2048*2,2048*2))
th = np.radians(25)
rotMtrx = np.array([[np.cos(th),np.sin(th)],[-np.sin(th),np.cos(th)]])

def fmt(x):
    x = x/simCountMax*100
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

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
        img[int(r[1])+xDiff][int(r[0])-xDiff] = float(r[2])

# print(dataList[0:])
# print("\n")

# with fits.open(fits_image) as hdul:
    # hdul.info()
    
fits.writeto(myData, data=img, overwrite = True)
# image_data = fits.getdata(inputImg)
image_data = fits.getdata(myData)

print("\n")

plot = plt.imshow(image_data, cmap='binary', origin='lower', norm='linear', vmin=0)
# plt.xlim(1000+35,1060+5)
# plt.ylim(960-30,1020-60)

if contour == True:
    plot.axes.set_title('Det12x6k_Szero',fontsize=20)
    cv = plt.contour(image_data, levels=[0])
    cb = plt.colorbar(plot, shrink=.6)
    # cs = plt.contour(image_data, levels=[-500, simCountMax*1/10, simCountMax*2/10, simCountMax*2/5, simCountMax*3/5, simCountMax*4/5, 500], 
                    #  colors=('#ffffff', '#ffe000','#ffbd00', '#ffb600', '#ff7f0e', '#d62728', '#ffffff'))
    # cs.set_linewidth(1.4)
    # cbc = plt.colorbar(cs, shrink=1)
    # cb = plt.colorbar(plot)
    # cbc.add_lines(levels=[-100, simCountMax*1/10, simCountMax*2/10, simCountMax*2/5, simCountMax*3/5, simCountMax*4/5, 100], 
    #               colors=('#ffffff', '#ffe000','#ffbd00', '#ffb600', '#ff7f0e', '#d62728', '#ffffff'), linewidths=6)

    # print(cbc.get_ticks())
    # cbc.set_ticks(ticks=cbc.get_ticks(),labels=['','10%','20%','40%','60%','80%',''], fontsize=12)
    # cbc.set_label('Simulation Percent of Max Intensity', fontsize=15)
    
else:
    plot.axes.set_title('SLTF(m4)',fontsize=20)
    plot.axes.set_title('csSim(m4)',fontsize=20)
    cb = plt.colorbar(plot)

plt.show()

# print(image_data)

"""
import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import matplotlib.cm as cm
from matplotlib import colors

#"C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/arc_m5.csv"
dataFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project3/simCsvs/Det12x6k_Sm6up.csv"
myData = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project3/simFits/myData_Det12x6k_Sm6up.fits"                #0
fits_image = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project2/fitsExamples/proc_focus+6_mask_5sec.fits"    #1
datRot = False
contour = True
dataSelect = 0
simCountMax = 50

if dataSelect == 0:
    inputImg = myData
    datRot = False
elif dataSelect == 1:
    inputImg = fits_image

img = np.zeros((2048*3,2048*6))
th = np.radians(25)
rotMtrx = np.array([[np.cos(th),np.sin(th)],[-np.sin(th),np.cos(th)]])

def fmt(x):
    x = x/simCountMax*100
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

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
        img[int(r[1])+xDiff][int(r[0])-xDiff] = float(r[2])

# print(dataList[0:])
# print("\n")

with fits.open(fits_image) as hdul:
    hdul.info()
    
fits.writeto(myData, data=img, overwrite = True)
# image_data = fits.getdata(inputImg)
image_data = fits.getdata(myData)

print("\n")

plot = plt.imshow(image_data, cmap='binary', origin='lower', norm='linear', vmin=0)
# plt.xlim(1000+35,1060+5)
# plt.ylim(960-30,1020-60)

if contour == True:
    plot.axes.set_title('Det12x6k_Szero',fontsize=20)
    cv = plt.contour(image_data, levels=[0])
    cb = plt.colorbar(plot, shrink=.6)
    # cs = plt.contour(image_data, levels=[-500, simCountMax*1/10, simCountMax*2/10, simCountMax*2/5, simCountMax*3/5, simCountMax*4/5, 500], 
                    #  colors=('#ffffff', '#ffe000','#ffbd00', '#ffb600', '#ff7f0e', '#d62728', '#ffffff'))
    # cs.set_linewidth(1.4)
    # cbc = plt.colorbar(cs, shrink=1)
    # cb = plt.colorbar(plot)
    # cbc.add_lines(levels=[-100, simCountMax*1/10, simCountMax*2/10, simCountMax*2/5, simCountMax*3/5, simCountMax*4/5, 100], 
    #               colors=('#ffffff', '#ffe000','#ffbd00', '#ffb600', '#ff7f0e', '#d62728', '#ffffff'), linewidths=6)

    # print(cbc.get_ticks())
    # cbc.set_ticks(ticks=cbc.get_ticks(),labels=['','10%','20%','40%','60%','80%',''], fontsize=12)
    # cbc.set_label('Simulation Percent of Max Intensity', fontsize=15)
    
else:
    plot.axes.set_title('SLTF(m4)',fontsize=20)
    plot.axes.set_title('csSim(m4)',fontsize=20)
    cb = plt.colorbar(plot)

plt.show()

# print(image_data)
"""