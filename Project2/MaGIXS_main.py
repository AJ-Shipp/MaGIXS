#!/usr/bin/env python
# June 2014, @milo & @Steven
''' Example 2 using the geometry of the Hyp and Par. This should show
perfect focusing with all rays falling in a point.'''

from foxsisim.module import Module
from foxsisim.shell import Shell
from foxsisim.detector import Detector
from foxsisim.source import Source
from foxsisim.mymath import reflect
from numpy.linalg import norm
from foxsisim.plotting import scatterHist, get3dAxes, plot
import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt

def distance(x,y):
        z = np.sqrt(np.square(x) + np.square(y))
        return z

if __name__ == '__main__':

    """
    Different Settings/Add-Ons

    passR = Passes the created rays through the module
    spotW = Function to find the width of the focus spot
    plot3D = Plots a 3D image of the source, module, rays, and detector
    plotDetector = Plots the detector
    plotScatHist = Plots the scatter histogram
    numRays = How many rays you want to create
    arcminOff = How far off axis you'd like the source to be
    arcminDiag = If you want the x & y axes to be equally off axis
    """
    passR = True
    spotW = False
    plot3D = True
    plotDetector = True
    plotScatHist = True
    numRays = 500000
    arcminOff = -5
    arcminDiag = False
    scatHistSize = 10
    binW = 0.01
    blockerSegment = True
    segmentPercents = True
    bSegPos = 0
    bSegAng = 17
    allRays2File = True
    fitsBool = True
    dataOutputFile = "C:/Users/antho/OneDrive/Documents/GitHub/MaGIXS/Project3/simCsvs/Det20x10k_Sm5upTest3.csv"
    detPosZ = 0.0              #[cm]

    ## Creating parameters for the shell
    bs1 = [0,0,0]               #[c
    focalLength = 109.01        #[cm]
    sLength = 12.50             #[cm]
    radii = [7.3061]
    angle = 0.00643732691573

    ## Initialization of the module, detector, source, and shell
    module = Module(base=[0, 0, 0], seglen=sLength, focal=focalLength, radii=radii, angles=None,
                 conic=False, shield=True, core_radius=0)
    """
    Parameters:
            base:       the center point of the wide end of the segment
            seglen:     the axial length of each segment
            focal:      the focal length, measured from the center of the module
            radii:      a list of radii, one for each shell from biggest to
                        smallest
            angles:     optional parameter to overwrite the shell angles computed
                        by constructor
            conic:      if True, use a conic approximation to the Wolter-I parabola-hyperbola.
            shield:       if True, create a shield at the core of the optic to block straight-thru
                        flux.
            core_radius If shield is True then use this value for the shield radius
    """
    
    #Setting variables for the module dimensions, prints wide section's radius
    modDims = module.getDims()
    modR_w = modDims[0]
    modR_s = modDims[1]
    modL = modDims[2]
    print(modDims[:])
    
    detector = Detector(center=[-5.736,0,120.84+detPosZ], normal=[0,0,1], height=1.3824, width=2.7648, reso=[1024,2048])
    """
    Parameters:
            center:    the center location of the detector
            width:     the width of the projection rectangle
            height:    the height of the projection rectangle
            normal:    direction the projection rectangle is facing
            reso:      resolution of pixels (W,H)
            pixels:    optional numpy array to hold pixel colors (W x H x 3)
            freqs:     optional array to count the number of times a pixel is
                       hit (W x H)
    """

    ## Creating the parameters for the source

    ##Spherical (r, theta, phi) to Cartesian
    ## r = c
    r = 50 #[cm]
    theta = arcminOff/60                      #Deg/60 = x'
    theta2 = 0/60                             #Deg/60 = x'
    phi = 0                                   #Deg/60 = x'
    phi2 = 0                                  #Deg/60 = x'
    x_coord = r*np.sin(theta)*np.cos(phi)
    y_coord = r*np.sin(theta2)*np.cos(phi2)
    # if arcminDiag == True:
    #     y_coord = x_coord
    # else:
    #     y_coord = 0
    z_coord = r*np.cos(theta)
    z_ang = r*np.cos(theta)/r

    print("[%.14f,%.14f,%.14f]"%(x_coord,y_coord,z_coord))

    source = Source(center=(x_coord,y_coord,-z_coord),width=r*np.tan(30/60), 
                    height=r*np.tan(30/60), normal=(-x_coord,-y_coord,z_coord))
    """
    Parameters:
            center:    the center location of the source
            width:     the width of the projection rectangle (atinf or nonpoint)
            height:    the height of the projection rectangle (atinf or nonpoint)
            normal:    direction the projection rectangle is facing
            type:      'atinf', 'point', or 'nonpoint'
            color:     color of projected rays
            pixels:    optional numpy array of pixel colors (W x H x 3)
            spectrum:  optional numpy array (2xN) of energy spectrum
            tag : 'Source' Tag all rays with 'Source', place where the came from.
    """
    
    shell = Shell(base=bs1, seglen=sLength, conic=False)
    """
    Parameters:
            base:    the center point of the wide end of the segment
            focal:   the focal length, measured from the center of the module
            seglen:  the axial length of each segment
            ang:     angle between the shell axis and the side of the front
                     segment
            r:       radius of the shell where the two segments meet    
    """
    
    if plot3D == True:
        # creating a 3D image of the shell to verify its creation
        fig3D = plt.figure(figsize=(5, 5))
        axes3D = get3dAxes(fig3D)
        plt.axes(axes3D)
        module.plot3D(axes3D, 'b')
        source.plot3D(axes3D, 'b')
        detector.plot3D(axes3D, 'b')

    # generate 50000 rays at source
    rays = source.generateRays(module.targetFront, numRays)

    # pass rays through shell
    surfaces = shell.getSurfaces() # each shell has two segments

    if passR == True:
        '''
        Takes an array of rays and passes them through the front end of
        the module.
        '''
        # print('Module: passing ',len(rays),' rays')
        robust = True

        # get all module surfaces
        allSurfaces = module.getSurfaces()
        allSurfaces.remove(module.coreFaces[0])  # we'll test these seperately
        allSurfaces.remove(module.coreFaces[1])

        # create regions consisting of adjacent shells
        regions = [None for shell in module.shells]
        for i, shell in enumerate(module.shells):
            # innermost shell
            if i == len(module.shells) - 1:
                regions[i] = shell.getSurfaces()
                # regions[i].append(self.core)
            else:
                # outer shell (reflective side facing region)
                regions[i] = shell.getSurfaces()
                # nested shell (non reflective)
                regions[i].extend(module.shells[i + 1].getSurfaces())

        for ray in rays:

            # move ray to the front of the optics
            ray.moveToZ(module.coreFaces[0].center[2])

            # reset surfaces
            surfaces = [s for s in allSurfaces]
            firstBounce = True  # used for optimization

            while True:
                
                # find nearest ray intersection
                bestDist = None
                bestSol = None
                bestSurf = None
                for surface in surfaces:

                    sol = surface.rayIntersect(ray)
                    if sol is not None:
                        dist = norm(ray.getPoint(sol[2]) - ray.pos)
                        if bestDist is None or dist < bestDist:
                            bestDist = dist
                            bestSol = sol
                            bestSurf = surface

                # if a closest intersection was found
                if bestSol is not None:

                    # update ray
                    ray.pos = ray.getPoint(bestSol[2])
                    ray.hist.append(ray.pos)
                    ray.bounces += 1
                    ray.update_tag(bestSurf.tag)
                    #print("%i ray bounce number %i" % (ray.num, ray.bounces))

                    x = reflect(ray.ori,
                                bestSurf.getNormal(bestSol[0], bestSol[1]),
                                ray.energy)
                    # print('x = ',x)
                    # if reflected
                    if x is not None:
                        # update ori to unit vector reflection
                        ray.ori = x / norm(x)
                    # otherwise, no reflection means ray is dead
                    else:
                        ray.dead = True
                        # print("%i ray killed by reflect" % ray.num)
                        break

                    # if plot3D == True:
                    #     """
                    #     If the ray is not dead, then the rays with 2 bounces are plotted with green
                    #         and the rays with only 1 bounce (ghost rays) are plotted in red
                    #     """
                    #     if ray.dead == False:
                    #         if ray.bounces == 2:
                    #             ray.plot3D(axes3D, 'g')
                    #         elif ray.bounces == 1:
                    #             ray.plot3D(axes3D, 'r')

                    # knowing the surface it has just hit, we can
                    # narrow down the number of surface to test

                    # remove shells the ray cannot even 'see'
                    if firstBounce:
                        firstBounce = False
                        for region in regions:
                            if bestSurf is region[0] or bestSurf is region[1]:
                                surfaces = [s for s in region]
                                break

                    # assuming each segment can be hit no more than once
                    # eliminate the surface from our list
                    if not robust:
                        surfaces.remove(bestSurf)

                # if no intersection, ray can exit module
                else:
                    break
                    # print(ray.hist)
                    # print(ray.des)

                """
                If the ray's x coordinate is below 0.5*module's wide end radius, or if the y coordinate 
                is below 0*module's wide end radius, then the ray is removed  
                """
                if blockerSegment == True:
                    bSeg = bSegPos + bSegAng 
                    if (
                        (ray.pos[0] < 0 or (ray.pos[0] > modR_w)) or       #x-axis
                        (ray.pos[1] < modR_w*-np.sin(np.radians(bSegAng)) or ray.pos[1] > modR_w*(np.sin(np.radians(bSegAng)))) #y-axis
                        ): 
                        ray.dead = True
                
                if plot3D == True:
                        """
                        If the ray is not dead, then the rays with 2 bounces are plotted with green
                            and the rays with only 1 bounce (ghost rays) are plotted in red
                        """
                        if ray.dead == False:
                            if ray.bounces == 2:
                                ray.plot3D(axes3D, 'g')
                            elif ray.bounces == 1:
                                ray.plot3D(axes3D, 'r')

            sol = module.coreFaces[1].rayIntersect(ray)
            if sol is not None:
                print("ray hit rear blocker")
                ray.pos = ray.getPoint(sol[2])
                #ray.bounces += 1
                ray.dead = True
                ray.des = ray.pos
                ray.hist.append(ray.pos)
                ray.update_tag(module.coreFaces[1].tag)
                continue
            else:
                ray.moveToZ(module.coreFaces[1].center[2])
                #ray.hist.append(ray.pos)

            if ray.bounces == 0:
                ray.dead = True
            if ray.num % 25000 == 0:
                print(ray.num, "rays passed")

        # catch rays at detector
        detector.catchRays(rays)

    if plotDetector == True:
        ## plot detector pixels
        plot(detector)
        plt.gca().invert_yaxis()

    if plotScatHist == True:
        # create scatter plot
        detectorRays = detector.rays
        fig = plt.figure(figsize=(scatHistSize,scatHistSize), dpi=50) #Default values of 'figsize=(5,5), dpi=100'
        scatterHist(detectorRays, fig, binwidth=binW) #binwidth = 1E-7 #-# 0.05 w/ default detector is wanted shape

    if spotW == True:
    
        k=0
        counter = int()
        distMax = 0
        #Writing to outputA File: f = open("outputA", "w")
        for i in rays:
            if rays[k].des[2] != 0.:
                counter += 1
                dist = distance(rays[k].des[0], rays[k].des[1])
                if dist > distMax:
                        distMax = dist
                #Writing to outputA File: f.write(str(distMin))
                #Writing to outputA File: f.write(", ")
                #Writing to outputA File: f.write(str(distMax))
                #Writing to outputA File: f.write("\n")
                print(counter, distMax)
            k += 1
        print("The radius of the spread is", distMax, "cm")

    if fitsBool == True:
        #.fits Files
        up1 = 0
        dataTemps = []
        dataPairs = []
        dataFinal = []

        for ray in rays:
            if rays[up1].des[2] != 0:
                rays[up1].des[0] = np.floor(((rays[up1].des[0] + 0.6912*2) / 0.00135))
                rays[up1].des[1] = np.floor(((rays[up1].des[1] + 1.3824) / 0.00135))
                duo = rays[up1].des[0], rays[up1].des[1]
                dataTemps.append(duo)

                #   Should append all rays' final destinations into a list with the third value in the tuple being equal to the 
                # number of times that specific final position occurs in the original list of destinations
                """
                # initializing the titles and rows list
                fields = []
                rows = []
                data_list = []
                pairs = []
                counters = []
                checked = []

                # reading csv file
                with open(filename, 'r') as csvfile:
                    # creating a csv reader object
                    csvreader = csv.reader(csvfile)

                    # extracting field names through first row
                    fields = next(csvreader)

                    # Iterate through each row in the CSV file
                    for row in csvreader:
                        currentX = row[5]
                        currentY = row[6]
                        pair = (currentX, currentY)
                        pairs.append(pair)

                    # get total number of rows
                    print("Total no. of rows: %d" % (csvreader.line_num))

                    for data in pairs:
                        need = data
                        print(need)
                        num = pairs.count(need) 
                        insert = need, num
                        if checked.count(insert) == 0:
                            checked.append(insert)
                """
                # print(rays[up1].des[0], rays[up1].des[1])
            up1 += 1

        for data in dataTemps:
            need = data
            num = dataTemps.count(need)
            insert = need, num
            if dataPairs.count(insert) == 0:
                dataPairs.append(insert)

        for data in dataPairs:
            temp = data
            stuff = temp[0][0], temp[0][1], temp[1]
            dataFinal.append(stuff)
        # print(dataFinal)

    if allRays2File == True:

        num = 0
        f = open(dataOutputFile, "a")

        for data in dataFinal:
            rowData = data
            print(rowData)
            f.write(str(int(rowData[0])))
            f.write(", ")
            f.write(str(int(rowData[1])))
            f.write(", ")
            f.write(str(int(rowData[2])))
            f.write("\n")
            num += 1
        
        f.close()

    if segmentPercents == True:
        count = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0

        for ray in rays:
            if (
                            (ray.pos[0] > 0 and (ray.pos[0] < modR_w)) and       #x-axis
                            (ray.pos[1] > modR_w*-np.sin(np.radians(bSegAng)) and ray.pos[1] < modR_w*(np.sin(np.radians(bSegAng)))) #y-axis
                            ): 
                if 'Hy' in ray.tag or 'Pa' in ray.tag:
                    count += 1
                if 'Hy' in ray.tag and 'Pa' not in ray.tag:
                    count2 += 1
                if 'Pa' in ray.tag and 'Hy' not in ray.tag:
                    count3 += 1

            if 'D' in ray.tag:
                if 'D' in ray.tag: 
                    if 'Hy' in ray.tag or 'Pa' in ray.tag:
                        count4 += 1
                    if 'Hy' in ray.tag and 'Pa' not in ray.tag:
                        count5 += 1
                    if 'Pa' in ray.tag and 'Hy' not in ray.tag:
                        count6 += 1
        
        if count > 0:
            print("\nSegment Percentages\n"
                "-------------------\n"
                "Overall number of individual ray collisions: {} \n"
                "{:.2f}% Paraboloid (n = {}) \n"
                "{:.2f}% Hyperboloid (n = {}) \n"
                .format(count,count2/count*100,count2,count3/count*100,count3)
                )
        
        if count4 > 0:         #:If any rays hit the detector, the collision percentages for these rays are shown 
            print("\nDetector Hit Percentages\n"
                    "-------------------\n"
                    "Overall number of individual ray's detected: {} \n"
                    "{:.2f}% Paraboloid Only (n = {}) \n"
                    "{:.2f}% Hyperboloid Only (n = {}) \n"
                    .format(count4,count5/count4*100,count5,count6/count4*100,count6)
                    )

    # print("Source Centriod Position is x={:.4f}, y={:.4f}, z={:.4f}\n"
    #       "Source Centriod Normal is x={:.4f}, y={:.4f}, z={:.4f}\n"
    #       "Source Angle is {:d} arcminutes off-axis\n"
    #       "Module Centriod Position is x={:.4f}, y={:.4f}, z={:.4f}\n"
    #       "Detector Centriod Position is x={:.4f}, y={:.4f}, z={:.4f}\n"
    #       "Detector Centriod Normal is x={:.4f}, y={:.4f}, z={:.4f}\n"
    #       "Detector Size = [1.3824,2.7648*2] with Resolution = [1024,2048]\n"          
    #       .format(x_coord,y_coord,-z_coord,-x_coord,-y_coord,z_coord,arcminOff,0,0,0,-5.736,0,121.51,0,0,1)
    #       )

    # for ray in rays:
    #     if 'D' in ray.tag:
    #         print("{}".format(ray.pos))

    # show
    plt.show()
