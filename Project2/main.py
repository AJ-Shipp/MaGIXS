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

import matplotlib.pyplot as plt

if __name__ == '__main__':

    passR = True
    ## Creating parameters for the shell
    bs1 = [0,0,0]               #[cm]
    focalLength = 200.0         #[cm]
    sLength = 30.0              #[cm]
    radii = [5.15100]
    angle = 0.00643732691573

    ## Initialization of the module, detector, source, and shell
    module = Module(base=[0, 0, 0], seglen=30.0, focal=200.0, radii=[5.151], angles=None,
                 conic=True, shield=True, core_radius=0)
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
    detector = Detector(center=[0,0,233.5], height=1.5, width=3, reso=[1024,2048])
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
    r = 200 #[cm]
    theta = 0/60/60                    #Deg/60/60 = x'
    phi = 0/60/60                      #Deg/60/60 = x'
    x_coord = r*np.sin(phi)*np.cos(theta)
    y_coord = r*np.sin(phi)*np.sin(theta)
    z_coord = r*np.cos(phi)
    z_ang = r*np.cos(theta)/r

    source = Source(center=(x_coord,x_coord,-z_coord),width=5.344122475091453*2, 
                    height=5.344122475091453*2, normal=(-x_coord,-x_coord,z_coord))
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
    shell = Shell(base=bs1, seglen=sLength, conic=True)
    """
    Parameters:
            base:    the center point of the wide end of the segment
            focal:   the focal length, measured from the center of the module
            seglen:  the axial length of each segment
            ang:     angle between the shell axis and the side of the front
                     segment
            r:       radius of the shell where the two segments meet    
    """
    
    # creating a 3D image of the shell to verify its creation
    fig2 = plt.figure(figsize=(5, 5))
    axes2 = get3dAxes(fig2)
    plt.axes(axes2)
    module.plot3D(axes2, 'b')
    source.plot3D(axes2, 'b')
    detector.plot3D(axes2, 'b')

    # generate 50000 rays at source
    rays = source.generateRays(module.targetFront, 50000)

    # pass rays through shell
    surfaces = shell.getSurfaces() # each shell has two segments

    modDims = module.getDims()
    modR_w = modDims[0]
    modR_s = modDims[1]
    modL = modDims[2]
    print(modDims[0])

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
            
            """
            If the ray's x coordinate is below 0.5*module's wide end radius, or if the y coordinate 
            is below 0*module's wide end radius, then the ray is removed  
            """
            #: if (ray.pos[0] < modR_w*0.5 or ray.pos[0] > modR_w) or (ray.pos[1] < modR_w*0 or ray.pos[1] > modR_w*(3**(1/2))/2): 
            #:     ray.dead = True

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

                    """
                    If the ray is not dead, then the rays with 2 bounces are plotted with green
                        and the rays with only 1 bounce (ghost rays) are plotted in red
                    """
                    if ray.dead == False:
                        if ray.bounces == 2:
                            ray.plot3D(axes2, 'g')
                        elif ray.bounces == 1:
                            ray.plot3D(axes2, 'r')

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
            if ray.num % 5000 == 0:
                print(ray.num, "rays passed")

    # catch rays at detector
    detector.catchRays(rays)

    ## plot detector pixels
    plot(detector)
    plt.gca().invert_yaxis()

    # create scatter plot
    detectorRays = detector.rays
    fig = plt.figure(figsize=(10,10), dpi=50) #Default values of 'figsize=(5,5), dpi=100'
    scatterHist(detectorRays, fig, binwidth=0.05) #binwidth = 1E? #-# 0.05 w/ default detector is wanted shape

    # show
    plt.show()
