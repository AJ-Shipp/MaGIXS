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

import matplotlib.pyplot as plt

if __name__ == '__main__':

    ## Creating parameters for the shell
    bs1 = [0,0,0]               #[cm]
    focalLength = 200.0         #[cm]
    sLength = 30.0              #[cm]
    radii = [5.15100]
    angle = 0.00643732691573

    ## Initialization of the module, detector, source, and shell
    module = Module(base=[0, 0, 0], seglen=30.0, focal=200.0, radii=[5.151], angles=None,
                 conic=True, shield=True, core_radius=None)
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
    detector = Detector(center=[2,2,235], height=10, width=10, reso=[512,512])
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
    source = Source()
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

    # generate 5000 rays at source
    rays = source.generateRays(module.targetFront, 5000)
    #: rays = source.generateRays(shell.targetFront, 5000)

    # pass rays through shell
    surfaces = shell.getSurfaces() # each shell has two segments
    for ray in rays:
        while True:
            #: ray.plot3D(axes2, 'k')              #Plots all rays that are created
            sol = None
            for surface in surfaces:
                
                # solve for intersection of ray with surface
                sol = surface.rayIntersect(ray)
                if sol is not None: break
            
            # if ray hits reflective surface
            if sol is not None:
                #: ray.plot3D(axes2, 'k')          #Plots the rays that will hit the module's initial position

                # update ray
                ray.pos = ray.getPoint(sol[2])
                ray.bounces += 1
                x = reflect(ray.ori,surface.getNormal(sol[0],sol[1]), None)
                # if reflected
                if x is not None:
                    ray.ori = x / norm(x) # update ori to unit vector reflection
                    ray.plot3D(axes2, 'g')      #Plots the rays' position after reflection on the module

                # otherwise, no reflection means ray is dead
                else:
                    ray.dead = True 
                    ray.plot3D(axes2, 'r')          #Plots the rays that will hit the module's initial position
                    break
                #: print(ray.ori)               #Prints the rays' origin
                #: print(ray.pos)               #Prints the rays' position
                
            else: break

    # catch rays at detector
    detector.catchRays(rays)

    ## plot detector pixels
    plot(detector)

    # create scatter plot
    detectorRays = detector.rays
    fig = plt.figure(dpi=100) #Default values of 'figsize=(5,5), dpi=100'
    scatterHist(detectorRays, fig, binwidth=0.05) #binwidth = 1E? #-# 0.05 w/ default detector is wanted shape

    # show
    plt.show()
