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
    base1 = [0,1,0]
    focalLength = 200.0
    segmentLength = 30.0
    radii = [5.15100]
    angle = 0.00643732691573

    ## Initialization of the module, detector, source, and shell
    #: module = Module(seglen=segmentLength, focal=focalLength, radii=radii, conic=True, core_radius=2.856)
    detector = Detector()
    source = Source()
    shell = Shell(base=base1, seglen=segmentLength, conic=True)
    
    # creating a 3D image of the shell to verify its creation
    #: fig2 = plt.figure(figsize=(5, 5))
    #: axes2 = get3dAxes(fig2)
    #: plt.axes(axes2)
    #: shell.plot3D(axes2, 'b')

    # generate 5000 rays at source
    rays = source.generateRays(shell.targetFront, 5000)

    # pass rays through shell
    surfaces = shell.getSurfaces() # each shell has two segments
    for ray in rays:
        while True:
            
            sol = None
            for surface in surfaces:
                
                # solve for intersection of ray with surface
                sol = surface.rayIntersect(ray)
                if sol is not None: break
            
            # if ray hits reflective surface
            if sol is not None:
                
                # update ray
                ray.pos = ray.getPoint(sol[2])
                ray.bounces += 1
                x = reflect(ray.ori,surface.getNormal(sol[0],sol[1]), None)
                # if reflected
                if x is not None:
                    ray.ori = x / norm(x) # update ori to unit vector reflection
                # otherwise, no reflection means ray is dead
                else:
                    ray.dead = True 
                    break
                
            else: break

    # catch rays at detector
    detector.catchRays(rays)

    ## plot detector pixels
    plot(detector)

    # create scatter plot
    detectorRays = detector.rays
    fig = plt.figure(figsize=(5,5))
    scatterHist(detectorRays, fig, binwidth=0.001) #binwidth = 1E-11

    # show
    plt.show()
