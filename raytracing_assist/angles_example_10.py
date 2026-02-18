from foxsisim.module import Module
from foxsisim.detector import Detector
from foxsisim.source import Source
import numpy as np


if __name__ == '__main__':

    # -----------------------------
    # Module parameters (MaGIXS)
    # -----------------------------
    focalLength = 109.01
    segmentLength = 12.50
    radii = [7.3061]

    module = Module(base=[0,0,0],
                    seglen=segmentLength,
                    focal=focalLength,
                    radii=radii,
                    conic=True,
                    core_radius=0)

    detector = Detector(center=[0,0,focalLength],
                        height=10,
                        width=10,
                        reso=[1024,1024],
                        tag='Det')

    moduleDimensions = module.getDims()
    moduleRadius_wide = moduleDimensions[0]

    source = Source(width=moduleRadius_wide*2,
                    height=moduleRadius_wide*2,
                    tag='Source')

    n_rays = 5000
    rays = source.generateRays(module.targetFront, n_rays)

    module.passRays(rays, robust=True)
    detector.catchRays(rays)

    # -----------------------------------------
    # Compute grazing angle for EACH ray
    # -----------------------------------------

    ray_index = 1

for ray in rays:

    if ray.bounces == 2 and len(ray.hist) >= 4:

        print("\n======================================")
        print(f"Ray {ray_index}")

        # ---- Hit positions ----
        P0 = np.array(ray.hist[0])
        P1 = np.array(ray.hist[1])   # First mirror hit
        P2 = np.array(ray.hist[2])   # Second mirror hit
        P3 = np.array(ray.hist[3])   # Detector hit

        print("P0 (source/front):   ", P0)
        print("P1 (1st mirror hit): ", P1)
        print("P2 (2nd mirror hit): ", P2)
        print("P3 (detector hit):   ", P3)

        # ---- First reflection ----
        v_in1  = P1 - P0
        v_out1 = P2 - P1

        v_in1  /= np.linalg.norm(v_in1)
        v_out1 /= np.linalg.norm(v_out1)

        dev1 = np.arccos(np.clip(np.dot(v_in1, v_out1), -1.0, 1.0))
        alpha1 = np.degrees(dev1 / 2.0)

        print("First mirror grazing angle (deg):", alpha1)

        # ---- Second reflection ----
        v_in2  = P2 - P1
        v_out2 = P3 - P2

        v_in2  /= np.linalg.norm(v_in2)
        v_out2 /= np.linalg.norm(v_out2)

        dev2 = np.arccos(np.clip(np.dot(v_in2, v_out2), -1.0, 1.0))
        alpha2 = np.degrees(dev2 / 2.0)

        print("Second mirror grazing angle (deg):", alpha2)

        ray_index += 1

