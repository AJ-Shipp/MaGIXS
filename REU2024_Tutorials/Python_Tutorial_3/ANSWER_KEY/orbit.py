import numpy as np
import matplotlib.pyplot as plt
import math as m

def orbit_func(xx, yy, vx, vy):
	G = 6.6742E-11 # m^3 kg^-1 s^2
	M_s = 1.9889E30 # kg
	r = np.sqrt(xx**2+yy**2)
	F_factor = -G*M_s / r**3
	ax = F_factor*xx; ay = F_factor*yy
	return ax, ay

au = 1.49598e11 # m
day = 60*60*24 # s
dt = 1 * day

# part a below
G = 6.6742E-11 # m^3 kg^-1 s^2
M_s = 1.9889E30 # kg
semi_major = 50.1 * au # convert
radius = 0.2 * au #convert
velocity = m.sqrt(G*M_s*(2/radius - 1/semi_major))
print(velocity, "m/s is the velocity of the comet at Perihelion")
#94104.26949202036 m/s is the velocity of the comet at Perihelion

#part b below
marsradius = 1.5 * au
marsv = m.sqrt((G*M_s)/marsradius)
print("The minimum launch speed needed to get to mars is", velocity - marsv, "m/s")
#The minimum launch speed needed to get to mars is 69782.36610748505

xx = 1 * au
yy = 0
vx = 0
vy = 29800

#part c
escape = m.sqrt((2*G*M_s)/ xx)
print("The escape speed for the solar system is", escape, "m/s")
#The escape speed for the solar system is 42126.77239879659 m/s

npoints = 1000
x_vect = np.zeros(npoints)
y_vect = np.zeros(npoints)

x_vect[0] = xx/au
y_vect[0] = yy/au
for k in range(1,npoints-1):
	ax, ay = orbit_func(xx, yy, vx, vy)
	xx = xx + vx*dt
	yy = yy + vy*dt
	vx = vx + ax*dt
	vy = vy + ay*dt
	x_vect[k] = xx/au
	y_vect[k] = yy/au

plt.figure()
plt.plot(x_vect, y_vect)
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
#plt.show()

