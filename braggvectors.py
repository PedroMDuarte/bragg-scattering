import numpy as np
import vec3
import pylab
from scipy import stats



# Q vector for bragg scattering
l671 = 671.
l1064 = 1064.

# Remember that for AFM there is a doubling of the unit cell, so the
# lattice spacing is lambda instead of lambda/2
Q = 2*np.pi/l1064 * vec3.vec3( -1., -1., 1.)  
Qunit = Q/abs(Q)

# Calculate angle for Bragg condition
# with respect to Q vector
braggTH = np.arccos( abs(Q) / 2. / (2.*np.pi/l671)  )  
print "Bragg angle wrt Q = ", braggTH * 180. / np.pi 

# Calculate angle for Bragg condition
# with respect to y axis, when coming from 
# under lattice beam 2. 
from scipy.optimize import root
def cond(x):
    return np.sin(x)-np.cos(x) + 3./2. * l671 / l1064
braggTH2 = (root(cond, 0.).x)[0]
print "Bragg angle wrt -y axis = ", braggTH2 * 180. / np.pi

# Incoming light vector
thi = np.pi/2 - braggTH2
phi = 90. * np.pi / 180.
kin = vec3.vec3()
kin.set_spherical( 2.*np.pi/l671, thi, phi )

# Default polarization of incoming light vector
kipol = [1.,0] 

# Unit vector that points prep to Bragg cone
kinperp = vec3.cross( kin, vec3.cross(Q,kin) ) 
kinperp = kinperp / abs(kinperp)

# Direction of A2 detector
kout = kin + Q  
a2 = kout / abs(kout)
kA2 = kout

# Unit vector perpendicular to plane of Q and A2
Qperp1 = vec3.cross( Q, a2 )
Qperp1 = Qperp1 / abs(Qperp1)

# Unit vector perpendicular to Q and in plane of kin,kout,Q
Qperp2 = vec3.cross( Q, Qperp1)
Qperp2 = Qperp2 / abs(Qperp2)
# Using Qunit and Qperp2 one can use the Bragg angle to 
# easily parameterize vectors near kin and kout

# Define direction of A1 camera
a1 = vec3.vec3()
a1.set_spherical( 1., np.pi/2., np.pi/3. ) 
kA1 = a1*abs(kin)



# The main vectors defined so far are plotted on the sphere here
import bsphere
b = bsphere.Bloch()
origin = vec3.vec3()
b.add_arrow( origin, Q/abs(kin) , 'blue')
b.add_arrow( -kin/abs(kin), origin, 'red')
b.add_arrow( origin, kout/abs(kout), 'red')
b.add_arrow( origin, a1, 'black')
#b.add_arrow( origin, Qperp2, 'green')
#b.show()

# Here two functions are defined that allow getting kin 
# and kout vectors as a function of their angle measured
# from the nominal bragg angle (useful for rocking curve)

def kinput( angle ):  # angle is in mrad
    kia =  -Qunit*np.cos(braggTH + angle/1000.) - Qperp2*np.sin(braggTH + angle/1000.)
    kia = abs(kin) * kia
    b.add_points( (-1*kia/abs(kia)).tolist()  )
    return kia

def koutput( angle ):  # angle is in mrad
    kfa =  Qunit*np.cos(braggTH + angle/1000.) - Qperp2*np.sin(braggTH + angle/1000.) 
    kfa = abs(kin) * kfa
    b.add_points( (kfa/abs(kfa)).tolist()  )
    return kfa 

