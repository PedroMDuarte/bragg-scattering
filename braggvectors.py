import numpy as np
import vec3
import pylab
from scipy import stats
import scipy



l671 = 671.
l1064 = 1064.

# Coordinate system:
# input Bragg light for HHH propagates almost along +Y 
# +Z is up


# Q vector for bragg scattering HHH
# Remember that for AFM there is a doubling of the unit cell, so the
# lattice spacing is lambda instead of lambda/2
Q = 2*np.pi/l1064 * vec3.vec3( -1., -1., 1.)  
Qunit = Q/abs(Q)

# Calculate angle for HHH Bragg conditiion
# with respect to Q vector
braggTH = np.arccos( abs(Q) / 2. / (2.*np.pi/l671)  )  
print "HHH Bragg angle wrt Q = ", braggTH * 180. / np.pi 

# Calculate angle for HHH Bragg condition
# with respect to y axis, when coming from 
# under lattice beam 2. 
from scipy.optimize import fsolve
def cond(x):
    return np.sin(x)-np.cos(x) + 3./2. * l671 / l1064
braggTH2 = fsolve(cond, 0.)
print "HHH Bragg angle wrt -y axis = ", braggTH2 * 180. / np.pi


# Q for 100 scattering
Q100 =  2*np.pi / (l1064/2) * vec3.vec3( 0., +1, 0.)
Q100unit = Q100/abs(Q100)

# Calculate angle for 100 Bragg condition
# with respect to Q vector
braggTH100 = np.arccos( abs(Q100) / 2. / (2.*np.pi/l671) ) 
print
print "100 Bragg angle wrt Q = ", braggTH100 * 180. / np.pi 


# Incoming and outgoing light vector for 100 
kin100 = vec3.vec3()
kin100.set_spherical( 2.*np.pi/l671, np.pi/2 - braggTH100, 3* np.pi / 2)
kout100 = vec3.vec3()
kout100.set_spherical( 2.*np.pi/l671, np.pi/2 - braggTH100, np.pi / 2)
kMANTA = kout100



# Incoming light vector HHH 
thi = np.pi/2 - braggTH2
phi = 90. * np.pi / 180.
kin = vec3.vec3()
kin.set_spherical( 2.*np.pi/l671, thi, phi )

# Default polarization of incoming light vector
kipol = [1.,0] 

# Unit vector that points perp to Bragg cone
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

#print kA1/abs(kA1)
#print kA2/abs(kA2)
#print kin/abs(kin)

# Angle between A1 and A2 
thetaA1A2 =  180. / np.pi * np.arccos(kA1*kA2 / abs(kA1) / abs(kA2))
print "Angle between A1 and A2 = %.2f" % thetaA1A2




# Here two functions are defined that allow getting kin 
# and kout vectors as a function of their angle measured
# from the nominal bragg angle (useful for rocking curve)

def kinput( angle ):  # angle is in mrad
    kia =  -Qunit*np.cos(braggTH + angle/1000.) - Qperp2*np.sin(braggTH + angle/1000.)
    kia = abs(kin) * kia
    #b.add_points( (-1*kia/abs(kia)).tolist()  )
    return kia

def koutput( angle ):  # angle is in mrad
    kfa =  Qunit*np.cos(braggTH + angle/1000.) - Qperp2*np.sin(braggTH + angle/1000.) 
    kfa = abs(kin) * kfa
    #b.add_points( (kfa/abs(kfa)).tolist()  )
    return kfa

# I can define a variation of kinput in the plane of the 
# chamber
def kinchamber( phi ):
    k = vec3.vec3()
    the_ = np.pi/2 
    phi_ = np.pi/2. + np.pi* phi/180.
    k.set_spherical( 2.*np.pi/l671, the_, phi_ )
    return k
     

# Here I can define a kinput angle by giving the polar and 
# azimuthal angles of the window that the Bragg beam is comming in
# Using this definition, our nominal Bragg position is
# polar = np.pi/2 + 3.0degrees
# azim  = 0.
def kinsph(theta, phi):
    k = vec3.vec3()
    the_ = np.pi/2 - ( np.pi*theta/180. - np.pi/2.)
    phi_ = np.pi/2. + np.pi* phi/180.
    k.set_spherical( 2.*np.pi/l671, the_, phi_ )
    return k

ksquad =[]
ksquad.append ( kinsph(91., -75.) )
ksquad.append ( kinsph(90.,-14.) )
ksquad.append ( kinsph(88., -4.) )
ksquad.append ( kinsph(93., 0.) )
ksquad.append ( kinsph(91., 15.) )
ksquad.append ( kinsph(90., 34.) )
# Extra points for plot 
#ksquad.append ( kinsph(90.,-2.))
#ksquad.append ( kinsph(90., 2.))
#ksquad.append ( kinsph(90., 4.))
#ksquad.append ( kinsph(90., 6.))

ksquadth = np.array( [-75.,-14.,-4.,0.,15.,34.] \
                  # + [-2.,2.,4.,6.] \
                   ) *1. * np.pi/180.

print 
print "Difference wrt Bragg Q,  |Q-K| * l1064 / (2)"
#print "  Nominal K =", 1./(abs(kout-kin - Q)*l1064/(4*np.pi))
print " -250mrad K =", abs(kout-kinsph(90.,-14.)-Q)*l1064/2
print "Same port K =", abs(kout-kinsph(88.,-4.)-Q)*l1064/2
print " +250mrad K =", abs(kout-kinsph(91.,15.)-Q)*l1064/2
print " +500mrad K =", abs(kout-kinsph(90.,34.)-Q)*l1064/2

# Here I can define four koutput vectors at the four quadrants
# of the Bragg lens 
kout_quadrants = []
kout_quadrants.append ( koutput(+50.))
kout_quadrants.append ( koutput(-50.))

def koutput_perp(angle):
    kfa = kout + abs(kout)*np.sin(angle/1000.)* Qperp1
    kfa = abs(kin) *  kfa / abs(kfa)
    return kfa 

kout_quadrants.append( koutput_perp(+50.))
kout_quadrants.append( koutput_perp(-50.))

# Next, a function is defined that creates a list of 
# koutput vectors in a circular solid angle, with a given
# diameter in mrad, centered around a given angle
def k2aperture( angle, aperture): 
    step = 20. #space in mrad of points in list
    nstep = np.ceil(aperture/step)
    avals = np.linspace( -nstep*step, nstep*step, 2*nstep+1,  endpoint=True)
    arrays = [avals,avals]
    arr = np.empty([avals.size]*2+[2])
    for i, a in enumerate( np.ix_(*arrays) ):
        arr[...,i] = a
    Kset = arr.reshape(-1,2)
    #for K in kset:
        

    #return arr.reshape(-1,2)
   
   
    print [avals.size]*2 + [2]
 
    print avals.shape
    print arr.shape
    
    #print arr
    #print avals 


###### VISUALIZATION #####
    
# The main vectors defined so far are plotted on the sphere here
import bsphere
b = bsphere.Bloch()
origin = vec3.vec3()

#b.add_arrow( origin, Q/abs(kin) , 'blue')
#b.add_arrow( -kin/abs(kin), origin, 'red')
#b.add_arrow( origin, kA2/abs(kA2), 'red')
#b.add_arrow( origin, kA1/abs(kA1), 'green')

b.add_arrow( origin, Q100/abs(kin100) , 'blue')
b.add_arrow( -kin100/abs(kin100), origin, 'red')
b.add_arrow( origin, kA2/abs(kA2), 'orange')
b.add_arrow( origin, kA1/abs(kA1), 'green')
b.add_arrow( origin, kMANTA/abs(kMANTA), 'red')

#b.show()

if __name__ == "__main__":
    verbose = True
else:
    verbose = False

##### VERTICAL LATTICE BEAM TILTED BY 30 mrad +Y, 20 mrad -X #####

# Direct lattice (For AFM lattice spacing is doubled ) 
a1 = 2. * l1064/2 * vec3.vec3( 1., 0., 0.)
a2 = 2. * l1064/2 * vec3.vec3( 0., 1., 0.)

a3 = vec3.vec3()
# Deviations observed in mirror mount
dy = -1.  # inch
dx = -0.5 
# inch  # Be careful with the sign of the arctangent
L = 88.0 / 2.54
dphi = np.arctan( dy/ dx) 
dtheta = np.sqrt(dy**2 + dx**2 ) / L

a3.set_spherical(2. * l1064/2 ,  dtheta, dphi)

# Reciprocal lattice
b1 = 2 * np.pi * vec3.cross( a2, a3) /  ( a1 * vec3.cross(a2, a3) ) 
b2 = 2 * np.pi * vec3.cross( a3, a1) /  ( a1 * vec3.cross(a2, a3) ) 
b3 = 2 * np.pi * vec3.cross( a1, a2) /  ( a1 * vec3.cross(a2, a3) )

Qtilt =  -b1 - b2 + b3

##print a3/abs(a3)
##print (a3/abs(a3)).get_spherical() 

btilt = bsphere.Bloch()
###btilt.add_arrow( origin, a1/abs(a1) , 'blue')
###btilt.add_arrow( origin, a2/abs(a1) , 'blue')
btilt.add_arrow( origin, a3/abs(a1) , 'blue')
###btilt.add_arrow( origin, b1/abs(b1) , 'red')
###btilt.add_arrow( origin, b2/abs(b1) , 'red')
###btilt.add_arrow( origin, b3/abs(b1) , 'red')
###btilt.add_arrow( origin, Qtilt/ abs(kin), 'green')
###btilt.add_arrow( origin, Q/ abs(kin), 'black')
###btilt.show()


if verbose:
    print
    print "### TILTED TOP LATTICE BEAM ###\n"
    print "Qtilt spherical coords.:"
    print (Qtilt/abs(Qtilt)).get_spherical()
    print "Q spherical coords.:"
    print (Q/abs(Q)).get_spherical()
    # Calculate angle for HHH Bragg conditiion
    # with respect to Q vector
    print "Percent difference between Q and Qtilt:"
    print 100*(abs(Qtilt)-abs(Q))/abs(Q)
    print 100*(Q-Qtilt)/abs(Q)

braggTHtilt = np.arccos( abs(Qtilt) / 2. / (2.*np.pi/l671)  ) 
if verbose:
    print 
    print "HHH Bragg angle wrt Qtilt = ", braggTHtilt * 180. / np.pi
    print "Delta Bragg angle (tilt/notilt) = ",
    print (braggTH - braggTHtilt) * 1000. , "mrad"

# Find the actual kinTilt that satisfies exactly the Bragg condition
# First find two vectors that are perpendicular to Qtilt
Qtilt_p1 = vec3.cross( Qtilt, kA2 ) 
Qtilt_p1 = Qtilt_p1 / abs(Qtilt_p1) * abs(Qtilt/2) * np.tan(braggTHtilt)
Qtilt_p2 = vec3.cross( Qtilt_p1, Qtilt )
Qtilt_p2 = Qtilt_p2 / abs(Qtilt_p2) * abs(Qtilt/2) * np.tan(braggTHtilt)

# This plots them on the sphere for checking
###tharray = np.linspace( 0., 2*np.pi, 30 )
###kTilt_p = [ Qtilt/2. + Qtilt_p1*np.sin(th) + Qtilt_p2*np.cos(th) for th in tharray  ] 
###for kt in kTilt_p:
###  btilt.add_arrow( origin, kt/abs(kin), 'purple' )
###btilt.show()

# I want to find the Bragg output vector that is closest to kA2
def delta_kA2 (  theta ):
   kt = Qtilt/2. + Qtilt_p1*np.sin(theta) + Qtilt_p2*np.cos(theta)
   return np.arccos( kt * kA2 / ( abs(kt) * abs(kA2) )  )

# Here it can be verified graphically that the minimum is indeed at theta=0
###thX = np.linspace( -np.pi/16, np.pi/16, 100)
###thY = np.array([ delta_kA2(th)  for th in thX])
###import matplotlib.pyplot as plt
###plt.plot( thX, thY)
###plt.show()

# The same theta=0 minimum is obtained using a numerical minimization
###th_min = scipy.optimize.brent( delta_kA2)
###print th_min

kOutTilt = Qtilt/2. + Qtilt_p2 
kinTilt = kOutTilt - Qtilt 



if verbose:
    print
    print "Angle between current output and Qtilt = ",
    thetaA2tilt = np.arccos( Qtilt * kA2 / ( abs(Qtilt) * abs(kA2) ) ) 
    print thetaA2tilt * 180./np.pi, "deg"
    print "Deviation of current kA2 from Bragg condition =",
    print (braggTHtilt - thetaA2tilt ) * 1000. , "mrad"
    
    #kinTilt = kA2 - Qtilt
    print 
    print "Angle between current input and kinTilt =",
    print np.arccos( kinTilt * kin / ( abs(kinTilt) * abs(kin) ) ) *1000., "mrad"
    
    kinS = kin.get_spherical()
    kinTiltS = kinTilt.get_spherical()
    print " dTheta = ",  (kinS[1] - kinTiltS[1])*1000.
    print " dPhi   = ",  (kinS[2] - kinTiltS[2])*1000.


# Here I printed out a short description of the system to send
# to the theorists
def printsph( l, k,  ):
    sph   = (k/abs(kin)).get_spherical()
    cartU = (k/abs(kin)) 
    cart  = (k/abs(kin)) * 532./671. 

    cartA1 = (kA1/abs(kin)) * 532./671.
    cartA2 = (kA2/abs(kin)) * 532./671.
    cartM = (kMANTA/abs(kin)) * 532./671.
  
    QA1 = cartA1 - cart
    QA2 = cartA2 - cart
    QM = cartM - cart

    cstr = '(%+.3f, %+.3f, %+.3f)'

    print ('%16s = (%+.3f*pi, %+.3f*pi) =  '+cstr+' =  '+cstr+' ==>  '+cstr +'  '+cstr+'  '+cstr) % \
          (l,sph[1]/np.pi,sph[2]/np.pi, cartU[0],cartU[1],cartU[2], cart[0], cart[1], cart[2], \
          QA1[0],QA1[1],QA1[2], QA2[0],QA2[1],QA2[2], QM[0],QM[1],QM[2]) 


if verbose:
    print
    print "##### SYSTEM DESCRIPTION #####\n"
    print "Optical lattice original design has three input beams which propagate in directions:\n"
    print "1. +x  (0.500*pi, 0.000*pi)"
    print "2. -y  (0.500*pi, -0.500*pi)" 
    print "3. -z  (1.000*pi, 0.000*pi)"
    print
    print "These three beams are retro reflected to form the lattice.\n" 
    
    
    print "For the beams on the xy plane we are confident that they point along the intended direction, however the beam along z is tilted."  
    print "As a result, in our actual setup the input beams propagate in the following directions:\n" 
    
    printsph( '1.       +x', a1)
    printsph( '2.       -y', -a2)
    printsph( '3. tilted z', -a3)
    
    print
    print "List of available input k vectors."
    print "The pair represents polar and azimuthal angle."
    print "Example: the HHH Input light propagates along +y\n" 
    
    print "\t\t   Spherical \t\t     Unit Cartesian \t\t Normed : |k671|==532/671 \tk_A1 - k_Input \t\t  k_A2 - k_Input \t    k_M - k_Input"
    #printsph('HHH Input', kin)
    printsph('100 Input', kin100)
    for i, k in enumerate(ksquad):
        printsph('Input #%d'%i, k)
    
    print 
    print "List of available output k vectors."
    print "The pair represents polar and azimuthal angle."
    print "Example: The ANDOR1 camera is on the xy plane, at the line y=x*tan(60deg) "
    printsph('ANDOR1', kA1)
    printsph('ANDOR2 (HHH)', kA2)
    printsph('MANTA  (100)', kMANTA)
     
    ###print "kin100",kin100/abs(kin100)
    ###print (kin100/abs(kin100)).get_spherical()
    ###print "kA1",kA1/abs(kA1)
    ###print (kA1/abs(kA1)).get_spherical()
    ###print "kA2",kA2/abs(kA2)
    ###print (kA2/abs(kA2)).get_spherical()
    ###print "kMANTA",kMANTA/abs(kMANTA)
    ###print (kMANTA/abs(kMANTA)).get_spherical()
