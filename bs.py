import numpy as np
import vec3
import pylab



# Q vector for bragg scattering
l671 = 671.
l1064 = 1064.
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



import afm
 



# Saturation intensity 
PBragg = 250. # microWatt 
sat = 2.*(PBragg/1000.) / np.pi / (0.05**2)  / 5.21 
print "PBragg = %.2f --> sat = %.2f" % (PBragg, sat)

# Parameters
N = 40 

class Kdet:
    def __init__(crystal, kin, kout):

        self.c = crystal

        v0 = [20.,20.,20.]
        kipol = [1., 0]
        PBragg = 250. 
        self.sat = 2.*(PBragg/1000.) / np.pi / (0.05**2)  / 5.21 
      
        self.sunits = 9. * (671e-7**2) / 16 / (np.pi**2) 

        self.pol = self.c.pol(kin, kout, kipol) 
        self.dw  = self.c.debyewaller( kin, kout, v0)

        Nr = 50 # Number of random crystal realizations 
        self.C, self.S = self.c.CS_random( kin, kA1, Nr)
        
 
    def alpha_beta(det)
        self.alpha, self.beta = self.c.alpha_beta_Pbroad( det, self.sat )
 
    def elastic(det):
        return self.sunits * self.polsum * self.dw * ( self.alpha * self.C + self.beta * self.S )

    def inelastic(det ):
        return self.sunits * self.polsum * self.dw * ( self.alpha * crystal.spin.size \
                                                     + self.beta * 1/4. * crystal.spin.size )

# Detuning calculation for validation of equations
for nafm in [4,6,8,10,12,16,20,24]:
    
    # Initialize crystal 
    c = afm.crystal(N, nafm, l1064/2)
        
    A1 = Kdet(c, kin, kA1)
    A2 = Kdet(c, kin, kA2)  

    for det in np.linspace(-10,10,51):
        e1 = A1.elastic(det) 
        i1 = A1.inelastic(det)
        
        e2 = A2.elastic(det)
        i2 = A2.inelastic(det) 
 




#Nafm = 12
#N = 40
#Nr   = 20
#v0 = [7.,7.,7.]
#d12 = 76./5.9
#aA1 = 1000.*(np.arccos( Qunit*kA1 / abs(kA1) ) - braggTH)
#kipol = [1., 0]
# 
#detsA2=[]
#detsA1=[]
#c = afm.crystal(N, Nafm, l1064/2.)
#for det in np.linspace(-600, 600,7):
#    A2sigma = c.cross_section_random( kin, kout, kipol, v0, det, 0., 0., Nr)
#    detsA2.append([det,A2sigma[0],A2sigma[1]])
#    A1sigma = c.cross_section_random( kin, kA1, kipol, v0, det, 0., aA1, Nr)
#    detsA1.append([det,A1sigma[0],A1sigma[1]])
#    print "det = %.2f --> A2/A1 = %.2f" % (det, A2sigma[0]/A1sigma[0])

#detsA2 = np.array(detsA2) 
#detsA1 = np.array(detsA1) 
#pylab.plot(detsA2[:,0],detsA2[:,1]/detsA1[:,1])
#pylab.show()
#exit()



def calculate( N, Nafm, Nr, mrad, Npts, det, v0, kipol):
    print
    print "...Calculating for %d sites, %d afm, %d spin realizations" % (N, Nafm, Nr)

    c = afm.crystal(N, Nafm, l1064/2.)
    c.set_vectors( kin, kout, Q)
    #c.show()

    # Use for calculation at an array of input/output angles
    kiangles = np.linspace( -1.* mrad, mrad, Npts)
    kfangles = np.linspace( -1.*mrad, mrad, Npts*4)

    # Use for calculation at exact Bragg angles:
    kiangles = np.array([0.])
    kfangles = np.array([0.])

    print "...kin from %.1f to %.1f,  %d pts" % (kiangles.min(), kiangles.max(), kiangles.size)

    print "...detuning = %.2f gamma" % det 
    print "...lattice depth = ",v0
    print "...kin polarization = ",kipol

    #print "...STARTING CALCULATION FOR kIN"
    for ai in kiangles:
	kia =  -Qunit*np.cos(braggTH + ai/1000.) - Qperp2*np.sin(braggTH + ai/1000.)
	kia = abs(kin) * kia
	b.add_points( (-1*kia/abs(kia)).tolist()  )

	# Obtain cross section along A1 
	aA1 = 1000.*(np.arccos( Qunit*kA1 / abs(kA1) ) - braggTH)
	A1sigma = c.cross_section_random( kia, kA1 , kipol, v0, det, ai, aA1, Nr)

        # Obtain cross section along all kfangles
	for af in kfangles:
	    kfa =  Qunit*np.cos(braggTH + af/1000.) - Qperp2*np.sin(braggTH + af/1000.) 
	    kfa = abs(kin) * kfa
	    b.add_points( (kfa/abs(kfa)).tolist()  ) 
		
	    A2sigma = c.cross_section_random( kia, kfa, kipol, v0, det, ai, af, Nr)



    
# Number of sites with AFM ordering at center of sample
Nafm = 36 
# Number of sites of entire sample
N = 40
# Number of spin realizations per k point
Nr   = 10 #14   

mrad = 150. # kin range in mrad
Npts = 40   # number of kin  points

# Spacing between states 1 and 2 in units of the linewidth
d12 = 76./5.9
 
# For the input polarization the first
# component is along x, and the second is along
# x cross kin
kipol = [1., 0]

# Lattice depth 
v0 = [7.,7.,7.]

# Do detuning data:
#for nafm in [4,6,8,10,12,16,20,24]:
#    for det in np.linspace(-800,800,25):
#        calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
for nafm in [4,6,8,10,12,16,20,24]:
    for det in np.linspace(-10,10,51):
        calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
#for nafm in [32,34,38,40]:
#    for det in np.linspace(-800,800,25):
#        calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
for nafm in [32,34,38,40]:
    for det in np.linspace(-10,10,51):
        calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
#for nafm in [40]:
#    for det in np.linspace(-8,8,25):
#        calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
exit()

nafm = 4
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 6
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 8
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 10
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 12
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 16
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 20
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)
nafm = 24
for det in np.linspace(-800, 800,25):
    calculate( N, nafm, Nr, mrad, Npts, det, v0, kipol)






