
import numpy as np
import scipy.constants as spC
import scipy as sp 
from scipy.special import sph_harm

class kvector():
    """Contains information about the light in bragg scattering."""
    def __init__(self, *args):
        if len(args) == 2:  self.__initSPH(*args)
        if len(args) == 3:  self.__initXYZ(*args)
    def __initSPH(self, theta, phi):
        # The k vector has a magnitude of k = 2pi/lambda
        self.x = np.sin(theta)*np.cos(phi) * 2 *np.pi / 0.671
        self.y = np.sin(theta)*np.sin(phi) * 2 *np.pi / 0.671
        self.z = np.cos(theta) * 2*np.pi / 0.671 
    def __initXYZ(self, kx, ky, kz):
        norm = np.linalg.norm( np.array([kx,ky,kz]) )
        self.x =  kx / norm * 2 *np.pi / 0.671
        self.y =  ky / norm * 2 *np.pi / 0.671
        self.z =  kz / norm * 2 *np.pi / 0.671
        
    

class crystal():
    """Contains information about crystal. 
       Used to calculate scattering of light"""
    def __init__(self, ordering, size, spacing):
        self.ordering = ordering
        self.size = size # number of lattice sites on the side
        self.a = spacing # lattice spacing
        #self.x, self.y, self.z = np.ogrid[ 0:size, 0:size, 0:size ]
        self.x, self.y, self.z = np.mgrid[ -size/2:size/2:size*1j, \
                                           -size/2:size/2:size*1j, \
                                            -size/2:size/2:size*1j]


    def S( self, ki, kf ):
        phase = np.exp( 1j * (ki.x-kf.x)*self.x*self.a + \
                        1j * (ki.y-kf.y)*self.y*self.a + \
                        1j * (ki.z-kf.z)*self.z*self.a )
        spin = (self.x + self.y + self.z)%2 - 0.5
 
        return np.abs( np.sum( phase*spin ) )**2

    def getS( self, Q):

        th = np.linspace(0., np.pi, 40) 
        ph = np.linspace(0., 2*np.pi, 40) 
        TH, PH = np.meshgrid( th, ph)

        THravel = TH.ravel()
        PHravel = PH.ravel()
    
        sums=[]
        for i in range( THravel.size ):
            kIN = kvector( THravel[i], PHravel[i] ) 
            kOUT = kvector( kIN.x + Q.x , kIN.y + Q.y, kIN.z + Q.z ) 
            sums.append( self.S( kIN, kOUT) )

        return TH, PH, np.array( sums).reshape( TH.shape) 
            

    def get_phase_diff(self, thetaIN, phiIN, thetaOUT, phiOUT):

        kxIN = np.sin(thetaIN)*np.cos(phiIN) * 2 *np.pi / self.wl
        kyIN = np.sin(thetaIN)*np.sin(phiIN) * 2 *np.pi / self.wl
        kzIN = np.cos(thetaIN) * 2*np.pi / self.wl         
        #phaseIN = self.a*( np.outer(kxIN,self.x) + \
        #                   np.outer(kyIN,self.y) + \
        #                   np.outer(kzIN,self.z) ) 
        phaseIN = self.a*( kxIN*self.x + \
                           kyIN*self.y + \
                           kzIN*self.z ) 

        kxOUT = np.sin(thetaOUT)*np.cos(phiOUT) * 2 *np.pi / self.wl
        kyOUT = np.sin(thetaOUT)*np.sin(phiOUT) * 2 *np.pi / self.wl
        kzOUT = np.cos(thetaOUT) * 2*np.pi / self.wl 
        #phaseOUT = self.a*( np.outer(kxOUT,self.x) + \
        #                    np.outer(kyOUT,self.y) + \
        #                    np.outer(kzOUT,self.z) ) 
        phaseOUT = self.a*( kxOUT*self.x + \
                            kyOUT*self.y + \
                            kzOUT*self.z ) 
 
        return phaseIN - phaseOUT
    
    def get_int(self, thetaIN, phiIN):
        kxIN = np.sin(thetaIN)*np.cos(phiIN) * 2 *np.pi / self.wl
        kyIN = np.sin(thetaIN)*np.sin(phiIN) * 2 *np.pi / self.wl
        kzIN = np.cos(thetaIN) * 2*np.pi / self.wl 

        Qx = 2* np.pi / self.a
        Qy = 0.
        Qz = 0.

        kxOUT = kxIN + Qx
        kyOUT = kyIN + Qy
        kzOUT = kzIN + Qz

        return  np.abs( np.sum( np.exp(1j*\
                       self.a*( (kxIN-kxOUT)*self.x + \
                                (kyIN-kyOUT)*self.y + \
                                (kzIN-kzOUT)*self.z ) \
                               ) \
                               ))**2

    def get_int_vec( self):

        th = np.linspace(0., np.pi, 40) 
        ph = np.linspace(0., 2*np.pi, 40) 
        TH, PH = np.meshgrid( th, ph)

        print TH.shape

        THravel = TH.ravel()
        PHravel = PH.ravel()
    
        sums=[]
        for i in range( THravel.size ):
            sums.append( self.get_int( THravel[i], PHravel[i]) )

        return TH, PH, np.array( sums).reshape( TH.shape) 
            

       

    
 


sample = crystal( 'afm', 20. , 0.532)
kIN = kvector( np.pi/2 + 1./8., 0. ) 
kOUT = kvector( 

braggIN  = light( np.pi/2 + 0.75/8., 0. )
from mayavi import mlab
X, Y, Z = sample.x, sample.y, sample.z



#PH = sample.get_intensity_out( np.pi/2 + 1./8., 0.)
PH = sample.get_phase_diff( np.pi/2 - 1./8., 0.,  \
                            np.pi/4, np.pi/4)
print PH.shape

r = 0.3
pi = np.pi
cos = np.cos
sin = np.sin
th, ph = np.mgrid[0:pi:21j, 0:2*pi:21j]

kx = r*np.sin(th)*np.cos(ph)
ky = r*np.sin(th)*np.sin(ph)
kz = r*np.cos(th)

s = sph_harm(1, 1, ph, th).real
print s.shape
print "Hello"
#get_int_vec = np.vectorize( sample.get_int) 
#INT = get_int_vec( kx, ky, kz, np.pi/4,  np.pi/4 )
INT = sample.get_int( np.pi/2-1./8., 0. )

print s.shape
print INT

#mlab.clf()
#mlab.mesh(kx, ky, kz, scalars=s, colormap='jet')

th, ph, INTVEC =  sample.get_int_vec( )
kx = r*np.sin(th)*np.cos(ph)
ky = r*np.sin(th)*np.sin(ph)
kz = r*np.cos(th)
print INTVEC
print kx.shape
print INTVEC.shape

mlab.clf()
mlab.mesh(kx, ky, kz, scalars=INTVEC, colormap='jet')
mlab.colorbar(title='BraggIntensity', orientation='vertical')

#mlab.contour3d(X,Y,Z,PH)

