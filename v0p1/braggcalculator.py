
import numpy as np
import scipy.constants as spC
import scipy as sp 
from scipy.special import sph_harm

from scipy import weave

def rotation_matrix_numpy(axis, theta):
    mat = np.eye(3,3)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.)
    b, c, d = -axis*np.sin(theta/2.)

    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])




class vector():
    def __init__(self, x, y , z):
        self.x = x
        self.y = y
        self.z = z

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

        th = np.linspace(50.*np.pi/180., 60.*np.pi/180., 20) 
        ph = np.linspace(0., 2*np.pi, 80) 
        TH, PH = np.meshgrid( th, ph)

        THravel = TH.ravel()
        PHravel = PH.ravel()
    
        sums=[]

        rotationaxis = np.cross( [0,0,1], [Q.x, Q.y, Q.z] ) 
        rotationangle = np.dot( [0,0,1], [Q.x, Q.y, Q.z]) / np.linalg.norm( [Q.x, Q.y, Q.z] )

        print rotationaxis
        print rotationangle

        mat1 = rotation_matrix_numpy( rotationaxis, rotationangle)
 
        for i in range( THravel.size ):
            kIN = kvector( THravel[i], PHravel[i] )
            #Rotate kIN so that TH and PH are with respect to the
            #Bragg normal vector 
            krot = np.dot( mat1, [kIN.x, kIN.y, kIN.z] )
            kIN = kvector( krot[0], krot[1], krot[2] )
             
            kOUT = kvector( kIN.x + Q.x , kIN.y + Q.y, kIN.z + Q.z ) 
            sums.append( self.S( kIN, kOUT) )

        return TH, PH, np.array( sums).reshape( TH.shape) 
            

    
 


sample = crystal( 'afm', 20. , 0.532)
Q = vector( -np.pi/0.532, np.pi/0.532, np.pi/0.532)

#Find output vector for the Bragg condition that we are using



th, ph, S = sample.getS( Q )

from mayavi import mlab


#mlab.clf()
#mlab.mesh(kx, ky, kz, scalars=s, colormap='jet')
r=1.0
kx = -r*np.sin(th)*np.cos(ph)
ky = -r*np.sin(th)*np.sin(ph)
kz = -r*np.cos(th)

mlab.clf()
mlab.mesh(kx, ky, kz, scalars=S, colormap='jet')
mlab.colorbar(title='BraggIntensity', orientation='vertical')
mlab.savefig('bragg.png')
#mlab.contour3d(X,Y,Z,PH)
