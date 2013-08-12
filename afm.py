
import numpy as np

import scipy.constants as spC
import scipy as sp 
from scipy.special import sph_harm

import vec3
import bsphere

import braggvectors

from pylab import figure, plot, show, savefig, close
from mpl_toolkits.mplot3d import Axes3D

# Pickle will be used to store and retrieve 
# previously calculated C and S sums
import cPickle as pickle
CSdictfile = 'CSafmdict.pck'
import hashlib
# A hash is created to uniquely identify a set of kvectors
# and crystal parameters that define C,S sums
def khash( N, nafm, Nr, kin, kout, kipol):
    kh = hashlib.sha224( str(kin) + str(kout) + str(kipol) ).hexdigest() 
    tag = '%d_%d_%d_%s'%(N, nafm, Nr, kh)
    return tag

from scipy import stats
from uncertainties import unumpy, ufloat



verbose = False

class crystal():
    """Contains information about crystal. 
       Used to calculate scattering of light"""
    def __init__(self, latsize, afmsize, spacing, kin_out=None):
        # crystal settings
        self.latsize = latsize # number of lattice sites on the side
        self.afmsize = afmsize # number of lattice sites with afm
        self.a = spacing # lattice spacing
        self.mass = 6.
        #self.x, self.y, self.z = np.ogrid[ 0:size, 0:size, 0:size ]
        self.x, self.y, self.z = np.mgrid[ 0:self.latsize, \
                                           0:self.latsize, \
                                           0:self.latsize ]

        # default values for scattering parameters
        self.set_detuning( 0. )
        self.set_pbragg( 250. )
        self.lBragg = braggvectors.l671
        self.sunits( ) 

        if kin_out == None:
            self.set_kvectors( braggvectors.kin, braggvectors.kout ,\
                               braggvectors.kipol)
        else:
            self.set_kvectors( kin_out[0], kin_out[1],\
                               braggvectors.kipol)
        
        self.set_v0( [20.,20.,20.] ) 
        
        

        # crystal visualization settings
        self.fig = None
        self.axes = None
        self.user_fig = None
        self.user_axes = None
        # The size of the figure in inches, default = [7,7].
        self.size = [7, 7]
        # Azimuthal and Elvation viewing angles, default = [-60,30].
        self.view = [-60, 30]
        # Labels for x-axis (in LaTex), default = ['$x$','']
        self.xlabel = ['$x$', '']
        # Position of x-axis labels, default = [1.2,-1.2]
        self.xlpos = [1.2, -1.2]
        # Labels for y-axis (in LaTex), default = ['$y$','']
        self.ylabel = ['$y$', '']
        # Position of y-axis labels, default = [1.1,-1.1]
        self.ylpos = [1.1, -1.1]
        # Labels for z-axis (in LaTex),
        # default = ['$\left|0\\right>$','$\left|1\\right>$']
        self.zlabel = ['$z$', '']
        # Position of z-axis labels, default = [1.2,-1.2]
        self.zlpos = [1.2, -1.2]
        #---font options---
        # Color of fonts, default = 'black'
        self.font_color = 'black'
        # Size of fonts, default = 20
        self.font_size = 20
        #---scatter options---
        self.point_size = 216
   
 
    def sunits(self ): 
        # lBragg is the wavelength of the Bragg probe in nm
        self.sunits = 9. * (self.lBragg**2.) / 16 / (np.pi**2) 
        

    def set_pbragg( self, pbragg):
        self.pbragg = 250.
        w0 = 0.05 #cm 
        Isat = 5.21 # mW/cm^2
        self.isat = 2.*(pbragg/1000.) / np.pi / (w0**2) / Isat

    def set_kvectors( self, kin, kout, kipol):
        self.kin = kin
        self.kout = kout
        self.kipol = kipol
        self.pol()
        self.Crandom = None
        self.Srandom = None

    def set_v0( self, v0):
        self.v0 = v0

    def set_timeofflight( self, t):
        self.timeofflight = t

    def set_detuning( self, det):
        self.det = det
        self.d12 = 76./5.9 
        # det is with respect to in between 1 and 2 
        self.det2 = self.det + self.d12/2.
        self.det1 = self.det - self. d12/2.


    def shuffle_spins(self):
        # First make a completely afm ordered lattice
        self.spin = (self.x + self.y + self.z)%2 - 0.5
        # Then isolate the central afm core by using a mask
        afm0 = (self.latsize-self.afmsize) /2 
        afm1 = afm0 + self.afmsize
        mask =   (self.x >= afm0) & (self.x < afm1) \
               & (self.y >= afm0) & (self.y < afm1) \
               & (self.z >= afm0) & (self.z < afm1)  
        maspin = np.ma.array( self.spin , mask = mask)
        # Get the spins on the metallic core and 
        # randomize their arrangement
        metal = maspin.compressed()
        metal.reshape( metal.size )
        np.random.shuffle(metal) 
     
        # Then assign the randomized spins to the metallic
        # shell of the crystal 
        shell =  np.where(~mask) 
        for i in range(shell[0].size):
            maspin[ shell[0][i], shell[1][i], shell[2][i] ] = metal[i]
            self.spin[ shell[0][i], shell[1][i], shell[2][i] ] = metal[i]

    ###############
    #  METHODS USED TO VISUALIZE THE SPIN DISTRIBUTION         
    def make_plot(self):
        # setup plot
        # Figure instance for Bloch sphere plot
        if self.user_axes:
            self.axes = self.user_axes
        else:
            if self.user_fig:
                self.fig = self.user_fig
            else:
                self.fig = figure(figsize=self.size)
            self.axes = Axes3D(self.fig, azim=self.view[0],
                               elev=self.view[1])
        self.axes.clear()
        #self.axes.grid(False)

        spinup = np.where( self.spin == 0.5 )
        spindn = np.where( self.spin == -0.5 )

        self.axes.scatter( spinup[0], spinup[1], spinup[2], 
                           color = 'blue', marker='o', s=self.point_size)

        self.axes.scatter( spindn[0], spindn[1], spindn[2],
                           color = 'blue', marker='o', s=self.point_size)
    def show(self):
        self.make_plot()
        if self.fig:
            show(self.fig)
    ###############

    def pol(self):
        try:
          kin = self.kin
          kout = self.kout
          kipol = self.kipol
        except Exception as e:
          print "k vectors have not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        # unit vectors for input polarization
        zvec = vec3.vec3(0.,0.,1.) 
        in1 = vec3.cross( zvec, kin )
        in1 = in1 / abs(in1) 
        in2 = vec3.cross( kin, in1)
        in2 = in2 / abs(in2)
        # input polarization
        inpol = kipol[0]*in1 + kipol[1]*in2
        
        # unit vectors for output polarization 
        out1 = vec3.cross( zvec, kout )
        out1 = out1 / abs(out1) 
        out2 = vec3.cross( kout, out1)
        out2 = out2 / abs(out2)

        # polarization of transition
        splus =  1/np.sqrt(2.) *( vec3.vec3(1.,0.,0.) + 1j * vec3.vec3(0.,1.,0.))
        sminus = 1/np.sqrt(2.) *( vec3.vec3(1.,0.,0.) - 1j * vec3.vec3(0.,1.,0.))

        # sum over output polarizations
        
        polsum = 0.
        for i,pout in enumerate([out1, out2]):
            term =  abs((pout*splus) * (sminus*inpol))**2.
            polsum = polsum + term
        if verbose:
            print "\tSum over output pol (vectors) = ", polsum
        self.polsum = polsum
   
        # also obtain this result using the angles 
        # with respect to the magnetic field
        cosIN  = kin/abs(kin) * zvec 
        cosOUT = kout/abs(kout) * zvec
        polsumANGLE =  1./4. * (1 + cosIN**2.) * (1 + cosOUT**2.)
        #print "Sum over output pol (angles)  = ", polsumANGLE
   
        # optional draw polarization vectors on bloch sphere
        if False: 
            b = bsphere.Bloch()
            origin = vec3.vec3()
            b.add_arrow( -kin/abs(kin), origin, 'red')
            b.add_arrow( -kin/abs(kin), -kin/abs(kin) + in1/2., 'black')
            b.add_arrow( -kin/abs(kin), -kin/abs(kin) + in2/2., 'black')

            b.add_arrow( origin, kout/abs(kout), 'red')
            b.add_arrow( origin, out1/2., 'black')
            b.add_arrow( origin, out2/2., 'black')
        
            b.show()
        return self.polsum

    def debyewaller(self):
        try:
          kin = self.kin
          kout = self.kout
          kipol = self.kipol
        except Exception as e:
          print "k vectors have not been defined in the crystal"
          print "program will stop"
          print e
          exit()
       
        try:
          v0 = self.v0
        except Exception as e:
          print "v0 has not been defined in the crystal"
          print "program will stop"
          print e
          exit()
 
        K = kin - kout
       
        dwx = np.exp( - (K.x * self.a / np.pi)**2. / 2. / np.sqrt(v0[0])) 
        dwy = np.exp( - (K.y * self.a / np.pi)**2. / 2. / np.sqrt(v0[1]))
        dwz = np.exp( - (K.z * self.a / np.pi)**2. / 2. / np.sqrt(v0[2]))

        dw = dwx * dwy * dwz
        if verbose:
            print "\tDebye-Waller factor = ", dw
        self.dw = dw
        return self.dw

    def debyewaller_time( self):
        try:
          kin = self.kin
          kout = self.kout
          kipol = self.kipol
        except Exception as e:
          print "k vectors have not been defined in the crystal"
          print "program will stop"
          print e
          exit()
       
        try:
          v0 = self.v0
        except Exception as e:
          print "v0 has not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        try:
          t = self.timeofflight
        except Exception as e:
          print "timeofflight has not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        K = kin - kout

        # The debye-waller factor goes down as a function
        # of time due to the initial spread in momentum of 
        # the wave-packet

        # The ground state is assumed to be that of a harmonic
        # oscillator at the lattice frequency, it therefore has 
        # a spread in position equal to 
        dx = self.a / np.pi / np.sqrt(2.) / (v0[0]**(1./4.)) 
        dy = self.a / np.pi / np.sqrt(2.) / (v0[1]**(1./4.)) 
        dz = self.a / np.pi / np.sqrt(2.) / (v0[2]**(1./4.)) 

        # In the end we will need the ratio h/m.  If at this point
        # we choose to express m in units of AMU then can write
        hplanck = 3.99e5 # nm^2 / us  

        # The ground state of the harmonic oscillator is a gaussian
        # wave packet and it is Heisenberg limited so we can obtain
        # the spread in momentum from the uncertainty relation
        dpxH = hplanck / 4. / np.pi / dx
        dpyH = hplanck / 4. / np.pi / dy
        dpzH = hplanck / 4. / np.pi / dz 

        # The spread in the momentum of the harmonic oscillator ground
        # state can also be expresed straight up as a function of the
        # lattice depth 
        dpx = hplanck * (v0[0]**(1./4.)) / self.a / np.sqrt(8.) 
        dpy = hplanck * (v0[1]**(1./4.)) / self.a / np.sqrt(8.) 
        dpz = hplanck * (v0[2]**(1./4.)) / self.a / np.sqrt(8.) 
        
        # Both ways of calculation the momentum spread should agree
        if not np.allclose( np.array([dpx,dpy,dpz]), \
                        np.array([dpxH,dpyH,dpzH])):
            print "Discrepancy in calculation of momentum spread of wave-packet"
            print "using Heisenberg uncertainty or direct expression"     
            print " --- HEISENBERG --- "
            print dpxH 
            print dpyH 
            print dpzH 
            print " --- DIRECT --- "
            print dpx
            print dpy 
            print dpz 

        mass = self.mass
 
        # The units of dp as epxpressed above are nm/us
        # With hplanck as defined above mass is unitless

        #print "Norm(K) = %e" % np.sqrt(K.x**2 + K.y**2 + K.z**2)
        #print "2pi/l671 = %e" % (2.*np.pi / 671.0)
        #print "Kx = %e" % K.x
        #print "dpx = %e" % dpx 
        #print "t = %.2f" % t 
   

        dwx = np.exp( - (K.x * dpx * t / mass )**2 )
        dwy = np.exp( - (K.y * dpy * t / mass )**2 )
        dwz = np.exp( - (K.z * dpz * t / mass )**2 )

        dw = dwx * dwy * dwz
        if verbose:
            print "\tDebye-Waller TOF factor = ", dw
        self.dwtof = dw
        return self.dwtof

    def alpha_beta( self):
        try:
          det = self.det
        except Exception as e:
          print "detuning has not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        f1 = -1. / ( 2 * self.det1 + 1j)
        f2 = -1. / ( 2 * self.det2 + 1j)

        self.alpha = 1/4. * np.abs( f1+f2 )**2.
        self.beta  = np.abs( f1-f2  )**2.
        
        if verbose:    
            print "\n\t--- Alpha and Beta ---"
            print "\tdet = ", det
            print "\talpha = ", self.alpha
            print "\tbeta = ", self.beta

    def alpha_beta_Pbroad( self):
        try:
          det = self.det
        except Exception as e:
          print "detuning has not been defined in the crystal"
          print "program will stop"
          print e
          exit()
        try:
          isat = self.isat 
        except Exception as e:
          print "isat (from Pbragg) has not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        f1 = -1. / ( 2 * self.det1/np.sqrt(1+isat) + 1j)
        f2 = -1. / ( 2 * self.det2/np.sqrt(1+isat) + 1j)

        self.alpha_Pbroad = 1/4. * np.abs( f1+f2 )**2.
        self.beta_Pbroad  = np.abs( f1-f2  )**2.
        
        if verbose:    
            print "\n\t--- Alpha and Beta Power Broadened---"
            print "\tdet = ", det
            print "\talpha = ", self.alpha_Pbroad
            print "\tbeta = ", self.beta_Pbroad
        return self.alpha_Pbroad, self.beta_Pbroad

    def coherent(self):
        '''Returns fraction of scattering that is elastic'''
        try:
          det = self.det
        except Exception as e:
          print "detuning has not been defined in the crystal"
          print "program will stop"
          print e
          exit()
        try:
          isat = self.isat 
        except Exception as e:
          print "isat (from Pbragg) has not been defined in the crystal"
          print "program will stop"
          print e
          exit()

        # Saturation parameter
        SP = 2. * isat * ( 0.5 / (1.+4.*(self.det1**2)) + 0.5 / (1.+4.*(self.det2**2))) 
        self.coh =  1. / (1. + SP)
        return self.coh

    def ODfactor(self):

        #Average distance for photon to travel to outside of sample (nm)
        dz = self.latsize * self.a / 2. 
 
        #On-resonace cross section of single atom (nm^2)
        sigm = 3.*np.pi*(self.lBragg/2./np.pi)**2 

        #Lorentzian lineshape
        sigm = sigm * ( 0.5/ (1+4*self.det1**2) + 0.5/ (1+4*self.det2**2) ) 
 
        #Density (1 per site)
        nc = 1. / (self.a**3) 

        #Optical density 
        OD = nc * sigm * dz
        #print "Optical density at det = %.2f --> OD = %.3f . e^(-OD) = %e" % (det,OD,np.exp(-OD)) 
         
        self.odfactor = np.exp(-1.*OD) 
        return self.odfactor

    ############# CRYSTAL AND MAGNETIC STRUCTURE SUMS

    def S( self, ki, kf ):
        if verbose and False:
            print (ki.x - kf.x)*self.a
            print (ki.y - kf.y)*self.a
            print (ki.z - kf.z)*self.a
        phase = np.exp( 1j * (ki.x-kf.x)*self.x*self.a + \
                        1j * (ki.y-kf.y)*self.y*self.a + \
                        1j * (ki.z-kf.z)*self.z*self.a )
        sfactor = np.abs( np.sum( phase*self.spin ) )**2
        if verbose:
            print "\tS = ", sfactor
        return sfactor

    def C( self, ki, kf ):
        if verbose and False:
            print (ki.x - kf.x)*self.a
            print (ki.y - kf.y)*self.a
            print (ki.z - kf.z)*self.a
        
        phase = np.exp( 1j * (ki.x-kf.x)*self.x*self.a + \
                        1j * (ki.y-kf.y)*self.y*self.a + \
                        1j * (ki.z-kf.z)*self.z*self.a )
        cfactor = np.abs( np.sum( phase ) )**2
        if verbose:
            print "\tC = ", cfactor
        return cfactor


    def CS_random( self, Nr): 
        try:
          kin = self.kin
          kout = self.kout
        except Exception as e:
          print "k vectors have not been defined in the crystal"
          print "program will stop"
          print e
          exit()
 
        Csum =[]
        Ssum =[]
        for i in range(Nr):
            self.shuffle_spins()
            Ssum.append( self.S(kin, kout))
            Csum.append( self.C(kin, kout))
     
        return np.array( Csum ), np.array( Ssum)


        
    def saveCS(self,Nr, C, S):
        tag =  khash( self.latsize, self.afmsize, Nr, self.kin, self.kout, self.kipol) 
        try:
            CSdict = pickle.load( open(CSdictfile,"rb") ) 
        except:
            CSdict = {} 
        CSdict['C'+tag]=C
        CSdict['S'+tag]=S
        pickle.dump( CSdict, open(CSdictfile,"wb") ) 
    
    
    def loadCS(self,Nr):
        tag =  khash( self.latsize, self.afmsize, Nr, self.kin, self.kout, self.kipol) 
        try:
            CSdict = pickle.load( open(CSdictfile,"rb") )
            Cr, Sr =  CSdict['C'+tag], CSdict['S'+tag] 
            #print "\nLoaded C, S values from pickle file:\n\t%s" % tag
            self.Crandom = Cr
            self.Srandom = Sr 
        except:
            print "\nError loading C, S from pickle file:\n\t%s" % tag
            print "C,S will be calculated"
            C, S = self.CS_random( Nr)
            self.saveCS(Nr, C, S)
            self.Crandom = C
            self.Srandom = S

        return self.Crandom,self.Srandom


    ############# FUNCTIONS TO FACILITATE THE CREATION OF PLOTS

    def dw_time( self, time): 
        '''This function is used to plot the time dependence of the
           Debye-Waller factor''' 
        self.set_timeofflight(time)
        return self.debyewaller_time() 

    def sigma_incoh_det( self, Nr, det, time):
        '''This function is used to plot the inelastic cross section 
           dependence on detuning, and time'''

        self.set_detuning(det)
        self.set_timeofflight(time)
        if self.Crandom == None or self.Srandom==None:
            self.loadCS(Nr)
        self.alpha_beta_Pbroad()

        crystalpart = self.alpha_Pbroad * ( self.x.size )

        magneticpart = self.beta_Pbroad * ( self.x.size )

        vals =  self.sunits * self.polsum \
               * ( crystalpart + magneticpart)
        return ufloat( vals, 0. ) 

    def sigma_coh_det( self, Nr, det, time):
        '''This function is used to plot the elastic cross section 
           dependence on detuning, and time'''

        self.set_detuning(det)
        self.set_timeofflight(time)
        if self.Crandom == None or self.Srandom==None:
            self.loadCS(Nr)
        self.alpha_beta_Pbroad()
       
        crystalpart = self.alpha_Pbroad * ( self.x.size + \
                      self.debyewaller_time() * self.Crandom * \
                      self.coherent() * self.ODfactor() )

        magneticpart = self.beta_Pbroad * ( self.x.size + \
                      self.debyewaller_time() * self.Srandom * \
                      self.coherent() * self.ODfactor() )

        vals = self.sunits * self.polsum \
               * ( crystalpart + magneticpart)
        return ufloat( np.mean(vals), stats.sem(vals) ) 

            

    
 
