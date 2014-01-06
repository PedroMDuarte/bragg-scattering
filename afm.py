
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

import os
path = os.path.dirname ( os.path.abspath(__file__) )
CSdictfile = path + '/CSafmdict.pck'

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

        if kin_out == None:
            self.set_kvectors( braggvectors.kin, braggvectors.kout ,\
                               braggvectors.kipol)
        else:
            self.set_kvectors( kin_out[0], kin_out[1],\
                               braggvectors.kipol)
        self.set_pbragg( 250. )
        self.set_detuning( 0. )
        self.lBragg = braggvectors.l671

        self.set_v0( [20.,20.,20.] ) 
        self.set_timeofflight( 0. ) 
        self.pol()
        self.debyewaller()

        self.SumsDone = False
        self.SpinsInit = False
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

        # Helper objects for faster spin-shuffling:
        # Mask for isolating the AFM domain 
        afm0 = (self.latsize-self.afmsize) /2 
        afm1 = afm0 + self.afmsize
        self.mask =   (self.x >= afm0) & (self.x < afm1) \
               & (self.y >= afm0) & (self.y < afm1) \
               & (self.z >= afm0) & (self.z < afm1)
        self.shell = np.where(~self.mask)  
 
        # A completely afm ordered lattice
        self.spin0 = (self.x + self.y + self.z)%2 - 0.5
        maspin = np.ma.array( self.spin0 , mask = self.mask)
        # Get the spins on the metallic core 
        self.metal = maspin.compressed()
        self.metal.reshape( self.metal.size )

    def shuffle_spins(self):
        '''Calling this method creates a random distribution of the spins 
           in the lattice.  The spin values are +/- 0.5
            ''' 
        # A completely afm ordered lattice
        self.spin = (self.x + self.y + self.z)%2 - 0.5
        # Randomize the spins on the metalic shell
        np.random.shuffle(self.metal) 
        # Then assign the randomized spins to the metallic
        # shell of the crystal 
        self.spin[ self.shell ] = self.metal

    def init_spins(self, Nr):
        '''This function precalculates random distributions of spins 
           for use later on by the functions that calculate the sums'''
        self.RandomSpins = []
        for i in range(Nr):
            self.shuffle_spins()
            self.RandomSpins.append( self.spin )
        self.SpinsInit = True
   
    def test_spin_distribution(self):
        '''This method verfies that the AFM core is indeed AFM
           ordered''' 
        sdiff = self.spin - self.spin0 
        core  = np.ma.array( sdiff, mask = ~self.mask ) 
        shell = np.ma.array( sdiff, mask = self.mask )
 
        #print "CORE:"
        #print core
        #print 
        #print "SHELL:"
        #print shell

        print "CORE  deviation from AFM = %.2f" % np.sum(core.compressed()**2)
        print "SHELL deviation from AFM = %.2f" % np.sum(shell.compressed()**2)
        

    #####################################################
    #  METHODS USED TO VISUALIZE THE SPIN DISTRIBUTION         
    #####################################################
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
                           color = 'red', marker='o', s=self.point_size)
    def show(self):
        self.make_plot()
        if self.fig:
            show(self.fig)

    
    #####################################################
    #  METHODS USED TO SET AND CALCULATE GENERAL SCATTERING PROPERTIES         
    #####################################################
    def set_kvectors( self, kin, kout, kipol):
        self.kin = kin
        self.kout = kout
        self.kipol = kipol
        self.pol()
        self.Crandom = None
        self.Srandom = None

    def set_pbragg( self, pbragg):
        self.pbragg = pbragg
        w0 = 0.05 #cm 
        Isat = 5.21 # mW/cm^2
        self.isat = 2.*(pbragg/1000.) / np.pi / (w0**2) / Isat

    def set_detuning( self, det):
        self.det = det
        self.d12 = 76./5.9 
        # det is with respect to in between 1 and 2
        # The units are linewidths 
        self.det2 = self.det + self.d12/2.
        self.det1 = self.det - self. d12/2.
        self.SumsDone = False

    def set_v0( self, v0):
        self.v0 = v0

    def set_timeofflight( self, t):
        self.timeofflight = t

    def pol(self):
        '''This method calculates the sum over final polarizations'''
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
            term =  abs( (splus * pout)*(inpol*sminus) )**2.
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

        try:
          t = self.timeofflight
        except Exception as e:
          print "timeofflight has not been defined in the crystal"
          print "program will stop"
          print e
          exit()
 
        K = kin - kout
       
        dwx = np.exp( - (K.x * self.a / np.pi)**2. / 2. / np.sqrt(v0[0])) 
        dwy = np.exp( - (K.y * self.a / np.pi)**2. / 2. / np.sqrt(v0[1]))
        dwz = np.exp( - (K.z * self.a / np.pi)**2. / 2. / np.sqrt(v0[2]))

        # This is the Debye-Waller factor at TOF = 0 
        dw0 = dwx * dwy * dwz

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

        # The total Debye-Waller factor including 
        # initial value and time of flight deca
        DW = dw0 * dw 
        if verbose:
            print "\tDebye-Waller factor for v0=%.2 Er, TOF=%.2 us"%(v0,t), DW
        self.DW = DW
        return self.DW

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

    #####################################################
    #  CALCULATION OF THE FINAL EXPRESSION FOR THE INTENSITY 
    #####################################################
    def Sums_random(self, Nr):
        Q = self.kout-self.kin  
        phase = np.exp( 1j * Q.x * self.x*self.a + \
                        1j * Q.y * self.y*self.a + \
                        1j * Q.z * self.z*self.a )
   
        self.Nr = Nr 
        self.Sums = []
        for i in range(Nr):
            delta = self.det - self.RandomSpins[i] * self.d12 
            satparam = 2 * self.isat / (  1 + 4*np.power( delta, 2 ) ) 
        
            # Upsilon - 2i Phi 
            UpsPhi = np.absolute( \
                          np.sum( satparam / ( 1 + satparam ) \
                                  * phase * (1 - 2*1j*delta) ) )**2
        
            # Kappa
            kappa = np.sum( satparam / (1  + satparam ) )
            # Xi  
            xi = np.sum( satparam / (1 + satparam)**2 )
            self.Sums.append( (UpsPhi, kappa, xi) )
        self.SumsDone = True
 
    def Intensity_random(self, Nr):
        iarray=[]
        if not self.SpinsInit:
            self.init_spins(Nr)
        if not self.SumsDone:
            self.Sums_random(Nr)  
      
        for i in range(Nr):
            UpsPhi = self.Sums[i][0] 
            kappa  = self.Sums[i][1]
            xi     = self.Sums[i][2]
            iarray.append( self.polsum * ( kappa + \
                             self.DW * ( UpsPhi/(2*self.isat) - xi ))  ) 
        iarray=np.array( iarray ) 
        return ufloat( np.mean( iarray), stats.sem( iarray) )

    ############# FUNCTIONS TO FACILITATE THE CREATION OF PLOTS

    def dw_( self, v0, tof):
        '''This function is used to plot the lattice depth dependence 
           of the Debye-Waller factor''' 
        self.set_v0([v0,v0,v0])
        self.set_timeofflight( tof ) 
        return self.debyewaller()

    def I_tof( self, Nr, tof):
        self.set_timeofflight( tof )
        self.debyewaller() 
        return  self.Intensity_random( Nr )

    def I_( self, Nr=10, detuning=0., v0=20., tof=0., pbragg=250.):
        self.set_detuning( detuning )
        self.set_v0([v0,v0,v0])
        self.set_timeofflight( tof )
        self.set_pbragg( pbragg )
        self.debyewaller()
        self.SumsDone = False 
        return  self.Intensity_random( Nr )
        
        
def debye_waller_Q( innum, cam, v0, TOF ):
    N = 40. 
    nafm = 8.
    kinput = braggvectors.ksquad[ innum ]  
    if cam == 'A1': 
        koutput = braggvectors.kA1
    elif cam == 'A2':
        koutput = braggvectors.kA2 
    else:
         raise Exception("Camera is not defined")
    crys = crystal( N, nafm, braggvectors.l1064/2, (kinput, koutput) )
    Q = (koutput - kinput) / abs(koutput) *532. / 671. 
    #print 'in=%d, cam=%s'%(innum,cam), 'Q = ', Q
    return crys.dw_(v0, TOF), Q

def debye_waller_Q_kin( kinput, cam, v0, TOF ):
    N = 40. 
    nafm = 8.
    if cam == 'A1': 
        koutput = braggvectors.kA1
    elif cam == 'A2':
        koutput = braggvectors.kA2 
    else:
         raise Exception("Camera is not defined")
    crys = crystal( N, nafm, braggvectors.l1064/2, (kinput, koutput) )
    Q = (koutput - kinput) / abs(koutput) *532. / 671. 
    #print 'in=%d, cam=%s'%(innum,cam), 'Q = ', Q
    return crys.dw_(v0, TOF), Q
 
 
