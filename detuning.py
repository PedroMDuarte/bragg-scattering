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
 




#Pickle implementation to retrieve previously calculated C and S 
import cPickle as pickle


def saveCS(nafm, Nr, kout, C, S):
    if (kout - kA1)**2 < 1e-6: 
        k=1
    elif (kout - kA2)**2 < 1e-6: 
        k=2
    else:
        print "Invalid kout encountered in C,S pickle save"
        return
    try:
        CSdict = pickle.load( open("CSdict.pck","rb") ) 
    except:
        CSdict = {} 
    CSdict['C%d_%d_%d'%(nafm,Nr,k)]=C
    CSdict['S%d_%d_%d'%(nafm,Nr,k)]=S
    pickle.dump( CSdict, open("CSdict.pck","wb") ) 


def loadCS(nafm, Nr, kout):
    if (kout - kA1)**2 < 1e-6: 
        k=1
    elif (kout - kA2)**2 < 1e-6: 
        k=2
    else:
        print "Invalid kout encountered in C,S pickle load"
        return 0.,0.
    try:
        CSdict = pickle.load( open("CSdict.pck","rb") )
        Cr, Sr =  CSdict['C%d_%d_%d'%(nafm,Nr,k)], CSdict['S%d_%d_%d'%(nafm,Nr,k)] 
        print "Loaded C, S values from pickle file.  nafm=%d, Nr=%d, k=%d"%(nafm,Nr,k)
        return Cr, Sr
    except:
        print "Error loading C, S from pickle file.  nafm=%d, Nr=%d, k=%d"%(nafm,Nr,k)
        print "C,S will be calculated"

        return None,None


class Kdet:
    def __init__(self, N, nafm, kin, kout):

        # Parameters
        self.c = afm.crystal(N, nafm, l1064/2)

        self.kin = kin
        self.kout = kout

        self.v0 = [20.,20.,20.]
        self.kipol = [1., 0]
        PBragg = 250. 
        self.sat = 2.*(PBragg/1000.) / np.pi / (0.05**2)  / 5.21 
        #print "PBragg = %.2f --> sat = %.2f" % (PBragg, self.sat)
      
        self.sunits = 9. * (671e-7**2) / 16 / (np.pi**2) 

        self.pol = self.c.pol(self.kin, self.kout, self.kipol) 
        self.dw  = self.c.debyewaller( self.kin, self.kout, self.v0)

        Nr = 320 # Number of random crystal realizations 
        self.C, self.S = loadCS(nafm,Nr,kout) 
        if self.C == None or self.S == None:
           self.C, self.S = self.c.CS_random( self.kin, self.kout, Nr)
           saveCS(nafm,Nr,self.kout, self.C, self.S)
           
   
    def dw_time( self, time):
        return self.c.debyewaller_time( self.kin, self.kout, self.v0, time)
 
    def alpha_beta(self, det):
        self.alpha, self.beta = self.c.alpha_beta_Pbroad( det, self.sat )
 
    def elastic(self, det):
        self.alpha_beta(det)
        return self.sunits * self.pol * self.dw * ( self.alpha * self.C + self.beta * self.S )

    def inelastic(self, det ):
        self.alpha_beta(det)
        return self.sunits * self.pol * 1.0 * ( self.alpha * self.c.x.size \
                                            + self.beta * 1/4. * self.c.x.size )

    def coh(self,det):
        '''Returns fraction of scattering that is elastic'''
        # Spacing between states 1 and 2 in units of the linewidth
        d12 = 76./5.9 
        # det is with respect to in between 1 and 2 
        det2 = det + d12/2.
        det1 = det - d12/2.
    
        # Saturation parameter
        SP = 2. * self.sat * ( 0.5 / (1.+4.*(det1**2)) + 0.5 / (1.+4.*(det2**2))) 
        
        return 1. / (1. + SP)
  
    def cohOD(self,det):
        '''Returns fraction of scattering that is elastic.
           It takes into account optical density. 
           '''
        # Spacing between states 1 and 2 in units of the linewidth
        d12 = 76./5.9 
        # det is with respect to in between 1 and 2 
        det2 = det + d12/2.
        det1 = det - d12/2.
        
        #Average distance for photon to travel to outside of sample (nm)
        dz = self.c.latsize * self.c.a / 2. 
 
        #On-resonace cross section of single atom (nm^2)
        sigm = 3.*np.pi*(l671/2./np.pi)**2 

        #Lorentzian lineshape
        sigm = sigm * ( 0.5/ (1+4*det1**2) + 0.5/ (1+4*det2**2) ) 
 
        #Density (1 per site)
        nc = 1. / (self.c.a**3) 

        #Optical density 
        OD = nc * sigm * dz
        #print "Optical density at det = %.2f --> OD = %.3f . e^(-OD) = %e" % (det,OD,np.exp(-OD)) 

        return  self.coh(det) * np.exp(-4.*OD) 
        
        

    def combined(self, det):
        return self.elastic(det)*self.coh(det) + \
               self.inelastic(det)*(1.-self.coh(det))

    def combOD(self, det):
        return self.elastic(det)*self.cohOD(det) + \
               self.inelastic(det)*(1.-self.cohOD(det))

  
from uncertainties import unumpy, ufloat
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'}

)




# TESTS to check validity of code 
# TEST DEBYE-WALLER TOF
TEST = True 
if TEST:
    N = 40 
    nafm = 8
    A1 = Kdet(N, nafm, kin, kA1)
    A2 = Kdet(N, nafm, kin, kA2)

    # Calculate Bragg scattering 100 for MIT experiment
    # Q vector for bragg scattering
    lRb = 780. 
    l1064 = 1064.
    Qmit = 4*np.pi/l1064 * vec3.vec3( 0, -1., 0) # 2pi/a
    QUmit = Qmit/abs(Qmit)

    ThetaMIT = np.arccos( abs(Qmit)/2. / ( 2.*np.pi / lRb ) ) 
    print "Bragg angle for MIT 100 experiment = %.2f" % (ThetaMIT*180./np.pi)
    kinMIT = vec3.vec3()
    kinMIT.set_spherical( 2.*np.pi/lRb, np.pi/2., np.pi/2. - ThetaMIT)
    koutMIT = kinMIT + Qmit
   
    AMIT = Kdet(N, nafm, kinMIT, koutMIT)  
    AMIT.c.mass = 87.
    AMIT.v0 = [15.,15.,15.]

    

    # Plot the decay of the debye-waller factor 
    # as a function of time  
    figure3 = plt.figure(figsize=(6.,4.5))
    gs3 = matplotlib.gridspec.GridSpec( 1,1) 
    figure3.suptitle('')
    axDW = plt.subplot( gs3[0,0] )
 
    print "t = %.2f--> DW_t = %e" % (1.0, A2.dw_time(1.0))
    print "t = %.2f--> DW_t = %e" % (100.0, A2.dw_time(100.0))

    
    time =  np.linspace( 0., 30., 200) 
    A2DWtof = np.array([A2.dw_time(t)  for t in time ])
    A1DWtof = np.array([A1.dw_time(t)  for t in time ])
    mitDWtof = np.array([AMIT.dw_time(t) for t in time ])
    
    axDW.plot( time, A2DWtof ,\
                '.-', mec='green', c='green',\
                mew=1.0, ms=3.,\
                marker='o', mfc='limegreen', \
                label='A2 Debye-Waller TOF') 
    
    axDW.plot( time, A1DWtof ,\
                '.-', mec='blue', c='blue',\
                mew=1.0, ms=3.,\
                marker='o', mfc='lightblue', \
                label='A1 Debye-Waller TOF')
 
    axDW.plot( time, mitDWtof ,\
                '.-', mec='black', c='black',\
                mew=1.0, ms=3.,\
                marker='o', mfc='gray', \
                label='MIT Debye-Waller TOF') 
    axDW.grid()
    axDW.set_xlabel('time of flight ($\mu$s)')

    axDWR = axDW.twinx()
    axDWR.plot( time, A2DWtof/A1DWtof ,\
                '.-', mec='red', c='red',\
                mew=1.0, ms=3.,\
                marker='o', mfc='pink', \
                label='A2/A1 Debye-Waller TOF') 
    
    axDW.legend(loc='best',numpoints=1,\
             prop={'size':10}, \
             handlelength=1.1,handletextpad=0.5)
#    axDWR.legend(loc='best',numpoints=1,\
#             prop={'size':10}, \
#             handlelength=1.1,handletextpad=0.5)
    
    gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
    outfile = 'DWtof_plot.png'
    figure3.savefig(outfile, dpi=250)
    exit()
     



TEST = True
if TEST: 
    # Initialize crystal
    #nafms = [2,4,6,7,8,10,12,16,20,24,32,34,38]
    nafms = [2,4,6,7,8]
    S1 = []
    S2 = []
    C1 = []
    C2 = []
    O = []
    N = 40
    for nafm in nafms:
        c = afm.crystal(N, nafm, l1064/2)
        # Define detuning
        det = 0.  # In between
        A1 = Kdet(N,nafm, kin, kA1)
        A2 = Kdet(N,nafm, kin, kA2)  
        inelastic1 = A1.inelastic(det)
        inelastic2 = A2.inelastic(det)
 
        S1.append( [np.mean(A1.S), stats.sem(A1.S)] ) 
        C1.append( [np.mean(A1.C), stats.sem(A1.C)] ) 
        S2.append( [np.mean(A2.S), stats.sem(A2.S)] ) 
        C2.append( [np.mean(A2.C), stats.sem(A2.C)] )

        A10 = ufloat( np.mean(A1.combOD(det)), stats.sem(A1.combOD(det)))
        A20 = ufloat( np.mean(A2.combOD(det)), stats.sem(A2.combOD(det)))
        
        detR = 76/5.9 /2. # On resonance
        A1R = ufloat( np.mean(A1.combOD(detR)), stats.sem(A1.combOD(detR)))
        A2R = ufloat( np.mean(A2.combOD(detR)), stats.sem(A2.combOD(detR)))
 
        ODpt = (A20/A10) / (A2R/A1R)

        print "Nafm = %d --> A2/A1 norm = %.2f" % (nafm, ODpt.nominal_value)

        O.append( [ ODpt.nominal_value, ODpt.std_dev ] ) 
   
          

   
    #Plot the structure factors as a function of Nafm 
    figure2 = plt.figure(figsize=(12.,9.))
    gs2 = matplotlib.gridspec.GridSpec( 2,3) 
    figure2.suptitle('')
    axC1 = plt.subplot( gs2[0,0] )
    axC2 = plt.subplot( gs2[1,0] )
    axS1 = plt.subplot( gs2[0,1] )
    axS2 = plt.subplot( gs2[1,1] )
    axSR = plt.subplot( gs2[0,2] )
    axOnafm = plt.subplot( gs2[1,2] )

    afmdat = np.array(nafms) 
    S1 = np.array(S1)
    C1 = np.array(C1)
    S2 = np.array(S2)
    C2 = np.array(C2)
    O = np.array(O)

    nafmx = np.array(nafms)
    def plotCS( ax, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
        ax.errorbar( nafmx, dat[:,0], yerr=dat[:,1],\
                    capsize=0., elinewidth = 1. ,\
                    fmt='.-', ecolor=lcolor, mec=lcolor, \
                    mew=1.0, ms=3.,\
                    alpha = alpha, \
                    marker=marker, mfc=fcolor, \
                    label=labelstr) 
    plotCS( axS1, S1, 'blue', 'lightblue', 'o', 'S1')  
    plotCS( axC1, C1, 'blue', 'lightblue', 'o', 'C1')  
    plotCS( axS2, S2, 'blue', 'lightblue', 'o', 'S2')  
    plotCS( axC2, C2, 'blue', 'lightblue', 'o', 'C2')
    plotCS( axOnafm, O, 'blue', 'lightblue', 'o', 'Bragg/Diffuse')

    

    plotCS( axSR, S2/S1, 'blue', 'lightblue', 'o', 'S2/S1')

     

    axes = [axS1, axC1, axS2, axC2, axSR, axOnafm]
    for ax in axes:
      ax.grid()
      ax.set_xlabel('Nafm')
    
      ax.legend(loc='best',numpoints=1,\
               prop={'size':12}, \
               handlelength=1.1,handletextpad=0.5)
    
    gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
    outfile = 'CSfactors_plot.png'
    figure2.savefig(outfile, dpi=250)

   
    PRINT = False
    if PRINT: 
        print inelastic1
        print inelastic2
        
        print 
       
        print "A1" 
        print "polsum =",A1.pol
        print "    dw =",A1.dw
        print " alpha =",A1.alpha
        print "  beta =",A1.beta
        print "     S =",np.mean(A1.S)
        print "     C =",np.mean(A1.C)
        print "   coh =",A1.coh(det)
        print "elastic =",np.mean(A1.elastic(det))
        print "inelastic =",A1.inelastic(det)  
        print "combined =",np.mean(A1.combined(det))
        
        print 
        
        print "A2" 
        print "polsum =",A2.pol
        print "    dw =",A2.dw
        print " alpha =",A2.alpha
        print "  beta =",A2.beta
        print "     S =",np.mean(A2.S)
        print "     C =",np.mean(A2.C)
        print "   coh =",A2.coh(det)
        print "elastic =",np.mean(A2.elastic(det))
        print "inelastic =",A2.inelastic(det)  
        print "combined =",np.mean(A2.combined(det))
        
        exit()




cols = 3
rows = 4

figure = plt.figure(figsize=(12.,12.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42) 


lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
m  = [ 'o', 'D', '+', '^', '<', '>']

ax1e = plt.subplot( gs[ 0,0] )
ax2e = plt.subplot( gs[ 1,0] )
axRe = plt.subplot( gs[ 2,0] )

ax1i = plt.subplot( gs[ 0,1] )
ax2i = plt.subplot( gs[ 1,1] )
axR1ie = plt.subplot( gs[ 2,1] )


ax1c = plt.subplot( gs[ 0,2] )
ax2c = plt.subplot( gs[ 1,2] )
axRc = plt.subplot( gs[ 2,2] )

axRcOD = plt.subplot(gs[3,2])

#axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie, ax1c, ax2c, axRc] 
axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie, ax1c, ax2c, axRc, axRcOD] 
for ax in axlist:
  ax.set_xlabel('Detuning ($\Gamma$)')
  # Put vertical lines to show the states
  ax.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)
  ax.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)

ax1e.set_ylabel('ANDOR1')
ax2e.set_ylabel('ANDOR2')
axRe.set_ylabel('A2(elastic)/A1(elastic)')
ax1i.set_ylabel('ANDOR1')
ax2i.set_ylabel('ANDOR2')
axR1ie.set_ylabel('A1(inelastic)/A1(elastic)')

def plotunc( ax, x, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    print "plotting ", labelstr
    ax.errorbar( x, unumpy.nominal_values(dat),\
                yerr=unumpy.std_devs(dat),\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=5.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, \
                label=labelstr) 


#nafms = [4,6,8,10,12,16,20,24,32,34,38,40]
nafms = [2,4,6,7,8,10]

plotafmval = [6,7,8]
plotafmdat = []
for i,nafm in enumerate(nafms):
   
    print "\nWorking on nafm = %d" % nafm 
    # Initialize crystal 
    c = afm.crystal(N, nafm, l1064/2)
    
    print "...calculating C and S"     
    N = 40 
    A1 = Kdet(N,nafm, kin, kA1)
    A2 = Kdet(N,nafm, kin, kA2)  
    print "...calculating detuning curve"

    a1 = []
    a2 = []
    x = np.hstack(( np.linspace(-400,-20,50), \
                    np.linspace(-20,20,51), \
                    np.linspace(20,400,50)))
    for det in x:
        e1 = A1.elastic(det) 
        i1 = A1.inelastic(det)
        c1 = A1.combined(det)
        o1 = A1.combOD(det) 
        
        e2 = A2.elastic(det)
        i2 = A2.inelastic(det) 
        c2 = A2.combined(det)
        o2 = A2.combOD(det)

        a1.append( [np.mean(e1), stats.sem(e1), i1, \
                    np.mean(c1), stats.sem(c1), \
                    np.mean(o1), stats.sem(o1) ] ) 
        a2.append( [np.mean(e2), stats.sem(e2), i2, \
                    np.mean(c2), stats.sem(c2), \
                    np.mean(o2), stats.sem(o2) ] )
    a1 = np.array( a1) 
    a2 = np.array( a2) 
     
    # Make the array with uncertainty of the various cross sections
    # for the uncertainty use col=3 which is the standard error
    a1_e = unumpy.uarray( a1[:,0], a1[:,1] )  
    a1_i = unumpy.uarray( a1[:,2], np.zeros_like(a1[:,0]) )
    a1_c = unumpy.uarray( a1[:,3], a1[:,4] )  
    a1_o = unumpy.uarray( a1[:,5], a1[:,6] )  
    a2_e = unumpy.uarray( a2[:,0], a2[:,1] )  
    a2_i = unumpy.uarray( a2[:,2], np.zeros_like(a2[:,0]) )
    a2_c = unumpy.uarray( a2[:,3], a2[:,4] )  
    a2_o = unumpy.uarray( a2[:,5], a2[:,6] )  

    # PLOT ANDOR 1 
    plotunc( ax1e, x, a1_e, lc[i], 'None', m[i], 'A1 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax1i, x, a1_i, lc[i], 'None', m[i], 'A1 inelastic' + ' Nafm=%d'%nafm)
    plotunc( ax1c, x, a1_c, lc[i], 'None', m[i], 'A1 combined' + ' Nafm=%d'%nafm)
    print "max a1_i = %f" % unumpy.nominal_values(a1_i).max()
 
    # PLOT ANDOR2
    plotunc( ax2e, x, a2_e, lc[i], 'None', m[i], 'A2 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax2i, x, a2_i, lc[i], 'None', m[i], 'A2 inelastic' + ' Nafm=%d'%nafm)
    plotunc( ax2c, x, a2_c, lc[i], 'None', m[i], 'A2 combined' + ' Nafm=%d'%nafm)


    # PLOT ANDOR2/ANDOR1 ELASTIC RATIO
    max0 = unumpy.nominal_values(a2_e/a1_e).max()
    plotunc( axRe, x, a2_e/a1_e/max0 , lc[i], fc[i],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./max0))

    # PLOT ANDOR1 INELATIC / ANDOR1 ELASTIC
    max0 = unumpy.nominal_values(a1_i/a1_e).max()
    plotunc( axR1ie, x, a1_i/a1_e, lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)

    # PLOT ANDOR2/ANDOR1 COMBINED RATIO
    plotunc( axRc, x, a2_c/a1_c , lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)

    # PLOT ANDOR2/ANDOR1 COMBINED PLUS OD RATIO
    plotunc( axRcOD, x, a2_o/a1_o , lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)
    if nafm in plotafmval:
        plotafmdat.append((x,a2_o/a1_o,i,nafm)) 

 


ax1e.set_yscale('log')
ax2e.set_yscale('log')
ax1i.set_yscale('log')
ax2i.set_yscale('log')
ax1c.set_yscale('log')
ax2c.set_yscale('log')
#axR.set_yscale('log')

lims=100.
for ax in axlist:
  ax.set_xlim(-lims,lims)
  ax.legend(loc='best', numpoints=1, prop={'size':5})

Rlims=500.
axRe.set_xlim(-Rlims,Rlims)
Rlims=1000.
axR1ie.set_xlim(-Rlims,Rlims)
Rlims=20.
axRc.set_xlim(-Rlims,Rlims)
axRcOD.set_xlim(-Rlims,Rlims)

gs.tight_layout(figure, rect=[0.,0.,1.0,0.91])
#plt.show()
figure.savefig('detuning.png', dpi=200)
pylab.clf()



########### PLOT SIMULATION ALONGSIDE DATA
import statdat
#Fetch data
detdat={}
detdat['130729_detuning'] = { \
                  'label':'7$E_{r}$  July 29',\
                  'dir':'/lab/data/app3/2013/1307/130729/',\
                  'shots':'3243:3258,3259:3357,-3325,3364:3379',\
                  'ec':'black','fc':'black',\
                  }

datakeys = ['DIMPLELATTICE:imgdet', 'ANDOR1EIGEN:signal' , 'ANDOR2EIGEN:signal', 'HHHEIGEN:andor2norm', 'DIMPLELATTICE:force_lcr3' ]

print
print "STARTING DATA + SIMULATION FIGURE"
print "Fetching data..."
for k in detdat.keys():
    try:
        detdat[k]['data'] = np.loadtxt(k+'.dat')
        print "Loaded %s succesfully." % k 
    except:
        print k 
        data, errmsg, rawdat = qrange.qrange_eval( detdat[k]['dir'], detdat[k]['shots'], datakeys) 
        np.savetxt(k+'.dat', data)
        detdat[k]['data'] = data

print "Done.\n"
        

########### PREPARE FIGURE AND AXIS FOR DATA PLOT
figure1 = plt.figure(figsize=(4.8,3.6))
gs1 = matplotlib.gridspec.GridSpec( 1,1) 
figure1.suptitle('')
ax1 = plt.subplot( gs1[0,0] )

########### PLOT THE SIMULATION RESULTS
def plotsim( ax, x, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    #Find the value at resonance
    res =  abs(abs(x) - 6.44) < 0.1
    datres = dat[res]
    ratiores = np.mean( unumpy.nominal_values(dat[res]))
    normTOF = 0.4/0.66
    ax.errorbar( x, normTOF*unumpy.nominal_values(dat)/ratiores,\
                yerr=0.*normTOF*unumpy.std_devs(dat)/ratiores,\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, \
                label=labelstr) 
for d in plotafmdat:
    plotsim( ax1, d[0], d[1], lc[d[2]], fc[d[2]],\
             'o', 'Nafm=%d'%d[3], alpha=0.5)

# Base is 0.5 TOF
base = 0.656
base_err = 0.068
base1 = base
## Base is on resonance
#base = 0.45
#base_err = 0.03
#base1 = base

ax1.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor='limegreen', alpha=0.6, linewidth=0)
ax1.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor='limegreen', alpha=0.6, linewidth=0)
ax1.axhspan(  (base-base_err)/base1, (base+base_err)/base1, facecolor='gray', alpha=0.6, linewidth=0)

for k in sorted(detdat.keys()):
    dat = detdat[k]['data'][:,(0,3,4)]
    
    # Remove baslines from the data
    dat = dat[ dat[:,2] == -1. ] 
    

    dat = statdat.statdat( dat, 0, 1 )
    ax1.errorbar( (dat[:,0]+131.)/5.9, dat[:,1]/base1, yerr=dat[:,3]/base1,\
                  capsize=0., elinewidth=1.,\
                  fmt='.', ecolor=detdat[k]['ec'], mec=detdat[k]['ec'],\
                  mew=1.0, marker='o', mfc=detdat[k]['fc'],\
                  alpha = 1.0,\
                  #label=detdat[k]['label']) 
                  label='$7E_{r}$, $a_{s}=190a_{0}$') 

#ax1.errorbar( dat[:,0], dat[:,1]/zero, yerr=dat[:,3]/zero, fmt='.', ecolor=edgecolor, mec=edgecolor, mew=1.0, marker='o', mfc=fillcolor) 


ax1.grid()
ax1.set_xlabel('Detuning ($\Gamma$)')
ax1.set_ylabel('Bragg / Diffuse',ha='center',labelpad=20)
ax1.set_xlim(-15.,15.)
#ax1.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(0.1) )
#ax1.set_ylim(0.92, 1.48)

#ax1.legend(loc='upper left',numpoints=1,\
#           prop={'size':7}, \
#           handlelength=1.1,handletextpad=0.5)

ax1.legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':7}, \
           handlelength=1.1,handletextpad=0.5)


gs1.tight_layout(figure1, rect=[0,0.0,0.82,0.96])
outfile = 'detuning_queue.png'
figure1.savefig(outfile, dpi=250)

