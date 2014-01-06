
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import numpy as np
import vec3
import pylab
from scipy import stats

import braggvectors as bv
import afm
from uncertainties import unumpy, ufloat

import statdat

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})

N = 40 
nafm = 8

A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
Er = 20. 
A1.set_v0([Er, Er, Er])
A2.set_v0([Er, Er, Er])

#A1.set_pbragg(0.0002)
#A2.set_pbragg(0.0002)

# Calculate Bragg scattering 100 for MIT experiment
# Q vector for bragg scattering
lRb = 780. 
l1064 = 1064.
Qmit = 4*np.pi/l1064 * vec3.vec3( 0, -1., 0) # 2pi/a
QUmit = Qmit/abs(Qmit)

ThetaMIT = np.arccos( abs(Qmit)/2. / ( 2.*np.pi / lRb ) )
print 
print "Bragg angle for MIT 100 experiment = %.2f" % (ThetaMIT*180./np.pi)
kinMIT = vec3.vec3()
kinMIT.set_spherical( 2.*np.pi/lRb, np.pi/2., np.pi/2. - ThetaMIT)
koutMIT = kinMIT + Qmit

AMIT = afm.crystal(N, nafm, bv.l1064/2, (kinMIT,koutMIT))
AMIT.mass = 87.
AMIT.set_v0( [15.,15.,15.])



# Plot the decay of the debye-waller factor 
# as a function of time  
figure3 = plt.figure(figsize=(6.,4.5))
gs3 = matplotlib.gridspec.GridSpec( 1,1) 
figure3.suptitle('')
axDW = plt.subplot( gs3[0,0] )

time =  np.linspace( 0., 30., 200) 
A2DWtof = np.array([A2.dw_(Er,t)  for t in time ])
A1DWtof = np.array([A1.dw_(Er,t)  for t in time ])
mitDWtof = np.array([AMIT.dw_(15.,t) for t in time ])

axDW.plot( time, A2DWtof ,\
            '.-', mec='green', c='green',\
            mew=1.0, ms=3.,\
            marker='o', mfc='limegreen', \
            label='A2 Debye-Waller TOF, 20$E_{r}$') 

axDW.plot( time, A1DWtof ,\
            '.-', mec='blue', c='blue',\
            mew=1.0, ms=3.,\
            marker='o', mfc='lightblue', \
            label='A1 Debye-Waller TOF, 20$E_{r}$')

axDW.plot( time, mitDWtof ,\
            '.-', mec='black', c='black',\
            mew=1.0, ms=3.,\
            marker='o', mfc='gray', \
            label='MIT Debye-Waller TOF, 15$E_{r}$') 
axDW.grid()
axDW.set_xlabel('time of flight ($\mu$s)')

axDW.legend(loc='best',numpoints=1,\
         prop={'size':10}, \
         handlelength=1.1,handletextpad=0.5)
#axDWR.legend(loc='best',numpoints=1,\
#         prop={'size':10}, \
#         handlelength=1.1,handletextpad=0.5)
    
gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'DWtof_plot_%.1fEr_wMIT.png' % Er
figure3.savefig(outfile, dpi=250)
     


