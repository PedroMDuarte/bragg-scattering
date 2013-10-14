
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

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

A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA2))
MANTA = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kMANTA))

A1h = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
A2h = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
MANTAh = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kMANTA))

A1n1 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[1],bv.kA1))
A2n1 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[1],bv.kA2))
MANTAn1 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[1],bv.kMANTA))

A1n5 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[5],bv.kA1))
A2n5 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[5],bv.kA2))
MANTAn5 = afm.crystal(N, nafm, bv.l1064/2, (bv.ksquad[5],bv.kMANTA))

# Plot the decay of the debye-waller factor 
# as a function of time  
figure3 = plt.figure(figsize=(5.,4.))
gs3 = matplotlib.gridspec.GridSpec( 1,1) 
figure3.suptitle('')
axA1 = plt.subplot( gs3[0,0] ) 

angle  = np.linspace( -180., 180., 200)

dipA1 = []
dipA2 = []
crys = afm.crystal( N, nafm, bv.l1064/2, (bv.kin, bv.kA2) ) 
for a in angle: 
    crys.set_kvectors( bv.kinput(a), bv.kA1, bv.kipol ) 
    dipA1.append( crys.pol() ) 
    
    crys.set_kvectors( bv.kinput(a), bv.kA2, bv.kipol ) 
    dipA2.append( crys.pol() ) 
print dipA2
dipA1 = np.array( dipA1)    
dipA2 = np.array( dipA2)    

ax.plot( angle, dipA1, c='red', lw=1.5, label='Cam = A1' )
ax.plot( angle, dipA2, c='green', lw=1.5, label='Cam = A2' )

ax.grid()
ax.set_xlabel('Angle along chamber (degrees)')
ax.set_ylabel('Dipole pattern strength')
#ax.set_xlim(0., 50.)
#ax.set_ylim(0., 1.)
   


gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'dipole.png' 
figure3.savefig(outfile, dpi=250)
     


