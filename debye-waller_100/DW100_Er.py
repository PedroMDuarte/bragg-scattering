
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
Er = 20. 
A1.set_v0([Er, Er, Er])
A2.set_v0([Er, Er, Er])
MANTA.set_v0([Er, Er, Er])

#A1.set_pbragg(0.0002)
#MANTA.set_pbragg(0.0002)




# Plot the decay of the debye-waller factor 
# as a function of time  
figure3 = plt.figure(figsize=(6.,4.5))
gs3 = matplotlib.gridspec.GridSpec( 1,1) 
figure3.suptitle('')
axDW = plt.subplot( gs3[0,0] )

MANTAdw = []
A1dw = []
A2dw = []

v0 = np.linspace( 0.1, 50., 200)
MANTAdw = np.array( [ MANTA.dw_( v, 0. ) for v in v0 ] ) 
A1dw = np.array( [ A1.dw_( v, 0. ) for v in v0 ] ) 
A2dw = np.array( [ A2.dw_( v, 0. ) for v in v0 ] ) 

axDW.plot( v0, MANTAdw ,\
            '-', c='green', lw=2.,\
            label='Input = 100, Cam = MANTA') 

axDW.plot( v0, A1dw ,\
            '-', c='blue', lw=2.,\
            label='Input = 100, Cam = A1')
 
axDW.plot( v0, A2dw ,\
            '-', c='magenta', lw=2.,\
            label='Input = 100, Cam = A2') 

axDW.grid()
axDW.set_xlabel('Lattice Depth ($E_{R}$)')
axDW.set_ylabel('Debye-Waller factor')

axDW.legend(loc='best',numpoints=1,\
         prop={'size':10}, \
         handlelength=1.1,handletextpad=0.5)
#axDWR.legend(loc='best',numpoints=1,\
#         prop={'size':10}, \
#         handlelength=1.1,handletextpad=0.5)
    
gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'DW100v0_plot.png' 
figure3.savefig(outfile, dpi=250)
     


