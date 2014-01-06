
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


# Plot the decay of the debye-waller factor 
# as a function of time  
figure3 = plt.figure(figsize=(10.,8.))
gs3 = matplotlib.gridspec.GridSpec( 2,2) 
figure3.suptitle('')
axHHH = plt.subplot( gs3[0,0] )
ax100 = plt.subplot( gs3[0,1] )
axN1  = plt.subplot( gs3[1,0] )
axN2  = plt.subplot( gs3[1,1] )

v0 = np.linspace( 0.1, 50., 200)

axes = [ axHHH, ax100, axN1, axN2 ] 
kins = [ bv.kin, bv.kin100, bv.ksquad[1], bv.ksquad[5] ]  
instr = [ 'HHH', '100', '#1', '#5' ]  
kouts = [ bv.kA1, bv.kA2, bv.kMANTA ]
camstr = ['A1', 'A2', 'MANTA'] 

lw=1.5 
def plotdw( ax, kin, kout, c='black', lt='-', label=None ): 
    crys = afm.crystal( N, nafm, bv.l1064/2, (kin, kout)) 
    y = np.array( [crys.dw_( v , 0.) for v in v0 ] ) 
    ax.plot( v0, y ,\
            lt, c=c, lw=lw,\
            label=label)

colors = ['red', 'green', 'blue']

for i,ax in enumerate(axes):
    for j,ko in enumerate(kouts):
        plotdw( ax , kins[i], ko, c=colors[j], label='Input = %s, Cam=%s'%(instr[i],camstr[j]) ) 
    ax.legend(loc='best',numpoints=2,\
             prop={'size':8}, \
             handlelength=3.1,handletextpad=0.5)
    ax.grid()
    ax.set_xlabel('Lattice Depth ($E_{R}$)')
    ax.set_ylabel('Debye-Waller factor')
    ax.set_title('kInput = %s'%instr[i])
    ax.set_xlim(0., 50.)
    ax.set_ylim(0., 1.)
   


gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'DWv0_plot.png' 
figure3.savefig(outfile, dpi=250)
     


