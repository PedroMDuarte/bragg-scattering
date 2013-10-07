
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
    y = np.array( [crys.dw_v0( v ) for v in v0 ] ) 
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
   


##MANTAdw = np.array( [ MANTA.dw_v0( v ) for v in v0 ] ) 
##A1dw = np.array( [ A1.dw_v0( v ) for v in v0 ] ) 
##A2dw = np.array( [ A2.dw_v0( v ) for v in v0 ] ) 
##
##MANTAdwH = np.array( [ MANTAh.dw_v0( v ) for v in v0 ] ) 
##A1dwH = np.array( [ A1h.dw_v0( v ) for v in v0 ] ) 
##A2dwH = np.array( [ A2h.dw_v0( v ) for v in v0 ] )
##
##MANTAdwN1 = np.array( [ MANTAn1.dw_v0( v ) for v in v0 ] ) 
##A1dwN1 = np.array( [ A1n1.dw_v0( v ) for v in v0 ] ) 
##A2dwN1 = np.array( [ A2n1.dw_v0( v ) for v in v0 ] )
##
##MANTAdwN5 = np.array( [ MANTAn5.dw_v0( v ) for v in v0 ] ) 
##A1dwN5 = np.array( [ A1n5.dw_v0( v ) for v in v0 ] ) 
##A2dwN5 = np.array( [ A2n5.dw_v0( v ) for v in v0 ] )
##
##
##
##axDW.plot( v0, MANTAdw ,\
##            ':', c='green', lw=lw,\
##            label='Input = 100, Cam = MANTA') 
##
##axDW.plot( v0, A1dw ,\
##            ':', c='blue', lw=lw,\
##            label='Input = 100, Cam = A1')
## 
##axDW.plot( v0, A2dw ,\
##            ':', c='red', lw=lw,\
##            label='Input = 100, Cam = A2') 
##
##axDW.plot( v0, MANTAdwH ,\
##            '-', c='green', lw=lw,\
##            label='Input = HHH, Cam = MANTA') 
##
##axDW.plot( v0, A1dwH ,\
##            '-', c='blue', lw=lw,\
##            label='Input = HHH, Cam = A1')
## 
##axDW.plot( v0, A2dwH ,\
##            '-', c='red', lw=lw,\
##            label='Input = HHH, Cam = A2') 
##
##axDW.plot( v0, MANTAdwN1 ,\
##            '--', c='green', lw=lw,\
##            label='Input = #1, Cam = MANTA') 
##axDW.plot( v0, A1dwN1 ,\
##            '--', c='blue', lw=lw,\
##            label='Input = #1, Cam = A1')
##axDW.plot( v0, A2dwN1 ,\
##            '--', c='red', lw=lw,\
##            label='Input = #1, Cam = A2') 
##
##axDW.plot( v0, MANTAdwN5 ,\
##            '-.', c='green', lw=lw,\
##            label='Input = #5, Cam = MANTA') 
##axDW.plot( v0, A1dwN5 ,\
##            '-.', c='blue', lw=lw,\
##            label='Input = #5, Cam = A1')
##axDW.plot( v0, A2dwN5 ,\
##            '-.', c='red', lw=lw,\
##            label='Input = #5, Cam = A2') 
##


#axDWR.legend(loc='best',numpoints=1,\
#         prop={'size':10}, \
#         handlelength=1.1,handletextpad=0.5)
    
gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'DWv0_plot.png' 
figure3.savefig(outfile, dpi=250)
     


