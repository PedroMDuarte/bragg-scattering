
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


axes = [ axHHH, ax100, axN1, axN2 ] 
kins = [ bv.kin, bv.kin100, bv.ksquad[1], bv.ksquad[5] ]  
instr = [ 'HHH', '100', '#1', '#5' ]  
kouts = [ bv.kA1, bv.kA2, bv.kMANTA ]
camstr = ['A1', 'A2', 'MANTA']

v0  = 50. # lattice depth in Er
tof = 0.  # tof in usec 
print "\nv0  = %.1f Er"%v0
print "tof = %.2f us\n"%tof

lw=1.5 
def printdw( kin, kout, c='black', lt='-', label=None ):
    crys = afm.crystal( N, nafm, bv.l1064/2, (kin, kout)) 
    dw = crys.dw_( v0 , tof)
    print "%s ==> DW  = %.2f"%(label,dw )

colors = ['red', 'green', 'blue']

for i,ax in enumerate(axes):
    for j,ko in enumerate(kouts):
        printdw( kins[i], ko, c=colors[j], label='Input = %s, Cam=%s'%(instr[i],camstr[j]) ) 
     


