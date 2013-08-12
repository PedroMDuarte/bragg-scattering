import numpy as np
import vec3
import pylab


import braggvectors as bv
import afm

 
N = 40 

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})
cols = 2
rows = 2
figure = plt.figure(figsize=(7.,7.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42)

axinp = plt.subplot( gs[ 0,0] )
axinp.set_title("Input rocking curve")
axout = plt.subplot( gs[ 1,0] )
axout.set_title("Output rocking curve")

def plotRock(ax, xdat, ydat, yerrdat, labelstr, lc, fc):
    xdat = np.array(xdat)
    ydat = np.array(ydat)
    ydat = ydat - 1.
    yerrdat = np.array(yerrdat)
    norm = ydat.max()
    ax.errorbar( xdat, ydat/norm, yerr=yerrdat/norm, \
               capsize=0., elinewidth = 1. ,\
               fmt='.', ecolor=lc, mec=fc, \
               mew=1., ms=5.,\
               marker='o', mfc='lightblue', \
               label=labelstr) 
 
lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
nc = len(lc)

#nafms=[2,4,6,7,8,9,10,12,16,20,24,32,34,38]
nafms=[8,16,24,32]
for i,nafm in enumerate(nafms): 
    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
    Npts = 20 # was 20
    inangles = np.linspace(-6000./nafm**1.5,6000./nafm**1.5,Npts)
    outangles = np.linspace(-6000./nafm**1.5,6000./nafm**1.5,Npts)

    inrock =[]
    inrockerr = []
    for angle in inangles:
        A1.set_kvectors( bv.kinput(angle), bv.kA1, bv.kipol ) 
        A2.set_kvectors( bv.kinput(angle), bv.kA2, bv.kipol ) 
        Nr = 20
        sig =  A2.sigma_coh_det( Nr, 0., 0.) / A1.sigma_coh_det( Nr, 0., 0.) 
        inrock.append( sig.nominal_value) 
        inrockerr.append( sig.std_dev)
    plotRock(axinp, inangles, inrock, inrockerr, 'nafm=%d'%nafm, lc[i%nc], fc[i%nc])

    outrock =[]
    outrockerr =[]
    for angle in outangles:
        A1.set_kvectors( bv.kin, bv.kA1, bv.kipol ) 
        A2.set_kvectors( bv.kin, bv.koutput(angle), bv.kipol ) 
        Nr = 20
        sig = A2.sigma_coh_det( Nr, 0., 0.) / A1.sigma_coh_det( Nr, 0., 0.) 
        outrock.append( sig.nominal_value)
        outrockerr.append( sig.std_dev)
    plotRock(axout, outangles, outrock, outrockerr, 'nafm=%d'%nafm, lc[i%nc], fc[i%nc])


axinp.grid()
axout.grid()


axinp.legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':7}, \
           handlelength=1.1,handletextpad=0.5)
axout.legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':7}, \
           handlelength=1.1,handletextpad=0.5)

gs.tight_layout(figure, rect=[0,0.0,0.82,0.96])
outfile = 'Rocking.png'
figure.savefig(outfile, dpi=250)
        







