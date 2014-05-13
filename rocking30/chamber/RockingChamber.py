
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import numpy as np
import vec3
import pylab


import braggvectors as bv
import afm

from uncertainties import unumpy, ufloat

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})
rows = 3
cols = 3
figure = plt.figure(figsize=(12.,8.))
figure.text(0.5,0.95,'CAUTION, THESE ARE SIMULATIONS. NOT REAL DATA',fontsize=15,ha='center')
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.7,\
                                   left=0.07,top=0.92,right=0.95,bottom=0.04)

ax = [] 
for i in range(rows):
    for j in range(cols):
        ax.append( plt.subplot( gs[i,j] ) ) 

# 0 is RATIO
# 1 is A1 
# 2 is A2 
nplots = 9

import fitlibrary
import uncertainties
from uncertainties import unumpy

import csv
# Load and Save arrays with uncertainties
def savetxtU( fname, X ):
    np.savetxt( fname, X, fmt='%r') 
def loadtxtU( fname ):
    with open(fname) as f:
        reader = csv.reader( f, delimiter=' ', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    converters = dict.fromkeys( range(num_cols), uncertainties.ufloat_fromstr)
    return np.loadtxt( fname, converters=converters, dtype=object) 



def plotRock(ax, xdat, ydat, labelstr, lc, fc, marker='o', Normalize=False):
    xdat = np.array(xdat)
    ydatval = unumpy.nominal_values(ydat) 
    ydaterr = unumpy.std_devs(ydat) 

    # We will plot how much above the norm is A2/A1
    # in units of the norm 
    if Normalize:
        ydatval = ydatval - 1.
        maxy = np.amax( ydatval)
        print "maxy =", maxy
    else:
        maxy = 1.0

    ax.errorbar( xdat, ydatval/maxy, yerr=ydaterr/maxy, \
               capsize=0., elinewidth = 1. ,\
               fmt='.', ecolor=lc, mec=lc, \
               mew=1., ms=5.,\
               marker=marker, mfc='None', \
               label=labelstr+', $B_{\mathrm{tof}}$=%.2g'%(maxy+1))
    # Fit data with a Gaussian
    fitdat = np.transpose( np.vstack( (xdat,ydat/maxy)))
    p0 = [1.0, 0., 10., 0.1]
    fun = fitlibrary.fitdict['Gaussian'].function
    ##pG, errorG = fitlibrary.fit_function( p0,  fitdat, fun)
    #print "Fitting with Gaussian:"
    #print pG
    #print errorG
    ##fitX, fitY = fitlibrary.plot_function(pG, \
    ##             np.linspace( xdat.min(), xdat.max(),120),fun) 
    ##ax.plot(  fitX, fitY, '-', c=lc, lw=1.0)
    ##return np.sqrt(2.)*pG[2], np.sqrt(2.)*errorG[2]
    return 1.,1.


################################################
#
#  CALCULATION USING FINITE SAMPLE WITH AFM CORE
#
################################################
 
lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
m  = [ 'o', 'D', '+', '^', '<', '>','o', 'D', '+', '^', '<', '>']
nc = len(lc)

#nafms=[2,4,6,7,8,9,10,12,16,20,24,32,34,38]
#nafms=[6,8,10,16,24]
N = 40
nafms=[6,7,8,9,10,11,12]
for i,nafm in enumerate(nafms): 
    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
   
    Nr = 260
    alim = 30.
    Npts = 33 
    inangles = np.linspace(-alim,alim,Npts) 
    print 'Nafm=%d' % nafm + '  alim = ',alim

    fname = 'CHAMBER__%.2f_%d_%d_%d'%(alim,Npts,nafm,Nr)
    try:
        rock = []
        for p in range(nplots):
            rock.append( loadtxtU( 'rockingdat/%02d'%p+fname )) 
        print "Loaded input ",fname," successfully."
    except:
        rock =  [ [] for p in range(nplots)] 
  
        for angle in inangles:
            print angle,
            sys.stdout.flush()
            k = bv.kinchamber( angle )
            A1.set_kvectors( k, bv.kA1, bv.kipol ) 
            A2.set_kvectors( k, bv.kA2, bv.kipol )

            tof = 100. # tof in us

            a2    = A2.I_(Nr=Nr, detuning=0., tof=0. )
            a2tof = A2.I_(Nr=Nr, detuning=0., tof=tof) 
            a1    = A1.I_(Nr=Nr, detuning=0., tof=0. )
            a1tof = A1.I_(Nr=Nr, detuning=0., tof=tof)
 
            vals = [ a2/a1 , a2, a1,  \
                     a2tof/a1tof, a2tof, a1tof, \
                     (a2/a1) / (a2tof/a1tof) , a2/a2tof, a1/a1tof ] 
            for p in range(nplots):
                rock[p].append( vals[p] ) 

        rock = [  np.array( rock[p] ) for  p in range(nplots) ] 
        for p in range(nplots):
            savetxtU( 'rockingdat/%02d'%p+fname, rock[p] ) 

    lab = '$N_{\mathrm{AFM}}$=%02d'%nafm
    for p in range(nplots):
        if p in [6,7,8]:
            Normalize = True
        else:
            Normalize = False
        plotRock( ax[p], inangles, rock[p], lab, lc[i%nc], fc[i%nc], \
                  marker=m[i%nc], Normalize=Normalize ) 


ylabels = ['a2/a1', 'a2', 'a1', 'a2tof/a1tof', 'a2tof', 'a1tof', \
           'a2/a1 / a2tof/a1tof', 'a2/a2tof', 'a1/a1tof' ,\
           '', 'a2_S', 'a2_C',\
           '', 'a1_S', 'a1_C'  ] 
for p in range(nplots):
    ax[p].set_ylabel( ylabels[p] )   
    ax[p].grid()
    
    ax[p].legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':5}, \
           handlelength=1.1,handletextpad=0.5)


#gs.tight_layout(figure, rect=[0,0.0,1.0,0.96])
outfile = 'RockingChamberN%d.png' % N
figure.savefig(outfile, dpi=250)
        







