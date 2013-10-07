
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
rows = 5
cols = 3
figure = plt.figure(figsize=(12.,17.5))
figure.text(0.5,0.95,'CAUTION, THESE ARE SIMULATIONS. NOT REAL DATA',fontsize=15,ha='center')
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols)
#wspace=0.6, hspace=0.42)

ax = [] 
for i in range(rows):
    for j in range(cols):
        ax.append( plt.subplot( gs[i,j] ) ) 

# 0 is RATIO
# 1 is A1 
# 2 is A2 
nplots = 15

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



def plotRock(ax, xdat, ydat, labelstr, lc, fc):
    xdat = np.array(xdat)
    ydatval = unumpy.nominal_values(ydat) 
    ydaterr = unumpy.std_devs(ydat) 

    # We will plot how much above the norm is A2/A1
    # in units of the norm 
    #ydat = ydat - 1.
    #maxy = ydat.max()
    maxy = 1.0
    ax.errorbar( xdat, ydatval/maxy, yerr=ydaterr/maxy, \
               capsize=0., elinewidth = 1. ,\
               fmt='.', ecolor=lc, mec=fc, \
               mew=1., ms=5.,\
               marker='o', mfc='lightblue', \
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
 
lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
nc = len(lc)

#nafms=[2,4,6,7,8,9,10,12,16,20,24,32,34,38]
#nafms=[6,8,10,16,24]
N = 40
nafms=[9,12]
rockingW=[]
#for i,nafm in enumerate(nafms): 
#    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
#    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
#
#    # Realizations for random spins
#    Nr = int(96  /nafm)
# 
#    alim = 30.
#    Npts = 32 
#    inangles = np.linspace(-alim,alim,Npts) + 2.0
#    print 'Nafm=%d' % nafm + ' +/-',alim
#
#    fname = 'CHAMBER_%.2f_%d_%d_%d_%d'%(alim,Npts,nafm,N,Nr)
#    try:
#        rock = []
#        for p in range(nplots):
#            rock.append( loadtxtU( 'rockingdat/%02d'%p+fname ) ) 
#        print "Loaded input ",fname," successfully."
#    except:
#        rock =  [ [] for p in range(nplots)]  
#  
#        for angle in inangles:
#            print "Working on angle = ",angle
#            k = bv.kinchamber( angle )
#            A1.set_kvectors( k, bv.kA1, bv.kipol ) 
#            A2.set_kvectors( k, bv.kA2, bv.kipol )
#
#            tof = 6. # tof in us
#
#            a1 = A1.sigma_coh_det( Nr, 0., 0. )
#            a1tof = A1.sigma_coh_det( Nr, 0., tof  ) 
#
#            a2 = A2.sigma_coh_det( Nr, 0., 0. )
#            a2tof = A2.sigma_coh_det( Nr, 0., tof  )
#  
#            vals = [ a2/a1 , a2, a1,  \
#                     a2tof/a1tof, a2tof, a1tof, \
#                     (a2/a1) / (a2tof/a1tof) , a2/a2tof, a1/a1tof ]
#            for p in range(nplots):
#                rock[p].append( vals[p] ) 
#
#        rock = [  np.array( rock[p] ) for  p in range(nplots) ] 
#        for p in range(nplots):
#            savetxtU( 'rockingdat/%02d'%p+fname, rock[p] ) 
#
#    lab = '$N_{\mathrm{AFM}}$=%02d'%nafm
#    for p in range(nplots):
#         plotRock( ax[p], inangles, rock[p], lab, lc[i%nc], fc[i%nc] ) 


################################################
#
#  CALCULATION USING FINITE CORRELATION LENGTH
#
################################################
nafms=[4,8,12]
for i,nafm in enumerate(nafms): 
    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))

    # Size of the sample (should be much larger than Lc)
    Nsize = 24
    alim = 120.
    Npts = 320 
    inangles = np.linspace(-alim,alim,Npts) 
    print 'Nafm=%d' % nafm + '  alim = ',alim

    fname = 'CHAMBER_LCORR_%.2f_%d_%d_%d'%(alim,Npts,nafm,Nsize)
    try:
        rock = []
        for p in range(nplots):
            rock.append( loadtxtU( 'rockingdat/%02d'%p+fname )) 
        print "Loaded input ",fname," successfully."
    except:
        rock =  [ [] for p in range(nplots)] 
  
        for angle in inangles:
            #print "Working on angle = ",angle
            k = bv.kinchamber( angle )
            A1.set_kvectors( k, bv.kA1, bv.kipol ) 
            A2.set_kvectors( k, bv.kA2, bv.kipol )

            tof = 6. # tof in us

            a2 = A2.sigma_coh_det_Lc( Nsize, 0., 0., nafm )
            a2tof = A2.sigma_coh_det_Lc( Nsize, 0., tof, nafm  ) 

            #print 'a2=',a2
            #print 'a2tof=',a2tof

            a1 = A1.sigma_coh_det_Lc( Nsize, 0., 0., nafm)
            a1tof = A1.sigma_coh_det_Lc( Nsize, 0., tof, nafm ) 
  
            vals = [ a2/a1 , a2, a1,  \
                     a2tof/a1tof, a2tof, a1tof, \
                     (a2/a1) / (a2tof/a1tof) , a2/a2tof, a1/a1tof, \
                     ufloat(1.,0.), ufloat(A2.Slcorr,0.), ufloat(A2.Clcorr,0.) ,\
                     ufloat(1.,0.), ufloat(A1.Slcorr,0.), ufloat(A2.Clcorr,0.) ] 
            for p in range(nplots):
                rock[p].append( vals[p] ) 

        rock = [  np.array( rock[p] ) for  p in range(nplots) ] 
        for p in range(nplots):
            savetxtU( 'rockingdat/%02d'%p+fname, rock[p] ) 

    lab = '$N_{\mathrm{AFM}}$=%02d'%nafm
    for p in range(nplots):
        plotRock( ax[p], inangles, rock[p], lab, lc[i%nc], fc[i%nc] ) 


ylabels = ['a2/a1', 'a2', 'a1', 'a2tof/a1tof', 'a2tof', 'a1tof', \
           'a2/a1 / a2tof/a1tof', 'a2/a2tof', 'a1/a1tof' ,\
           '', 'a2_S', 'a2_C',\
           '', 'a1_S', 'a1_C'  ] 
for p in range(nplots):
    ax[p].set_ylabel( ylabels[p] )   
    ax[p].grid()

# Plot the widths as a function of cristal size

##axSigW = plt.subplot( gs[ 0,1] )
##
##print nafm
##print rockingW[:,0]
##print rockingW[:,1]
##
##axSigW.errorbar( nafms, rockingW[:,0], yerr=rockingW[:,1], \
##           capsize=0., elinewidth = 1. ,\
##           fmt='.', ecolor='black', mec='black', \
##           mew=1., ms=5.,\
##           marker='o', mfc='lightblue')

# Plot the widths that you expect from beam diffraction
# half 1/e^2 angle = lambda / (pi w0) 
# w0 = nafm/2 * sitespacing = nafm/2 * lambda/2
# half 1/e^2 angle = 4 / (pi nafm)
##xdiff = np.linspace(min(nafms),max(nafms), 100)
##Wdiff = (1/np.sqrt(2.))*1000. *   4. / (np.pi*xdiff)
##axSigW.plot( xdiff, Wdiff, '-b',label=r'$4/(\sqrt{2}\pi N_{\mathrm{AFM}})$')
##axSigW.legend(loc='best',numpoints=1,prop={'size':11},handlelength=1.1,handletextpad=0.5)
##
##axSigW.set_xlabel('Nafm')
##axSigW.set_ylabel('Rocking curve\n$1/e^{2}$ half width (mrad)',ha='center',labelpad=16)


#axSig.set_ylim(0.,10)
#axout.set_ylim(0.,10)


#axS.legend(bbox_to_anchor=(1.0,1.0),loc='upper right',numpoints=1,prop={'size':5}, \
#           handlelength=1.1,handletextpad=0.5)

gs.tight_layout(figure, rect=[0,0.0,1.0,0.96])
outfile = 'RockingChamberN%d.png' % N
figure.savefig(outfile, dpi=250)
        







