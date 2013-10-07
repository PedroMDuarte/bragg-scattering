import numpy as np
import vec3
import pylab


import braggvectors as bv
import afm

 
N = 10

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})
cols = 1
rows = 3
figure = plt.figure(figsize=(7.5,9.))
figure.text(0.5,0.95,'CAUTION, THESE ARE SIMULATIONS. NOT REAL DATA',fontsize=15,ha='center')
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols)
#wspace=0.6, hspace=0.42)

axSig = plt.subplot( gs[ 0,0] )
axSig.set_title("Bragg chamber rocking curve, N=%d"%N)
axA1 = plt.subplot( gs[ 1,0] )
axA2 = plt.subplot( gs[ 2,0] )

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
 
lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
nc = len(lc)

#nafms=[2,4,6,7,8,9,10,12,16,20,24,32,34,38]
#nafms=[6,8,10,16,24]
nafms=[2,6,8]
rockingW=[]
for i,nafm in enumerate(nafms): 
    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
    # Four quadrants to consider lens apperture
    q = bv.kout_quadrants
    A2q0 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,q[0]))
    A2q1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,q[1]))
    A2q2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,q[2]))
    A2q3 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,q[3]))
   
    # Realizations for random spins
    Nr = int(48  /nafm)
 
    alim = 180.
    Npts = 92 # was 
    inangles = np.linspace(-alim,alim,Npts) + 2.0
    print 'Nafm=%d' % nafm + '  inangles = ',inangles

    fname = 'CHAMBER_%.2f_%d_%d_%d'%(alim,Npts,nafm,N)
    try:
        rock_sig = loadtxtU( 'rockingdat/rock_sig'+fname)
        rock_A1 = loadtxtU( 'rockingdat/rock_A1'+fname)
        rock_A2 = loadtxtU( 'rockingdat/rock_A2'+fname)
        print "Loaded input ",fname," successfully."
    except: 
        rock_sig =[]
        rock_A1 = []
        rock_A2 = []
        
        for angle in inangles:
            print "Working on angle = ",angle
            k = bv.kinchamber( angle )
            A1.set_kvectors( k, bv.kA1, bv.kipol ) 
            A2.set_kvectors( k, bv.kA2, bv.kipol )
            # Lens quadrants
            A2q0.set_kvectors( k, q[0], bv.kipol )
            A2q1.set_kvectors( k, q[1], bv.kipol )
            A2q2.set_kvectors( k, q[2], bv.kipol )
            A2q3.set_kvectors( k, q[3], bv.kipol )
        
            #A2_quad_average = ( A2.sigma_coh_det( Nr, 0., 0.) \
            A2_quad_average = (\
                                A2q0.sigma_coh_det( Nr, 0., 0.) \
                               + A2q1.sigma_coh_det( Nr, 0., 0.) \
                               + A2q2.sigma_coh_det( Nr, 0., 0.) \
                               + A2q3.sigma_coh_det( Nr, 0., 0.) \
                              ) / 4. 
                              
            A1sig = A1.sigma_coh_det( Nr, 0., 0.) 
            sig =  A2_quad_average / A1sig
            rock_sig.append( sig )
            rock_A1.append( A1sig )
            rock_A2.append( A2_quad_average) 
 
        rock_sig = np.array( rock_sig) 
        rock_A1 = np.array( rock_A1) 
        rock_A2 = np.array( rock_A2) 
        savetxtU( 'rockingdat/rock_sig'+fname, rock_sig)
        savetxtU( 'rockingdat/rock_A1'+fname, rock_A1)
        savetxtU( 'rockingdat/rock_A2'+fname, rock_A2)

    inW,inWerr=plotRock(axSig, inangles, rock_sig,'$N_{\mathrm{AFM}}$=%02d'%nafm, lc[i%nc], fc[i%nc])
    inW,inWerr=plotRock(axA1,  inangles, rock_A1,'$N_{\mathrm{AFM}}$=%02d'%nafm, lc[i%nc], fc[i%nc])
    inW,inWerr=plotRock(axA2,  inangles, rock_A2,'$N_{\mathrm{AFM}}$=%02d'%nafm, lc[i%nc], fc[i%nc])
    

    rockingW.append( [inW, inWerr])

rockingW = np.array(rockingW) 

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


allax = [axSig, axA1, axA2]#, axSigW ]
for ax in allax:
    ax.grid()

#axSig.set_ylim(0.,10)
#axout.set_ylim(0.,10)


axSig.legend(bbox_to_anchor=(1.0,1.0),loc='upper right',numpoints=1,prop={'size':5}, \
           handlelength=1.1,handletextpad=0.5)

gs.tight_layout(figure, rect=[0,0.0,1.0,0.96])
outfile = 'RockingChamberN%d.png' % N
figure.savefig(outfile, dpi=250)
        







