import numpy as np
import vec3
import pylab


import braggvectors as bv
import afm

 
N = 12 

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})
cols = 2
rows = 1
figure = plt.figure(figsize=(7.5,3.75))
figure.text(0.5,0.95,'CAUTION, THESE ARE SIMULATIONS. NOT REAL DATA',fontsize=15,ha='center')
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols)
#wspace=0.6, hspace=0.42)

axinp = plt.subplot( gs[ 0,0] )
axinp.set_title("Bragg squad rocking curve, N=%d"%N)

import fitlibrary

def plotRock(ax, xdat, ydat, yerrdat, norm, labelstr, lc, fc):
    xdat = np.array(xdat)
    ydat = np.array(ydat) /norm
    yerrdat = np.array(yerrdat)/norm

    # We will plot how much above the norm is A2/A1
    # in units of the norm 
    ydat = ydat - 1.
    maxy = ydat.max()
    ax.errorbar( xdat, ydat/maxy, yerr=yerrdat/maxy, \
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
nafms=[6]
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
    Nr = int(72 * 8 /nafm)
 
    # Long time-of-flight normalization factor
    A1.set_kvectors( bv.kin, bv.kA1, bv.kipol ) 
    A2.set_kvectors( bv.kin, bv.kA2, bv.kipol ) 
    normtof = A2.sigma_coh_det( Nr, 0., 1e4) / A1.sigma_coh_det( Nr, 0., 1e4)
    normtof = normtof.nominal_value
    if i == 0:
        print "\nA2/A1(t0f=10000us) = ",normtof
        print
  

    fname = 'BraggSQUAD_%d'%(nafm)
    try:
        inrock = np.loadtxt( 'rockingdat/inrock'+fname)
        inrockerr = np.loadtxt( 'rockingdat/inrockerr'+fname)
        print "Loaded input ",fname," successfully."
    except: 
        inrock =[]
        inrockerr = []
        for k in bv.ksquad:
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
                              
            
            sig =  A2_quad_average / A1.sigma_coh_det( Nr, 0., 0.)
            inrock.append( sig.nominal_value) 
            inrockerr.append( sig.std_dev)
        np.savetxt( 'rockingdat/inrock'+fname, inrock)
        np.savetxt( 'rockingdat/inrockerr'+fname, inrockerr)

    inW,inWerr=plotRock(axinp, bv.ksquadth, inrock, inrockerr, normtof,'$N_{\mathrm{AFM}}$=%02d'%nafm, lc[i%nc], fc[i%nc])
    

    rockingW.append( [inW, inWerr])

rockingW = np.array(rockingW) 

# Plot the widths as a function of cristal size
axinpW = plt.subplot( gs[ 0,1] )

print nafm
print rockingW[:,0]
print rockingW[:,1]

axinpW.errorbar( nafms, rockingW[:,0], yerr=rockingW[:,1], \
           capsize=0., elinewidth = 1. ,\
           fmt='.', ecolor='black', mec='black', \
           mew=1., ms=5.,\
           marker='o', mfc='lightblue')

# Plot the widths that you expect from beam diffraction
# half 1/e^2 angle = lambda / (pi w0) 
# w0 = nafm/2 * sitespacing = nafm/2 * lambda/2
# half 1/e^2 angle = 4 / (pi nafm)
xdiff = np.linspace(min(nafms),max(nafms), 100)
Wdiff = (1/np.sqrt(2.))*1000. *   4. / (np.pi*xdiff)
axinpW.plot( xdiff, Wdiff, '-b',label=r'$4/(\sqrt{2}\pi N_{\mathrm{AFM}})$')
axinpW.legend(loc='best',numpoints=1,prop={'size':11},handlelength=1.1,handletextpad=0.5)

axinpW.set_xlabel('Nafm')
axinpW.set_ylabel('Rocking curve\n$1/e^{2}$ half width (mrad)',ha='center',labelpad=16)


allax = [axinp, axinpW ]
for ax in allax:
    ax.grid()

#axinp.set_ylim(0.,10)
#axout.set_ylim(0.,10)


axinp.legend(bbox_to_anchor=(1.0,1.0),loc='upper right',numpoints=1,prop={'size':5}, \
           handlelength=1.1,handletextpad=0.5)

gs.tight_layout(figure, rect=[0,0.0,1.0,0.96])
outfile = 'RockingSquadN%d.png' % N
figure.savefig(outfile, dpi=250)
        







