
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import numpy as np
import vec3
import pylab
from scipy import stats

import braggvectors as bv
import afm
 

from uncertainties import unumpy, ufloat
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})


cols = 2
rows = 3

figure = plt.figure(figsize=(10.,9.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.7, hspace=0.2, \
                                   top=0.97,bottom=0.12) 


lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
m  = [ 'o', 'D', '+', '^', '<', '>','o', 'D', '+', '^', '<', '>']

ax1e = plt.subplot( gs[ 0,0] )
ax2e = plt.subplot( gs[ 0,1] )
axMe = plt.subplot( gs[ 1,1] )
axRe = plt.subplot( gs[ 1,0] )
axMScaled = plt.subplot( gs[2,1])

unscaled = True





figure.text( 0.,0.015, r'*We define $\frac{X}{X_{\mathrm{TOF}}} = X_{t}$ for any quantity $X$.',\
             fontsize=8)
ax1e.set_ylabel(r'$A1_{t}-1$',fontsize=14,rotation=0)
ax2e.set_ylabel(r'$A2_{t}-1$',fontsize=14,rotation=0)
axMe.set_ylabel(r'$M_{t}-1$',fontsize=14,rotation=0)
axRe.set_ylabel(r'$(A2/A1)_{t}$',fontsize=14,rotation=0)
axMScaled.set_ylabel(r'$\frac{ M_{t}-1 }{ (M_{t}-1)_{P_{\mathrm{B}}=0}}$',\
                      fontsize=14,rotation=0)

def plotunc( ax, x, dat, lcolor, fcolor, marker, labelstr, alpha=1.0, errscale=1.0):
    print "plotting ", labelstr
    ax.errorbar( x, unumpy.nominal_values(dat),\
                yerr=errscale*unumpy.std_devs(dat),\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=5.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, \
                label=labelstr)
    if ax == axRe:
        pdat = np.transpose( np.vstack(( x, unumpy.nominal_values(dat) )))
        np.savetxt( 'pbraggSim.dat', pdat)
    #plotunc( axRcOD, x, a2_o/a1_o , lc[i], fc[i],\


#nafms = [4,5,6,7,8,10,12,16,20,24,32,34,38,40]
nafms = [4,5,6,7,8,9]

A1pfits = []
A2pfits = []
ARpfits = []

plotafmval = [4,5,6,7,8]
plotafmdat = []
for i,nafm in enumerate(nafms):
    N = 40 
    Nr = 10
    dirname = 'pbraggsimdat/' 
    fname = '_N%02d_nafm%02d_Nr%04d.dat'%(N,nafm,Nr)
    x = np.linspace( 1.0, 2000., 16 ) 
    x = np.hstack(( np.linspace(1.0, 500., 6), \
                    np.linspace(600,6000,12), \
                    np.linspace(8000,32000,4))  )
    try:
        a1 = np.loadtxt( dirname + 'A1'+fname) 
        a2 = np.loadtxt( dirname + 'A2'+fname) 
        m0 = np.loadtxt( dirname + 'MANTA'+fname) 
    except:
        print "\nWorking on nafm = %d" % nafm 
        # Initialize crystal
        A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA1))
        A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA2))
        M  = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kMANTA))
       
        det = -13.  
       
        a1 = []
        a2 = []
        m0 = [] 
        for PBragg in x:
            print PBragg, 
            sys.stdout.flush()
            a1e = A1.I_( Nr=Nr, pbragg=PBragg, detuning=det ) \
                / A1.I_( Nr=Nr, pbragg=PBragg, tof=100., detuning=det) 
            a2e = A2.I_( Nr=Nr, pbragg=PBragg, detuning=det ) \
                / A2.I_( Nr=Nr, pbragg=PBragg, tof=100., detuning=det)
            aMe = M.I_( Nr=Nr, pbragg=PBragg, detuning=det ) \
                / M.I_( Nr=Nr, pbragg=PBragg, tof=100., detuning=det)
        
            a1.append( [a1e.nominal_value, a1e.std_dev])
            a2.append( [a2e.nominal_value, a2e.std_dev])
            m0.append( [aMe.nominal_value, aMe.std_dev])
        
        a1 = np.array( a1) 
        a2 = np.array( a2)
        m0 = np.array( m0)  
        np.savetxt( dirname + 'A1'+fname, a1) 
        np.savetxt( dirname + 'A2'+fname, a2) 
        np.savetxt( dirname + 'MANTA'+fname, m0) 
     
    # Make the array with uncertainty of the various cross sections
    # for the uncertainty use col=3 which is the standard error
    a1_e = unumpy.uarray( a1[:,0], a1[:,1] )  
    a2_e = unumpy.uarray( a2[:,0], a2[:,1] )  
    aM_e = unumpy.uarray( m0[:,0], m0[:,1] )  


    # PLOT ANDOR 1 
    y  = (a1_e - 1)
    plotunc( ax1e, x, y , lc[i%10], 'None', m[i%10], 'A1' + ' Nafm=%d'%nafm)
 
    # PLOT ANDOR2
    y  = (a2_e - 1)
    plotunc( ax2e, x, y , lc[i%10], 'None', m[i%10], 'A2' + ' Nafm=%d'%nafm)

    # PLOT MANTA
    y  = (aM_e - 1)
    plotunc( axMe, x, y , lc[i%10], 'None', m[i%10], 'MANTA' + ' Nafm=%d'%nafm)

    # PLOT MANTA RESCALED 
    y  = (aM_e - 1)/ (aM_e[0]-1)
    plotunc( axMScaled, x, y , lc[i%10], 'None', m[i%10], 'A2' + ' Nafm=%d'%nafm, errscale=0.)
    pfit =  np.polyfit( x[0:6], unumpy.nominal_values(y[0:6]), 1) 
    z = np.poly1d( pfit )
    axMScaled.plot( x, z(x), '-', color=lc[i%10])


    # PLOT ANDOR2/ANDOR1 ELASTIC RATIO
    y = a2_e / a1_e
    plotunc( axRe, x, y , lc[i%10], fc[i%10],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./1.))

    if nafm in plotafmval:
        plotafmdat.append((x,a2_e/a1_e,i,nafm)) 

def strfit(f):
    out = ''
    for p in f:
       out = out + '%8.2g'%p
    return out 
        

print "A1 PFITS"
for f in A1pfits:
    print strfit(f) 
print "A2 PFITS"
for f in A2pfits:
    print strfit(f) 
print "AR PFITS"
for f in ARpfits:
    print strfit(f) 

#ax1e.set_yscale('log')
#ax2e.set_yscale('log')

lims=50.
axMScaled.set_ylim(0.,1.05)
axlist = [ax1e, ax2e, axRe, axMe, axMScaled] 
for ax in axlist:
  ax.set_xlabel('PBragg ($\mu\mathrm{W}$)')
  ax.set_xlim(-10.,6010)
  ax.grid()
  #ax.set_xlim(-lims,lims)
  #ax.legend(loc='best', numpoints=1, prop={'size':5})
  ax.legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':5}, \
           handlelength=1.1,handletextpad=0.5)

#gs.tight_layout(figure, rect=[0.,0.,1.0,0.91])
#plt.show()
if unscaled:
    outfile = 'pbragg_unscaled.png'
else:
    outfile = 'pbragg.png'
figure.savefig(outfile, dpi=200)
pylab.clf()


