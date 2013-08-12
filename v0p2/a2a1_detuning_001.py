
import numpy as np
from scipy import stats
from statarray import statdat, bindat

#a2a1 = np.loadtxt('a2a1_130707_2300.dat')
#a2a1 = np.concatenate( (a2a1, np.loadtxt('a2a1_130708_1223.dat')), axis=0 )

#a2a1 = np.loadtxt('a2a1_130708_1654.dat')
#a2a1 = np.loadtxt('a2a1_130709_0030.dat')


import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc
rc('font',**{'family':'serif'})


# Data file
datfile = 'data001/a2a1_detuning_001.dat' 

# Values of nafm for which plots will be shown
#nafms = [4,6,8,10,12,16,20,24,32,34,38,40]
nafms = [12,24,38,40]

cols = 3
rows = 3

figure = plt.figure(figsize=(12.,12.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42) 

import fetchdata
from uncertainties import unumpy

lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
m  = [ 'o', 'D', '+', '^', '<', '>']

ax1e = plt.subplot( gs[ 0,0] )
ax2e = plt.subplot( gs[ 1,0] )
axRe = plt.subplot( gs[ 2,0] )

ax1i = plt.subplot( gs[ 0,1] )
ax2i = plt.subplot( gs[ 1,1] )
axR1ie = plt.subplot( gs[ 2,1] )


ax1c = plt.subplot( gs[ 0,2] )
#ax2c = plt.subplot( gs[ 1,2] )
#axRc = plt.subplot( gs[ 2,2] )

#axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie, ax1c, ax2c, axRc] 
axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie] 
for ax in axlist:
  ax.set_xlabel('Detuning ($\Gamma$)')
  # Put vertical lines to show the states
  ax.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)
  ax.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)

ax1e.set_ylabel('ANDOR1')
ax2e.set_ylabel('ANDOR2')
axRe.set_ylabel('A2(elastic)/A1(elastic)')
ax1i.set_ylabel('ANDOR1')
ax2i.set_ylabel('ANDOR2')
axR1ie.set_ylabel('A1(inelastic)/A1(elastic)')

def plotunc( ax, x, dat, lcolor, fcolor, marker, labelstr):
    print "plotting ", labelstr
    ax.errorbar( x, unumpy.nominal_values(dat),\
                yerr=unumpy.std_devs(dat),\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=5.,\
                alpha = 1.0, \
                marker=marker, mfc=fcolor, \
                label=labelstr) 

for i,nafm in enumerate(nafms):
    detuning = 6.44
    a1_elastic, a2_elastic = fetchdata.fetch_data_A1A2( {'afmsize':nafm, 'ai':0.}, 'det', datfile, ykey='sigma')
    a1_inelastic, a2_inelastic = fetchdata.fetch_data_A1A2( {'afmsize':nafm, 'ai':0.}, 'det', datfile, ykey='sigma_inelastic' )

    # Define the units for cross sections
    sunits = 9 * (671e-7**2) / 16 / ( np.pi**2)

 
    i % len(nafms)

    # Verify that the x array of the various cross sections are the same
    xlist = [a1_elastic[:,0], a2_elastic[:,0], a1_inelastic[:,0], a2_inelastic[:,0] ] 
    if not all( np.array_equal(x, xlist[0]) for x in xlist ):
        print "Error, cameras a2, a1 do not have the same detuning data" 
        exit() 
    x = xlist[0] 
    
    # Get the array with uncertainty of the various cross sections
    # for the uncertainty use col=3 which is the standard error
    a1_e = unumpy.uarray( sunits*a1_elastic[:,1], sunits*a1_elastic[:,3] )  
    a2_e = unumpy.uarray( sunits*a2_elastic[:,1], sunits*a2_elastic[:,3] ) 
    a1_i = unumpy.uarray( sunits*a1_inelastic[:,1], sunits*a1_inelastic[:,3] )  
    a2_i = unumpy.uarray( sunits*a2_inelastic[:,1], sunits*a2_inelastic[:,3] )

    #aR_e = Rscale * a1_e / a1_i
    #aR_i = Rscale * a1_i / a1_i 

    # PLOT ANDOR 1 
    plotunc( ax1e, x, a1_e, lc[i], 'None', m[i], 'A1 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax1i, x, a1_i, lc[i], 'None', m[i], 'A1 inelastic' + ' Nafm=%d'%nafm)
    print "max a1_i = %f" % unumpy.nominal_values(a1_i).max()
 
    # PLOT ANDOR2
    plotunc( ax2e, x, a2_e, lc[i], 'None', m[i], 'A2 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax2i, x, a2_i, lc[i], 'None', m[i], 'A2 inelastic' + ' Nafm=%d'%nafm)


    # PLOT ANDOR2/ANDOR1 ELASTIC RATIO
    rdat = np.transpose(np.vstack((x, unumpy.nominal_values(a2_e/a1_e))))
    rdat = bindat( rdat, 0, 1, -500., 500., 50)
    xdat = rdat[:,0]
    max0 = rdat[:,1].max() 
    rdat = unumpy.uarray( rdat[:,1], np.ones_like(rdat[:,1]))
    plotunc( axRe, xdat, rdat/max0, lc[i], fc[i],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./max0))

    # PLOT ANDOR2/ANDOR1 INELASTIC RATIO
    rdat = np.transpose(np.vstack((x, unumpy.nominal_values(a1_i/a1_e))))
    rdat = bindat( rdat, 0, 1, -500., 500., 50)
    xdat = rdat[:,0]
    max0 = rdat[:,1].max() 
    rdat = unumpy.uarray( rdat[:,1], np.ones_like(rdat[:,1]))
    plotunc( axR1ie, xdat, rdat/max0, lc[i], fc[i],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./max0))

    # PLOT ANDOR1 ELASTIC-INELATIC COMBINED  
    # The polarization sum and the Debye-Waller factor 
    # are always the same so they can be obtained from one element
    a1_polsum, a2_polsum = fetchdata.fetch_data_A1A2( {'afmsize':nafm, 'ai':0.}, 'det', datfile, ykey='polsum')
    a1_dw, a2_dw = fetchdata.fetch_data_A1A2( {'afmsize':nafm, 'ai':0.}, 'det', datfile, ykey='debye-waller')
    # Veryfy that indeed all entries for polsum and dw are equal
    if not (np.allclose( a1_polsum[:,1], a1_polsum[:,1][0]) and\
    np.allclose( a2_polsum[:,1], a2_polsum[:,1][0]) and\
    np.allclose( a1_dw[:,1], a1_dw[:,1][0]) and\
    np.allclose( a2_dw[:,1], a2_dw[:,1][0])):
        print "Error, polsum or debye-waller are NOT all equal for a given detuning"
    a1_polsum = a1_polsum[:,1][0]
    a2_polsum = a2_polsum[:,1][0]
    a1_dw = a1_dw[:,1][0]
    a2_dw = a2_dw[:,1][0]
   


ax1e.set_yscale('log')
ax2e.set_yscale('log')
ax1i.set_yscale('log')
ax2i.set_yscale('log')
#axR.set_yscale('log')

lims=100.
for ax in axlist:
  ax.set_xlim(-lims,lims)
  ax.legend(loc='best', numpoints=1, prop={'size':5})

Rlims=500.
axRe.set_xlim(-Rlims,Rlims)
Rlims=1000.
axR1ie.set_xlim(-Rlims,Rlims)

#ax1i.set_ylim(0.,0.000002)

gs.tight_layout(figure, rect=[0.,0.,1.0,0.91])
#plt.show()
figure.savefig('a2a1_detuning_001.png', dpi=200)
#pylab.clf()

