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


cols = 3
rows = 4

figure = plt.figure(figsize=(12.,12.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42) 


lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
m  = [ 'o', 'D', '+', '^', '<', '>','o', 'D', '+', '^', '<', '>']

ax1e = plt.subplot( gs[ 0,0] )
ax2e = plt.subplot( gs[ 1,0] )
axRe = plt.subplot( gs[ 2,0] )

ax1i = plt.subplot( gs[ 0,1] )
ax2i = plt.subplot( gs[ 1,1] )
axR1ie = plt.subplot( gs[ 2,1] )


ax1c = plt.subplot( gs[ 0,2] )
ax2c = plt.subplot( gs[ 1,2] )
axRc = plt.subplot( gs[ 2,2] )

axRcOD = plt.subplot(gs[3,2])

#axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie, ax1c, ax2c, axRc] 
axlist = [ax1e, ax2e, axRe, ax1i, ax2i, axR1ie, ax1c, ax2c, axRc, axRcOD] 
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

def plotunc( ax, x, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    print "plotting ", labelstr
    ax.errorbar( x, unumpy.nominal_values(dat),\
                yerr=unumpy.std_devs(dat),\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=5.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, \
                label=labelstr)
    if ax == axRcOD:
        detdat = np.transpose( np.vstack(( x, unumpy.nominal_values(dat) )))
        np.savetxt( 'detuningSim.dat', detdat)
    #plotunc( axRcOD, x, a2_o/a1_o , lc[i], fc[i],\


#nafms = [4,6,8,10,12,16,20,24,32,34,38,40]
nafms = [2,4,6,7,8,9,10,16]

plotafmval = [6,7,8,9,10]
plotafmdat = []
for i,nafm in enumerate(nafms):
   
    print "\nWorking on nafm = %d" % nafm 
    # Initialize crystal
    N = 40 
    A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
    A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))

    Nr = 320
    A1.loadCS(Nr)
    A2.loadCS(Nr)


    a1 = []
    a2 = []
    x = np.hstack(( np.linspace(-400,-20,50), \
                    np.linspace(-20,20,101), \
                    np.linspace(20,400,50)))
    for det in x:
        a1i = A1.sigma_incoh_det( Nr, det, 0.)
        a2i = A2.sigma_incoh_det( Nr, det, 0.)
        a1e = A1.sigma_coh_det( Nr, det, 0.)
        a2e = A2.sigma_coh_det( Nr, det, 0.)
    
        a1.append( [a1e.nominal_value, a1e.std_dev, a1i.nominal_value,\
                    a1e.nominal_value, a1e.std_dev, \
                    a1e.nominal_value, a1e.std_dev, ])
        a2.append( [a2e.nominal_value, a2e.std_dev, a2i.nominal_value,\
                    a2e.nominal_value, a2e.std_dev, \
                    a2e.nominal_value, a2e.std_dev, ])
    
    a1 = np.array( a1) 
    a2 = np.array( a2) 
     
    # Make the array with uncertainty of the various cross sections
    # for the uncertainty use col=3 which is the standard error
    a1_e = unumpy.uarray( a1[:,0], a1[:,1] )  
    a1_i = unumpy.uarray( a1[:,2], np.zeros_like(a1[:,0]) )
    a1_c = unumpy.uarray( a1[:,3], a1[:,4] )  
    a1_o = unumpy.uarray( a1[:,5], a1[:,6] )  
    a2_e = unumpy.uarray( a2[:,0], a2[:,1] )  
    a2_i = unumpy.uarray( a2[:,2], np.zeros_like(a2[:,0]) )
    a2_c = unumpy.uarray( a2[:,3], a2[:,4] )  
    a2_o = unumpy.uarray( a2[:,5], a2[:,6] )  

    # PLOT ANDOR 1 
    plotunc( ax1e, x, a1_e, lc[i], 'None', m[i], 'A1 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax1i, x, a1_i, lc[i], 'None', m[i], 'A1 inelastic' + ' Nafm=%d'%nafm)
    plotunc( ax1c, x, a1_c, lc[i], 'None', m[i], 'A1 combined' + ' Nafm=%d'%nafm)
    print "max a1_i = %f" % unumpy.nominal_values(a1_i).max()
 
    # PLOT ANDOR2
    plotunc( ax2e, x, a2_e, lc[i], 'None', m[i], 'A2 elastic' + ' Nafm=%d'%nafm)
    plotunc( ax2i, x, a2_i, lc[i], 'None', m[i], 'A2 inelastic' + ' Nafm=%d'%nafm)
    plotunc( ax2c, x, a2_c, lc[i], 'None', m[i], 'A2 combined' + ' Nafm=%d'%nafm)


    # PLOT ANDOR2/ANDOR1 ELASTIC RATIO
    max0 = unumpy.nominal_values(a2_e/a1_e).max()
    plotunc( axRe, x, a2_e/a1_e/max0 , lc[i], fc[i],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./max0))

    # PLOT ANDOR1 INELATIC / ANDOR1 ELASTIC
    max0 = unumpy.nominal_values(a1_i/a1_e).max()
    plotunc( axR1ie, x, a1_i/a1_e, lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)

    # PLOT ANDOR2/ANDOR1 COMBINED RATIO
    plotunc( axRc, x, a2_c/a1_c , lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)

    # PLOT ANDOR2/ANDOR1 COMBINED PLUS OD RATIO
    plotunc( axRcOD, x, a2_o/a1_o , lc[i], fc[i],\
             'o', 'Nafm=%d'%nafm)
    if nafm in plotafmval:
        plotafmdat.append((x,a2_o/a1_o,i,nafm)) 

 


ax1e.set_yscale('log')
ax2e.set_yscale('log')
ax1i.set_yscale('log')
ax2i.set_yscale('log')
ax1c.set_yscale('log')
ax2c.set_yscale('log')
#axR.set_yscale('log')

lims=100.
for ax in axlist:
  ax.set_xlim(-lims,lims)
  ax.legend(loc='best', numpoints=1, prop={'size':5})

Rlims=500.
axRe.set_xlim(-Rlims,Rlims)
Rlims=1000.
axR1ie.set_xlim(-Rlims,Rlims)
Rlims=40.
axRc.set_xlim(-Rlims,Rlims)
axRcOD.set_xlim(-Rlims,Rlims)

gs.tight_layout(figure, rect=[0.,0.,1.0,0.91])
#plt.show()
figure.savefig('detuning.png', dpi=200)
pylab.clf()



########### PLOT SIMULATION ALONGSIDE DATA
import statdat
#Fetch data
detdat={}
detdat['130729_detuning'] = { \
                  'label':'7$E_{r}$  July 29',\
                  'dir':'/lab/data/app3/2013/1307/130729/',\
                  'shots':'3243:3258,3259:3357,-3325,3364:3379',\
                  'ec':'black','fc':'black',\
                  }

datakeys = ['DIMPLELATTICE:imgdet', 'ANDOR1EIGEN:signal' , 'ANDOR2EIGEN:signal', 'HHHEIGEN:andor2norm', 'DIMPLELATTICE:force_lcr3' ]

print
print "STARTING DATA + SIMULATION FIGURE"
print "Fetching data..."
for k in detdat.keys():
    try:
        detdat[k]['data'] = np.loadtxt(k+'.dat')
        print "Loaded %s succesfully." % k 
    except:
        print k 
        data, errmsg, rawdat = qrange.qrange_eval( detdat[k]['dir'], detdat[k]['shots'], datakeys) 
        np.savetxt(k+'.dat', data)
        detdat[k]['data'] = data

print "Done.\n"
        

########### PREPARE FIGURE AND AXIS FOR DATA PLOT
figure1 = plt.figure(figsize=(4.8,3.6))
gs1 = matplotlib.gridspec.GridSpec( 1,1) 
figure1.suptitle('')
ax1 = plt.subplot( gs1[0,0] )

########### PLOT THE SIMULATION RESULTS
def plotsim( ax, x, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    #Find the value at resonance
    res =  abs(abs(x) - 6.44) < 0.1
    datres = dat[res]
    ratiores = np.mean( unumpy.nominal_values(dat[res]))
    normTOF = 0.4/0.66
    normTOF = 1.0
    ax.errorbar( x, normTOF*unumpy.nominal_values(dat)/ratiores,\
                yerr=0.*normTOF*unumpy.std_devs(dat)/ratiores,\
                capsize=0., elinewidth = 1. ,\
                fmt='.', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, \
                label=labelstr) 
for d in plotafmdat:
    plotsim( ax1, d[0], d[1], lc[d[2]], fc[d[2]],\
             'o', 'Nafm=%d'%d[3], alpha=0.5)

# Base is 0.5 TOF
base = 0.656
base_err = 0.068
base1 = base
## Base is on resonance
base = 0.45
base_err = 0.03
base1 = base

reslinec = 'limegreen'
ax1.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor=reslinec, alpha=0.6, linewidth=0)
ax1.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor=reslinec, alpha=0.6, linewidth=0)
ax1.axhspan(  (base-base_err)/base1, (base+base_err)/base1, facecolor='gray', alpha=0.6, linewidth=0)

for k in sorted(detdat.keys()):
    dat = detdat[k]['data'][:,(0,3,4)]
    
    # Remove baslines from the data
    dat = dat[ dat[:,2] == -1. ] 
    

    dat = statdat.statdat( dat, 0, 1 )
    ax1.errorbar( (dat[:,0]+131.)/5.9, dat[:,1]/base1, yerr=dat[:,3]/base1,\
                  capsize=0., elinewidth=1.,\
                  fmt='.', ecolor=detdat[k]['ec'], mec=detdat[k]['ec'],\
                  mew=1.0, marker='o', mfc=detdat[k]['fc'],\
                  alpha = 1.0,\
                  #label=detdat[k]['label']) 
                  label='$7E_{r}$, $a_{s}=190a_{0}$') 

#ax1.errorbar( dat[:,0], dat[:,1]/zero, yerr=dat[:,3]/zero, fmt='.', ecolor=edgecolor, mec=edgecolor, mew=1.0, marker='o', mfc=fillcolor) 


ax1.grid()
ax1.set_xlabel('Detuning ($\Gamma$)')
ax1.set_ylabel('Bragg / Diffuse',ha='center',labelpad=20)
ax1.set_xlim(-15.,15.)
#ax1.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(0.1) )
#ax1.set_ylim(0.92, 1.48)

#ax1.legend(loc='upper left',numpoints=1,\
#           prop={'size':7}, \
#           handlelength=1.1,handletextpad=0.5)

ax1.legend(bbox_to_anchor=(1.0,1.0),loc='upper left',numpoints=1,prop={'size':7}, \
           handlelength=1.1,handletextpad=0.5)


gs1.tight_layout(figure1, rect=[0,0.0,0.82,0.96])
outfile = 'detuning_queue.png'
print outfile
figure1.savefig(outfile, dpi=250)




########### PLOT INSIDE/OUTSIDE RATIO AS A FUNCTION OF Nafm
#figure1 = plt.figure(figsize=(4.8,3.6))
#gs1 = matplotlib.gridspec.GridSpec( 1,1) 
#figure1.suptitle('')
#ax1 = plt.subplot( gs1[0,0] )
