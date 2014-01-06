
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
rows = 2

figure = plt.figure(figsize=(8.,6.))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42) 


lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','violet']
m  = [ 'o', 'D', '+', '^', '<', '>','o', 'D', '+', '^', '<', '>']

ax1e = plt.subplot( gs[ 1,0] )
ax2e = plt.subplot( gs[ 1,1] )
axRe = plt.subplot( gs[ 0,0] )



axlist = [ax1e, ax2e, axRe] 
for ax in axlist:
  ax.set_xlabel('Detuning ($\Gamma$)')
  # Put vertical lines to show the states
  ax.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)
  ax.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor='gray', alpha=0.6, linewidth=0)

ax1e.set_ylabel('ANDOR1')
ax2e.set_ylabel('ANDOR2')
axRe.set_ylabel('A2/A1')

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
    if ax == axRe:
        detdat = np.transpose( np.vstack(( x, unumpy.nominal_values(dat) )))
        np.savetxt( 'detuningSim.dat', detdat)
    #plotunc( axRcOD, x, a2_o/a1_o , lc[i], fc[i],\


#nafms = [4,6,8,10,12,16,20,24,32,34,38,40]
nafms = [4,5,6,7,8]

plotafmval = [4,5,6,7,8]
plotafmdat = []
for i,nafm in enumerate(nafms):
    N = 40 
    Nr = 200
    dirname = 'detuningsimdat/' 
    fname = '_N%02d_nafm%02d_Nr%04d.dat'%(N,nafm,Nr)
    x = np.hstack(( np.linspace(-50,-25,5), \
                    np.linspace(-20,20,21), \
                    np.linspace(25,50,5)))
    try:
        a1 = np.loadtxt( dirname + 'A1'+fname) 
        a2 = np.loadtxt( dirname + 'A2'+fname) 
    except:
        print "\nWorking on nafm = %d" % nafm 
        # Initialize crystal
        A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
        A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
    
    
        a1 = []
        a2 = []
        for det in x:
            print det, 
            sys.stdout.flush()
            a1e = A1.I_( Nr=Nr, detuning=det, ) / A1.I_( Nr=Nr, detuning=det, tof=100.)
            a2e = A2.I_( Nr=Nr, detuning=det, ) / A2.I_( Nr=Nr, detuning=det, tof=100.)
        
            a1.append( [a1e.nominal_value, a1e.std_dev])
            a2.append( [a2e.nominal_value, a2e.std_dev])
        
        a1 = np.array( a1) 
        a2 = np.array( a2) 
        np.savetxt( dirname + 'A1'+fname, a1) 
        np.savetxt( dirname + 'A2'+fname, a2) 
     
    # Make the array with uncertainty of the various cross sections
    # for the uncertainty use col=3 which is the standard error
    a1_e = unumpy.uarray( a1[:,0], a1[:,1] )  
    a2_e = unumpy.uarray( a2[:,0], a2[:,1] )  

    # PLOT ANDOR 1 
    plotunc( ax1e, x, a1_e, lc[i], 'None', m[i], 'A1 elastic' + ' Nafm=%d'%nafm)
 
    # PLOT ANDOR2
    plotunc( ax2e, x, a2_e, lc[i], 'None', m[i], 'A2 elastic' + ' Nafm=%d'%nafm)


    # PLOT ANDOR2/ANDOR1 ELASTIC RATIO
    plotunc( axRe, x, a2_e/a1_e , lc[i], fc[i],\
             'o', 'Nafm=%d, (x%.2g)'%(nafm,1./1.))

    if nafm in plotafmval:
        plotafmdat.append((x,a2_e/a1_e,i,nafm)) 


#ax1e.set_yscale('log')
#ax2e.set_yscale('log')

lims=50.
for ax in axlist:
  ax.grid()
  ax.set_xlim(-lims,lims)
  ax.legend(loc='best', numpoints=1, prop={'size':5})

Rlims=50.
axRe.set_xlim(-Rlims,Rlims)

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
    ax.errorbar( x, unumpy.nominal_values(dat),\
                yerr=0.*unumpy.std_devs(dat),\
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
