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
rc('font',**{'family':'serif'}

)

#Plot the structure factors as a function of N and Nafm 
figure2 = plt.figure(figsize=(12.,9.))
gs2 = matplotlib.gridspec.GridSpec( 2,3) 
figure2.suptitle('')
axC1 = plt.subplot( gs2[0,0] )
axC2 = plt.subplot( gs2[1,0] )
axS1 = plt.subplot( gs2[0,1] )
axS2 = plt.subplot( gs2[1,1] )
axSR = plt.subplot( gs2[0,2] )
axOnafm = plt.subplot( gs2[1,2] )

def plotCS( ax, nafmx, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    ax.errorbar( nafmx, dat[:,0], yerr=dat[:,1],\
                capsize=0., elinewidth = 1. ,\
                fmt='.-', ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, label=labelstr) 

lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray']

Ns = [  20, 25, 30, 40, 50, 60]
for i,N in enumerate(Ns):
    #nafms = [2,4,6,7,8,10,12,16,20,24,32,34,38]
    nafms = [2,3,4,5,6,7,8,9,10,11,12]
    S1 = []
    S2 = []
    C1 = []
    C2 = []
    O = []
    pafms = []
    for nafm in nafms:
        if nafm + 3 >= N:
            continue
        pafms.append(nafm)
        A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
        A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
    
        Nr = 320
        A1sigma = A1.sigma_coh_det( Nr, 0., 0.) 
        A2sigma = A2.sigma_coh_det( Nr, 0., 0.) 
    
        detR = 76/5.9 /2. # On resonance
        A1sigmaR = A1.sigma_coh_det( Nr, detR, 0.) 
        A2sigmaR = A2.sigma_coh_det( Nr, detR, 0.) 
    
        norm = (A2sigma / A1sigma) / (A1sigmaR / A2sigmaR) 
        print "Nafm = %d --> A2/A1 norm = %.2f" % (nafm, norm.nominal_value)
        O.append( [norm.nominal_value, norm.std_dev] )
    
        C1.append( [np.mean(A1.Crandom), stats.sem(A1.Crandom)])
        C2.append( [np.mean(A2.Crandom), stats.sem(A2.Crandom)])
        S1.append( [np.mean(A1.Srandom), stats.sem(A1.Srandom)])
        S2.append( [np.mean(A2.Srandom), stats.sem(A2.Srandom)])
    
    
    S1 = np.array(S1)
    C1 = np.array(C1)
    S2 = np.array(S2)
    C2 = np.array(C2)
    O = np.array(O)
    
    Nstr = ', N=%d' % N
      
    nafmx = np.array(pafms)
    plotCS( axS1, nafmx, S1, lc[i], fc[i], 'o', 'S1'+Nstr)  
    plotCS( axC1, nafmx, C1, lc[i], fc[i], 'o', 'C1'+Nstr)  
    plotCS( axS2, nafmx, S2, lc[i], fc[i], 'o', 'S2'+Nstr)  
    plotCS( axC2, nafmx, C2, lc[i], fc[i], 'o', 'C2'+Nstr)
    plotCS( axOnafm, nafmx, O, lc[i], fc[i], 'o', 'Bragg/Diffuse'+Nstr)
    plotCS( axSR, nafmx, S2/S1, lc[i], fc[i], 'o', 'S2/S1'+Nstr)

 

axes = [axS1, axC1, axS2, axC2, axSR, axOnafm]
for ax in axes:
  ax.grid()
  ax.set_xlabel('Nafm')

  ax.legend(loc='best',numpoints=1,\
           prop={'size':12}, \
           handlelength=1.1,handletextpad=0.5)

axOnafm.set_ylim(0.,10)

gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
outfile = 'CrystalSize_plot.png'
figure2.savefig(outfile, dpi=250)




