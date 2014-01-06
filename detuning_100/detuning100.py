
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import numpy as np
import vec3
import pylab
from scipy import stats

import braggvectors as bv
import afm
from uncertainties import unumpy, ufloat

import statdat

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})

N = 40 
nafm = 8

A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kA2))
MANTA = afm.crystal(N, nafm, bv.l1064/2, (bv.kin100,bv.kMANTA))


figure2 = plt.figure(figsize=(8.,7.))
gs2 = matplotlib.gridspec.GridSpec( 2,2) 
ax = plt.subplot( gs2[0,0] )
axMANTA = plt.subplot( gs2[0,1] )
axA1 = plt.subplot( gs2[1,0] )
axA2 = plt.subplot( gs2[1,1] )

detuningX =  np.linspace( -30., 30., 50)
Nr = 20
Er = 20. 
PBragg = 500.

R_det  = []
A1_det = []
A2_det = []
MANTA_det = []
for det in detuningX:
    print det
    m  = MANTA.I_( Nr=Nr, detuning=det, tof=0., pbragg=PBragg, v0=Er)
    a1 = A1.I_( Nr=Nr, detuning=det, tof=0., pbragg=PBragg, v0=Er)
    a2 = A2.I_( Nr=Nr, detuning=det, tof=0., pbragg=PBragg, v0=Er)
    
    R_det.append( m/a1 )
    A1_det.append( a1 ) 
    A2_det.append( a2 ) 
    MANTA_det.append( m ) 



a2a1 = np.array( [s.nominal_value for s in R_det])
a2a1err = np.array( [s.std_dev for s in R_det] ) 
np.savetxt('det100_%.1fPBragg.dat' % PBragg, np.transpose( np.vstack( (detuningX,a2a1))))


lcolor='blue'
fcolor='lightblue'
labelstr='N=40, nafm=8, %.1f$E_{r}$ M/A1' % Er
ax.errorbar(detuningX, a2a1, yerr=a2a1err,\
            capsize=0., elinewidth = 1. ,\
            fmt='.-', ecolor=lcolor, mec=lcolor, \
            mew=1.0, ms=3.,\
            alpha = 1.0, \
            marker='o', mfc=fcolor, label=labelstr)

labelstr='N=40, nafm=8, %.1f$E_{r}$ M' % Er
y = np.array( [val.nominal_value for val in MANTA_det ] ) 
np.savetxt('det100_MANTA_%.1fEr.dat' % Er, np.transpose( np.vstack( (detuningX,y))))
axMANTA.plot(detuningX, y,\
            '.-', color=lcolor, mfc=fcolor, \
            mew=1.0, ms=3.,\
            alpha = 1.0, \
            marker='o', label=labelstr)

labelstr='N=40, nafm=8, %.1f$E_{r}$ A1' % Er
y = np.array( [val.nominal_value for val in A1_det ] ) 
np.savetxt('det100_A1_%.1fEr.dat' % Er, np.transpose( np.vstack( (detuningX,y))))
axA1.plot(detuningX, y,\
            '.-', color=lcolor, mfc=fcolor, \
            mew=1.0, ms=3.,\
            alpha = 1.0, \
            marker='o', label=labelstr)

labelstr='N=40, nafm=8, %.1f$E_{r}$ A2' % Er
y = np.array( [val.nominal_value for val in A2_det ] ) 
np.savetxt('det100_A2_%.1fEr.dat' % Er, np.transpose( np.vstack( (detuningX,y))))
axA2.plot(detuningX, y,\
            '.-', color=lcolor, mfc=fcolor, \
            mew=1.0, ms=3.,\
            alpha = 1.0, \
            marker='o', label=labelstr)

ax.set_ylabel(r'$\frac{S_{\mathrm{B}}}{S_{\mathrm{db}}}$',rotation=0,fontsize=24,labelpad=10)
axMANTA.set_ylabel('MANTA',labelpad=10,fontsize=14)	
axA1.set_ylabel('ANDOR1',labelpad=10,fontsize=14)
axA2.set_ylabel('ANDOR2',labelpad=10,fontsize=14)
titlestr='N=%d,  nafm=%d,  %.1f$E_{r}$,  %.1f uW' % (N,nafm,Er, PBragg)
figure2.suptitle(titlestr)


for axs in [ax, axMANTA, axA1, axA2]:
    axs.grid()
    axs.set_xlabel('detuning($\Gamma$)')
    #axs.legend(loc='best',numpoints=1,\
    #         prop={'size':10}, \
    #         handlelength=1.1,handletextpad=0.5)
    axs.set_xlim(-30.,30.)

gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
outfile = 'detuning_plot_%.1fEr_%.1fuW.png' % (Er,PBragg)
figure2.savefig(outfile, dpi=250)


