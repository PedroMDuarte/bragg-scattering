
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
Er = 20


# Plot the decay of A2/A1  
# as a function of time  
A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
A1.set_v0([Er, Er, Er])
A2.set_v0([Er, Er, Er])

time =  np.linspace( 0., 100., 100)
Nr = 320

#figure2 = plt.figure(figsize=(8.,7.))
#gs2 = matplotlib.gridspec.GridSpec( 2,2) 
#figure2.suptitle('')
#ax = plt.subplot( gs2[0,0] )
#axA1 = plt.subplot( gs2[1,0] )
#axA2 = plt.subplot( gs2[1,1] )
#
#normTOF1 = A1.I_tof(Nr, 100.)
#normTOF2 = A2.I_tof(Nr, 100.)
#normTOF  = normTOF2 / normTOF1 
#
#R_DWtof=[]
#A1_DWtof=[]
#A2_DWtof=[]
#for t in time:
#    if t%10==0: print t 
#    a1t = A1.I_tof(Nr, t)
#    a2t = A2.I_tof(Nr, t)
#    R_DWtof.append( a2t/a1t ) 
#    A1_DWtof.append( a1t )
#    A2_DWtof.append( a2t )
#
#print
#print normTOF
#
#detRES = A1.d12 /2. 
#normRES = A2.sigma_coh_det(Nr,detRES,0.) / A1.sigma_coh_det(Nr,detRES,0.) 
#print normRES
#normRES = normRES / normTOF
#ax.axhspan( normRES.nominal_value - 0.03,\
#             normRES.nominal_value + 0.03,\
#             facecolor='gray', alpha=0.6, linewidth=0)
#
#
#a2a1 = np.array( [s.nominal_value for s in R_DWtof])
#a2a1err = np.array( [s.std_dev for s in R_DWtof] ) 
#np.savetxt('DW_%.1fEr.dat' % Er, np.transpose( np.vstack( (time,a2a1))))
#
#lcolor='blue'
#fcolor='lightblue'
#labelstr='N=40, nafm=8, %.1f$E_{r}$' % Er
#ax.errorbar(time, a2a1, yerr=a2a1err,\
#            capsize=0., elinewidth = 1. ,\
#            fmt='.-', ecolor=lcolor, mec=lcolor, \
#            mew=1.0, ms=3.,\
#            alpha = 1.0, \
#            marker='o', mfc=fcolor, label=labelstr) 
#
#ax.grid()
#ax.set_xlabel('time of flight ($\mu$s)')
#labelstr='N=40, nafm=8, %.1f$E_{r}$ A1' % Er
#y = np.array( [val.nominal_value for val in A1_DWtof ] ) 
#np.savetxt('DW100_A1_%.1fEr.dat' % Er, np.transpose( np.vstack( (time,y))))
#axA1.plot(time, y,\
#            '.-', color=lcolor, mfc=fcolor, \
#            mew=1.0, ms=3.,\
#            alpha = 1.0, \
#            marker='o', label=labelstr)
#
#labelstr='N=40, nafm=8, %.1f$E_{r}$ A2' % Er
#y = np.array( [val.nominal_value for val in A2_DWtof ] ) 
#np.savetxt('DW100_A2_%.1fEr.dat' % Er, np.transpose( np.vstack( (time,y))))
#axA2.plot(time, y,\
#            '.-', color=lcolor, mfc=fcolor, \
#            mew=1.0, ms=3.,\
#            alpha = 1.0, \
#            marker='o', label=labelstr)
#ax.set_ylabel(r'$\frac{S_{\mathrm{B}}}{S_{\mathrm{db}}}$',rotation=0,fontsize=24,labelpad=10)
#
#
#axA1.set_ylabel('ANDOR1',labelpad=10,fontsize=14)
#axA2.set_ylabel('ANDOR2',labelpad=10,fontsize=14)
#titlestr='N=%d,  nafm=%d,  %.1f$E_{r}$' % (N,nafm,Er)
#for axs in [ax, axA1, axA2]:
#    axs.grid()
#    axs.set_xlabel('time of flight ($\mu$s)')
#    #axs.legend(loc='best',numpoints=1,\
#    #         prop={'size':10}, \
#    #         handlelength=1.1,handletextpad=0.5)
#    axs.set_xlim(0.,10.)
#
#gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
#outfile = 'tof_plot_%.1fEr.png' % Er
#figure2.savefig(outfile, dpi=250)
#
#
