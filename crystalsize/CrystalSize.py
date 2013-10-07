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
sc=1.2
figure1 = plt.figure(figsize=(sc*8.,sc*7.))
gs1 = matplotlib.gridspec.GridSpec( 2,2) 
figure1.suptitle('')

figure2 = plt.figure(figsize=(sc*8.,sc*7.))
gs2 = matplotlib.gridspec.GridSpec( 2,2) 
figure2.suptitle('')

axR = figure2.add_subplot(gs2[0,0])
axA1 = figure2.add_subplot(gs2[0,1])
axA2 = figure2.add_subplot(gs2[1,1])

axR.set_ylabel('A2/A1 scaled such that\nA2/A1 at 0.5 ms == 1',ha='center',labelpad=20)
axA1.set_ylabel('A1 scaled such that\nA1 at 0.5 ms == 1',ha='center',labelpad=20)
axA2.set_ylabel('A2 scaled such that\nA2 at 0.5 ms == 1',ha='center',labelpad=20)

#axA1 = plt.subplot( gs2[0,0] )
#axA2 = plt.subplot( gs2[1,0] )
axSR = figure1.add_subplot( gs1[0,0] )
axS1 = figure1.add_subplot( gs1[0,1] )
axS2 = figure1.add_subplot( gs1[1,1] )
#axR = plt.subplot( gs2[1,2] )

axSR.set_ylabel('S2/N^3',ha='center',labelpad=20)
axS1.set_ylabel('S1',ha='center',labelpad=20)
axS2.set_ylabel('S2',ha='center',labelpad=20)

def plotCS( ax, pafms, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    ax.errorbar( pafms, dat[:,0], yerr=dat[:,1],\
                capsize=0., elinewidth = 1. ,\
                fmt='.-', c=lcolor,ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, label=labelstr) 

lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']

Ns = [ 5, 8, 10, 12, 20, 25, 30, 40, 50, 60]
#Ns = [  40, 50, 60]
for i,N in enumerate(Ns):
    #nafms = [2,4,6,7,8,10,12,16,20,24,32,34,38]
    nafms = [2,3,4,5,6,7,8,9,10,11,12]
    #nafms = [6,8,10,12]
    dirname='crystalsizedat/'
    fname = '_N%02d.dat'%N
    try: 
        S1_ = np.loadtxt(dirname+ 'S1'+fname)    
        S2_ = np.loadtxt(dirname+ 'S2'+fname)    
        A1_ = np.loadtxt(dirname+ 'A1'+fname)    
        A2_ = np.loadtxt(dirname+ 'A2'+fname)    
        RATIO = np.loadtxt(dirname+ 'R'+fname)    
        pafms = np.loadtxt(dirname+ 'pafms'+fname)    
    except:
        S1_ = []
        S2_ = []
        A1_ = []
        A2_ = []
        RATIO = []
        pafms = []
        for nafm in nafms:
            if nafm  > N:
                continue
            pafms.append(nafm)
            A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
            A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
        
            #Nr = 320
            Nr = 36
            A1sigma = A1.sigma_coh_det( Nr, 0., 0.) 
            A2sigma = A2.sigma_coh_det( Nr, 0., 0.)
      
            
            #detR = 76/5.9 /2. # On resonance
            detR = 0. 
            tof = 0.5
            A1sigmaNorm = A1.sigma_coh_det( Nr, detR, tof) 
            A2sigmaNorm = A2.sigma_coh_det( Nr, detR, tof) 
        
            norm = (A2sigma / A1sigma) / (A2sigmaNorm / A1sigmaNorm) 
            RATIO.append( [norm.nominal_value, norm.std_dev] )
            #print "Nafm = %d --> A2/A1= %.4g" % \
            #      (nafm, (A2sigma/A1sigma).nominal_value)
            #print "Nafm = %d --> A2/A1(tof)= %.4g" % \
            #      (nafm, (A2sigmaNorm/A1sigmaNorm).nominal_value)
        
            
            A1norm = A1sigma/ A1sigmaNorm
            A2norm = A2sigma/ A2sigmaNorm
            A1_.append( [A1norm.nominal_value, A1norm.std_dev] )
            A2_.append( [A2norm.nominal_value, A2norm.std_dev] )
    
            S1_.append( [np.mean(A1.Srandom), stats.sem(A1.Srandom)])
            S2_.append( [np.mean(A2.Srandom), stats.sem(A2.Srandom)])

        S1_ = np.array(S1_)
        S2_ = np.array(S2_)
        A1_ = np.array(A1_)
        A2_ = np.array(A2_)
        RATIO = np.array(RATIO)
        pafms = np.array(pafms)
        
        np.savetxt(dirname+ 'S1'+fname, S1_)    
        np.savetxt(dirname+ 'S2'+fname, S2_)    
        np.savetxt(dirname+ 'A1'+fname, A1_)    
        np.savetxt(dirname+ 'A2'+fname, A2_)    
        np.savetxt(dirname+ 'R'+fname, RATIO)    
        np.savetxt(dirname+ 'pafms'+fname, pafms)    
    
    
    Nstr = 'N=%d' % N
    plotCS( axR, pafms, RATIO, lc[i], fc[i], 'o', 'N=%d'%N)
    plotCS( axA1, pafms, A1_, lc[i], fc[i], 'o', Nstr)  
    plotCS( axA2, pafms, A2_, lc[i], fc[i], 'o', Nstr)

    plotCS( axSR, pafms, S2_/N**3, lc[i], fc[i], 'o', Nstr)
    plotCS( axS1, pafms, S1_, lc[i], fc[i], 'o', Nstr)  
    plotCS( axS2, pafms, S2_, lc[i], fc[i], 'o', Nstr)  

 

axes = [axR, axA1, axA2, axSR, axS1, axS2]
for ax in axes:
  ax.grid()
  ax.set_xlabel('Nafm (~correlation length)')
  ax.set_xlim(0.,12.5)

  ax.legend(loc='best',numpoints=1,\
           prop={'size':8}, \
           handlelength=1.1,handletextpad=0.5)

axA2.set_ylim(0.,6.)
axR.set_ylim(0.,6.)
axR.axhspan(  1.8,2.4, facecolor='gray', alpha=0.6, linewidth=0)


gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
outfile = 'CrystalSize_plot.png'
figure2.savefig(outfile, dpi=250)

gs1.tight_layout(figure1, rect=[0,0.0,1.0,0.96])
outfile = 'CrystalSize_plot_Sfactor.png'
figure1.savefig(outfile, dpi=250)



