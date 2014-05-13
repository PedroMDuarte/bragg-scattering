
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
rc('font',**{'family':'serif'}

)

#Plot the intensity in A1, A2 as a function of N and Nafm 
sc=1.2
figure2 = plt.figure(figsize=(sc*8.,sc*6.))
gs2 = matplotlib.gridspec.GridSpec( 2,2, hspace=0.3, wspace=0.3,
                                   top=0.95,left=0.1,right=0.95) 
figure2.suptitle('')

axR = figure2.add_subplot(gs2[1,1])
axA1 = figure2.add_subplot(gs2[0,0])
axA2 = figure2.add_subplot(gs2[0,1])

axR.set_ylabel(r'$\frac{A2/A1}{(A2/A1)_{\mathrm{TOF}}}$',ha='center',\
               fontsize=20, labelpad=20)
axA1.set_ylabel(r'$A1/A1_{\mathrm{TOF}}$',ha='center',\
               fontsize=20, labelpad=20)
axA2.set_ylabel(r'$A2/A2_{\mathrm{TOF}}$',ha='center',\
               fontsize=20, labelpad=20)

def plotAVE( ax, pafms, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    print "plotting...  ", dat.shape
    x = pafms
    y = np.mean(dat[:,0])
    yerr = np.std( dat[:,0] )
    ax.errorbar( x, y, yerr=yerr,\
                #yerr=dat[:,1],\
                capsize=0., elinewidth = 1. ,\
                fmt='.-', c=lcolor,ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, label=labelstr) 

def plotCS( ax, pafms, dat, lcolor, fcolor, marker, labelstr, alpha=1.0):
    print "plotting...  ", dat.shape
    x = pafms*np.ones_like( dat[:,0] ) + 0.1*np.arange(len(dat[:,0])) 
    ax.errorbar( x, dat[:,0], \
                #yerr=dat[:,1],\
                capsize=0., elinewidth = 1. ,\
                fmt='.-', c=lcolor,ecolor=lcolor, mec=lcolor, \
                mew=1.0, ms=3.,\
                alpha = alpha, \
                marker=marker, mfc=fcolor, label=labelstr) 

lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']

Ns = [ 5, 8, 10, 12, 20, 25, 30, 40]
#Ns = [  40, 50, 60]
for i,N in enumerate(Ns):
    print "\n",N
    #nafms = [2,4,6,7,8,10,12,16,20,24,32,34,38]
    #nafms = [2,3,4,5,6,7,8,9,10,12,16,20,24,32,38]
    #nafms = [6,8,10,12]
    nafms = [2]
    dirname='noisedat/'
    fname = '_N%02d.dat'%N
    try: 
        A1_ = np.loadtxt(dirname+ 'NOFILE'+fname)    
        A1_ = np.loadtxt(dirname+ 'A1'+fname)    
        A2_ = np.loadtxt(dirname+ 'A2'+fname)    
         
    except:
        A1_ = []
        A2_ = []
        RATIO = []
        nafm = 1 
        print nafm,",",
        sys.stdout.flush()
        if nafm  > N:
            continue
        A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
        A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
        
        Nr = 120
        #Nr = 36
        A1sigma = A1.I_tof( Nr, 0., return_array = True) 
        A2sigma = A2.I_tof( Nr, 0., return_array = True)

        # TOF norm 
        A1sigmaNorm = A1.I_tof( Nr, 100., return_array = True) 
        A2sigmaNorm = A2.I_tof( Nr, 100., return_array = True) 
      

        A1_ = np.column_stack(( A1sigma/A1sigmaNorm, A1sigma, A1sigmaNorm ))
        A2_ = np.column_stack(( A2sigma/A2sigmaNorm, A2sigma, A2sigmaNorm )) 
        
        np.savetxt(dirname+ 'A1'+fname, A1_)    
        np.savetxt(dirname+ 'A2'+fname, A2_)   

         
    
    
    Nstr = '$N=%d$' % N
    plotAVE( axA1, N, A1_, lc[i], fc[i], 'o', Nstr)  
    plotAVE( axA2, N, A2_, lc[i], fc[i], 'o', Nstr)



axes = [axR, axA1, axA2]
for ax in axes:
  ax.grid()
  ax.set_xlabel('$L_{\mathrm{AFM}}$',fontsize=16)

  ax.legend(loc='best',numpoints=1,\
           prop={'size':8}, \
           handlelength=1.1,handletextpad=0.5)

#axA2.set_ylim(0.,6.)
#axR.set_ylim(0.,6.)

axR.set_xlim(0.,12)
axA2.set_xlim(0.,50.)
axA1.set_xlim(0.,50.)

#gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
outfile = 'Noise_plot.png'
figure2.savefig(outfile, dpi=250)




