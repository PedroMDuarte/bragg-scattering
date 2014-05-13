
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
s = 0.8
figure = plt.figure(figsize=(16.*s,12.*s))
gs = matplotlib.gridspec.GridSpec( 4,5, hspace=0.25, wspace=0.25,
                                   top=0.98,bottom=0.07,left=0.05,right=0.98) 

axL = [[ figure.add_subplot( gs[j,i] ) for i in range(5)] for j in range(4) ] 

def plotHIST( ax, x, dat, xmax, **kwargs ):
    n, bins, patches = ax.hist(dat[:,0], 75, (-1., xmax), normed=False, \
        alpha=0.75, **kwargs)
    mu = np.mean(dat[:,0]) 
    sig = np.std(dat[:,0])
    ax.axvline( mu, color='black' )  
    ax.axvspan( mu-sig, mu+sig,  facecolor='lightgray', alpha=0.4)
    textstr = '$N^{1/3}=%d$'%N+'\n'+'$L_{AFM}=%d$'%nafm + '\n' + \
               '$\mu=%.3g$'% mu + '\n' + '$\sigma=%.2g$'%sig + '\n'\
               + '$\sigma/\mu=%.3g$'%(sig/mu)
    ax.text( 0.95,0.95, textstr, \
        va='top', ha='right', fontsize=10,transform=ax.transAxes,\
        bbox = dict(boxstyle='round', facecolor='lightgrey', alpha=1.) ) 
     


lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', 'purple', 'gray','orange','green']

xmax = [ [9.]*4,  [12.]*4, [12.]*4, [80., 55., 40., 25.], [300., 190., 130., 75.] ] 


Ns = [ 20, 25, 30, 40]
for i,N in enumerate(Ns):
    print "\n",N
    nafms = [1,5,6,9,12]
    dirname='noisedat/'
    Nstr = '$N=%d$' % N
    for j,nafm in enumerate(nafms):
        Nr = 520
        fname = '_N%02d_%02d_%04d.dat'%(N,nafm,Nr)
        print nafm,",",
        sys.stdout.flush()
        if nafm  > N:
            continue
        try: 
            #A1_ = np.loadtxt(dirname+ 'NOFILE'+fname)    
            A1_ = np.loadtxt(dirname+ 'A1'+fname)    
            A2_ = np.loadtxt(dirname+ 'A2'+fname)   
        except:
        
            A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
            A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))
            
            # In-situ 
            A1sigma = A1.I_tof( Nr, 0., return_array = True) 
            A2sigma = A2.I_tof( Nr, 0., return_array = True)
            # TOF norm 
            A1sigmaNorm = A1.I_tof( Nr, 100., return_array = True) 
            A2sigmaNorm = A2.I_tof( Nr, 100., return_array = True) 
    

            A1_ = np.column_stack(( A1sigma/A1sigmaNorm, \
                                    A1sigma, A1sigmaNorm ))
            A2_ = np.column_stack(( A2sigma/A2sigmaNorm, \
                                    A2sigma, A2sigmaNorm )) 
            
            np.savetxt(dirname+ 'A1'+fname, A1_)    
            np.savetxt(dirname+ 'A2'+fname, A2_)   
   
        plotHIST( axL[i][j], N, A2_, xmax[j][i], color=lc[i], facecolor=fc[i], label=Nstr)


xlims = [(0., 5.), (0.,5.)] +  [(0., None)]*8

for i,axR in enumerate(axL):
    for j,ax in enumerate(axR):
        if j <= 1 :
            ax.set_xlim( xlims[i] ) 
        else:
            ax.set_xlim( 0., None) 
        ax.grid()
        if j==0: 
            ax.set_ylabel('counts/bin', fontsize=14)
        if i==len(axL)-1:
            ax.set_xlabel('$S_{\mathbf{\pi}}$',fontsize=18)
    
  #ax.set_xlabel('$L_{\mathrm{AFM}}$',fontsize=16)

  #ax.legend(loc='best',numpoints=1,\
  #         prop={'size':8}, \
  #         handlelength=1.1,handletextpad=0.5)

#axA2.set_ylim(0.,6.)
#axR.set_ylim(0.,6.)


outfile = 'Noise_hist.png'
figure.savefig(outfile, dpi=250)




