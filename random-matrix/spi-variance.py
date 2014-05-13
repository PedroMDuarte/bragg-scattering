
import sys
import numpy as np
import random

def sig( n = 40**3, p=0.5): 
    q=1.-p
    return np.sqrt( 16*p*q*( n*(2*p-1.)**2. + 2*p*q ) ) 

def mu( n = 40**3, p=0.5): 
    q=1.-p
    return n*(2*p-1.)**2. + 4*p*q 


import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})

figure = plt.figure(figsize=(10.,3.3))
gs = matplotlib.gridspec.GridSpec( 1,3, hspace=0.25, wspace=0.4,
                                   top=0.86,bottom=0.16,left=0.05,right=0.98) 

axMean = figure.add_subplot( gs[0,0] ) 
axSig = figure.add_subplot( gs[0,1] ) 
axRel = figure.add_subplot( gs[0,2] ) 

Ns = [ 20, 25, 30, 40]
for i,N in enumerate(Ns):
    p = np.linspace(0.5, 0.505, 100)
    lstr = '$N^{1/3}=%02d$'%N 
    axMean.plot( (p-0.5)*1000., mu( N**3, p ) , label=lstr)  
    axSig.plot( (p-0.5)*1000., sig( N**3, p ) , label=lstr)  
    axRel.plot( (p-0.5)*1000.,sig( N**3, p ) / mu(N**3,p), label=lstr )  

figure.suptitle('$S_{\mathbf{Q}}$',fontsize=16)
axMean.set_ylabel('$\mu$', rotation=0, labelpad=15, fontsize=14)
axSig.set_ylabel('$\delta$', rotation=0, labelpad=15, fontsize=14)
axRel.set_ylabel('$\delta/\mu$', rotation=0, labelpad=15, fontsize=14)

axMean.legend( bbox_to_anchor=(0.05,0.98), loc='upper left', numpoints=1,\
               prop={'size':9}, handlelength=1.1, handletextpad=0.5) 
axSig.legend( bbox_to_anchor=(0.05,0.98), loc='upper left', numpoints=1,\
               prop={'size':9}, handlelength=1.1, handletextpad=0.5) 
axRel.legend( bbox_to_anchor=(0.03,0.03), loc='lower left', numpoints=1,\
               prop={'size':9}, handlelength=1.1, handletextpad=0.5) 

for ax in [axMean, axSig, axRel]:
    ax.xaxis.set_major_locator( matplotlib.ticker.MaxNLocator(5) )
    ax.grid()
    ax.set_xlabel(r'$(p-0.5)\times 10^{3}$', fontsize=14)
    ax.set_xlim(0., 5.)
 
outfile = 'spi-variance.png'
figure.savefig(outfile, dpi=250)
