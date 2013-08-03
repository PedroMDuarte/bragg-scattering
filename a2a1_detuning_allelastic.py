
import numpy as np
from scipy import stats
from statarray import statdat

#a2a1 = np.loadtxt('a2a1_130707_2300.dat')
#a2a1 = np.concatenate( (a2a1, np.loadtxt('a2a1_130708_1223.dat')), axis=0 )

#a2a1 = np.loadtxt('a2a1_130708_1654.dat')
#a2a1 = np.loadtxt('a2a1_130709_0030.dat')


import matplotlib.pyplot as plt
import matplotlib

from matplotlib import rc
rc('font',**{'family':'serif'})


# Data file
datfile = 'data001/a2a1_detuning_allelastic.dat' 

# Values of nafm for which plots will be shown
nafms = [4,6,8,10,12,16,20,24,32,34,38,40]

cols = 2
rows = len(nafms)/2+len(nafms)%2

figure = plt.figure(figsize=(10.8,3.6*rows))
#figure.suptitle('Bragg')
gs = matplotlib.gridspec.GridSpec( rows,cols, wspace=0.6, hspace=0.42) 

import fetchdata
from uncertainties import unumpy

for i,nafm in enumerate(nafms):
    detuning = 6.44
    a1, a2 = fetchdata.fetch_data_A1A2( {'afmsize':nafm, 'ai':0.}, 'det', datfile )

    # Put the units in the cross section
    sunits = 9 * (671e-7**2) / 16 / ( np.pi**2)
    a1[:,1] = sunits*a1[:,1]
    a1[:,2] = sunits*a1[:,2]
    a2[:,1] = sunits*a2[:,1]
    a2[:,2] = sunits*a2[:,2]
 
    i % len(nafms) 
    ax = plt.subplot( gs[ i%rows, i/rows] )
    ax.set_title('AFM = %d sites' % nafm)

    a1s = unumpy.uarray( a1[:,1] , a1[:,2] ) 
    a2s = unumpy.uarray( a2[:,1] , a2[:,2] )

    a2a1 = a2s/ a1s
  
    a2a1_mean = unumpy.nominal_values( a2a1 )
    a2a1_std  = unumpy.std_devs( a2a1)
  
   
    #ax.errorbar( a1[:,0], a1[:,1], yerr=a1[:,2], \
    #           capsize=0., elinewidth = 1. ,\
    #           fmt='.', ecolor='red', mec='red', \
    #           mew=1., ms=5.,\
    #           marker='o', mfc='pink', \
    #           label="A1")  

    #ax.errorbar( a2[:,0], a2[:,1], yerr=a2[:,2], \
    #           capsize=0., elinewidth = 1. ,\
    #           fmt='.', ecolor='green', mec='green', \
    #           mew=1., ms=5.,\
    #           marker='o', mfc='limegreen', \
    #           label="A2") 

    #ax2 = ax.twinx() 
    ax.errorbar( a2[:,0], a2a1_mean , yerr=a2a1_std, \
               capsize=0., elinewidth = 1. ,\
               fmt='.', ecolor='blue', mec='blue', \
               mew=1., ms=5.,\
               marker='o', mfc='lightblue', \
               label="A2/A1")  
    #ax2.set_ylabel('A2/A1')    
    ax.set_ylabel('A2/A1')    

    ax.grid()
    ax.set_xlabel('Detuning from state 2 ($\Gamma$)')
    #ax.set_ylabel('Cross section (cm$^{2}$)')

    if nafm == 40:
        ax.set_xlim(-10,10)

#plt.show()
figure.savefig('a2a1_detuning.png', dpi=140)
#pylab.clf()

