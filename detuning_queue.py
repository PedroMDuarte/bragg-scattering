import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
sys.path.append('/lab/software/apparatus3/py')
import qrange, statdat, fitlibrary
from uncertainties import ufloat


# Setup data locations and fetch data

detdat={}
detdat['130729_detuning'] = { \
                  'label':'7$E_{r}$  July 29',\
                  'dir':'/lab/data/app3/2013/1307/130729/',\
                  'shots':'3243:3258,3259:3357,-3325,3364:3379',\
                  'ec':'blue','fc':'lightblue',\
                  }

datakeys = ['DIMPLELATTICE:imgdet', 'ANDOR1EIGEN:signal' , 'ANDOR2EIGEN:signal', 'HHHEIGEN:andor2norm', 'DIMPLELATTICE:force_lcr3' ]

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

print "Done."
        

 
# Get figure started

from matplotlib import rc
rc('font',**{'family':'serif'})
figure = plt.figure(figsize=(4.,3.))
gs = matplotlib.gridspec.GridSpec( 1,1) 
figure.suptitle('')
ax1 = plt.subplot( gs[0,0] )

base = 0.656
base_err = 0.068
base1 = base
ax1.axvspan(  37.9/5.9-0.5, 37.9/5.9+0.5, facecolor='lightblue', alpha=0.6, linewidth=0)
ax1.axvspan(  -37.9/5.9-0.5, -37.9/5.9+0.5, facecolor='lightblue', alpha=0.6, linewidth=0)
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
                  #label=detdat[k]['label']) 
                  label='$7E_{r}$\n$a_{s}=190a_{0}$') 

#ax1.errorbar( dat[:,0], dat[:,1]/zero, yerr=dat[:,3]/zero, fmt='.', ecolor=edgecolor, mec=edgecolor, mew=1.0, marker='o', mfc=fillcolor) 

ax1.grid()
ax1.set_xlabel('Detuning ($\Gamma$)')
ax1.set_ylabel('Bragg / Diffuse',ha='center',labelpad=20)
ax1.set_xlim(-15.,15.)
#ax1.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(0.1) )
#ax1.set_ylim(0.92, 1.48)
ax1.legend(loc='upper right',numpoints=1,prop={'size':9}, handlelength=1.1,handletextpad=0.5)
#ax1.legend(bbox_to_anchor=(0.85,0.3),loc='center',numpoints=1,prop={'size':11}, \
#           handlelength=1.1,handletextpad=0.5)


gs.tight_layout(figure, rect=[0,0.0,1.,0.91])
outfile = 'detuning_queue.png'
plt.savefig(outfile, dpi=250)

exit()

