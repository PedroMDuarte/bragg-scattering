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

A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))

#A1.set_pbragg(0.0002)
#A2.set_pbragg(0.0002)

# Calculate Bragg scattering 100 for MIT experiment
# Q vector for bragg scattering
lRb = 780. 
l1064 = 1064.
Qmit = 4*np.pi/l1064 * vec3.vec3( 0, -1., 0) # 2pi/a
QUmit = Qmit/abs(Qmit)

ThetaMIT = np.arccos( abs(Qmit)/2. / ( 2.*np.pi / lRb ) ) 
print "Bragg angle for MIT 100 experiment = %.2f" % (ThetaMIT*180./np.pi)
kinMIT = vec3.vec3()
kinMIT.set_spherical( 2.*np.pi/lRb, np.pi/2., np.pi/2. - ThetaMIT)
koutMIT = kinMIT + Qmit

AMIT = afm.crystal(N, nafm, bv.l1064/2, (kinMIT,koutMIT))
AMIT.mass = 87.
AMIT.set_v0( [15.,15.,15.])



# Plot the decay of the debye-waller factor 
# as a function of time  
figure3 = plt.figure(figsize=(6.,4.5))
gs3 = matplotlib.gridspec.GridSpec( 1,1) 
figure3.suptitle('')
axDW = plt.subplot( gs3[0,0] )

time =  np.linspace( 0., 30., 200) 
A2DWtof = np.array([A2.dw_time(t)  for t in time ])
A1DWtof = np.array([A1.dw_time(t)  for t in time ])
mitDWtof = np.array([AMIT.dw_time(t) for t in time ])

axDW.plot( time, A2DWtof ,\
            '.-', mec='green', c='green',\
            mew=1.0, ms=3.,\
            marker='o', mfc='limegreen', \
            label='A2 Debye-Waller TOF, 20$E_{r}$') 

axDW.plot( time, A1DWtof ,\
            '.-', mec='blue', c='blue',\
            mew=1.0, ms=3.,\
            marker='o', mfc='lightblue', \
            label='A1 Debye-Waller TOF, 20$E_{r}$')

axDW.plot( time, mitDWtof ,\
            '.-', mec='black', c='black',\
            mew=1.0, ms=3.,\
            marker='o', mfc='gray', \
            label='MIT Debye-Waller TOF, 15$E_{r}$') 
axDW.grid()
axDW.set_xlabel('time of flight ($\mu$s)')

axDW.legend(loc='best',numpoints=1,\
         prop={'size':10}, \
         handlelength=1.1,handletextpad=0.5)
#axDWR.legend(loc='best',numpoints=1,\
#         prop={'size':10}, \
#         handlelength=1.1,handletextpad=0.5)
    
gs3.tight_layout(figure3, rect=[0,0.0,1.0,0.96])
outfile = 'DWtof_plot.png'
figure3.savefig(outfile, dpi=250)
     


# Plot the decay of A2/A1  
# as a function of time  
figure2 = plt.figure(figsize=(6.,4.5))
gs2 = matplotlib.gridspec.GridSpec( 1,1) 
figure2.suptitle('')
ax = plt.subplot( gs2[0,0] )

time =  np.linspace( 0., 6., 50)
Nr = 320 
normTOF = A2.sigma_coh_det(Nr,0.,100.) / A1.sigma_coh_det(Nr,0.,100.)
R_DWtof = [A2.sigma_coh_det(Nr,0.,t)/A1.sigma_coh_det(Nr,0.,t) / normTOF  for t in time ]

print normTOF

detRES = A1.d12 /2. 
normRES = A2.sigma_coh_det(Nr,detRES,0.) / A1.sigma_coh_det(Nr,detRES,0.) 
print normRES
normRES = normRES / normTOF
ax.axhspan( normRES.nominal_value - 0.03,\
             normRES.nominal_value + 0.03,\
             facecolor='gray', alpha=0.6, linewidth=0)

a2a1 = np.array( [s.nominal_value for s in R_DWtof])
a2a1err = np.array( [s.std_dev for s in R_DWtof] ) 

lcolor='blue'
fcolor='lightblue'
labelstr='N=40, nafm=8'
ax.errorbar(time, a2a1, yerr=a2a1err,\
            capsize=0., elinewidth = 1. ,\
            fmt='.-', ecolor=lcolor, mec=lcolor, \
            mew=1.0, ms=3.,\
            alpha = 1.0, \
            marker='o', mfc=fcolor, label=labelstr) 

ax.grid()
ax.set_xlabel('time of flight ($\mu$s)')
ax.set_ylabel('A2/A1')
ax.legend(loc='best',numpoints=1,\
         prop={'size':10}, \
         handlelength=1.1,handletextpad=0.5)
#ax.set_ylim(0.,2.)

gs2.tight_layout(figure2, rect=[0,0.0,1.0,0.96])
outfile = 'A2A1tof_plot.png'
figure2.savefig(outfile, dpi=250)


# Plot the decay of A2/A1  
# as a function of time 
# alongside with the data we have taken 
figure1 = plt.figure(figsize=(6.,4.5))
gs1 = matplotlib.gridspec.GridSpec( 1,1) 
figure1.suptitle('')
ax = plt.subplot( gs1[0,0] )


lcolor='blue'
fcolor='lightblue'
labelstr='N=40, nafm=8'
ax.plot(time, a2a1,'-',\
            alpha = 1.0, lw=2.5,\
            c=lcolor, label=labelstr)

# FETCH DATA 
tofdat={}
tofdat['130726_tof'] = { \
                  'label':'7$E_{r}$  July 26',\
                  'dir':'/lab/data/app3/2013/1307/130726/',\
                  'shots':'3019:3034,-3026,-3029,3039:3052,3055:3102',\
                  'ec':'blue','fc':'lightblue',\
                  }
tofdat['130725_tof'] = { \
                  'label':'7$E_{r}$  July 25',\
                  'dir':'/lab/data/app3/2013/1307/130725/',\
                  'shots':'2699:2758',\
                  'ec':'green','fc':'limegreen',\
                  }

tofdat['130801_tof'] = { \
                  'label':'7$E_{r}$  August 01',\
                  'dir':'/lab/data/app3/2013/1308/130801/',\
                  'shots':'4497:4546,-4539,4547:4586',\
                  'ec':'red','fc':'pink',\
                  }

datakeys = ['DIMPLELATTICE:tof', 'ANDOR1EIGEN:signal' , 'ANDOR2EIGEN:signal', 'HHHEIGEN:andor2norm', 'DIMPLELATTICE:force_lcr3' ]

print "Fetching data..."
for k in tofdat.keys():
    try:
        tofdat[k]['data'] = np.loadtxt(k+'.dat')
        print "Loaded %s succesfully." % k 
    except:
        print k 
        data, errmsg, rawdat = qrange.qrange_eval( tofdat[k]['dir'], tofdat[k]['shots'], datakeys) 
        np.savetxt(k+'.dat', data)
        tofdat[k]['data'] = data

print "Done."


base = 0.656
base_err = 0.068
ax.axhspan(  (base-base_err)/base, (base+base_err)/base, facecolor='gray', alpha=0.6, linewidth=0)

for k in sorted(tofdat.keys()):
    dat = tofdat[k]['data'][:,(0,3,4)]
    # Remove baselines from the data
    dat = dat[ dat[:,2] == -1. ] 
    dat = statdat.statdat( dat, 0, 1 )
    ax.errorbar( 1000.*dat[:,0], dat[:,1]/base, yerr=dat[:,3]/base,\
                  capsize=0., elinewidth=1.,\
                  fmt='.', ecolor=tofdat[k]['ec'], mec=tofdat[k]['ec'],\
                  mew=1.0, marker='o', mfc=tofdat[k]['fc'],\
                  label=tofdat[k]['label']) 
ax.grid()
ax.set_xlabel('Time of flight ($\mu$s)')
ax.set_ylabel('Bragg / Diffuse',ha='center',labelpad=20)
ax.set_xlim(-1.,25.)
#ax.xaxis.set_major_locator( matplotlib.ticker.MultipleLocator(0.1) )
#ax.set_ylim(0.92, 1.48)
ax.legend(loc='best',numpoints=1,prop={'size':8},\
           handlelength=1.1,handletextpad=0.5)

gs1.tight_layout(figure1, rect=[0,0.0,1.0,0.96])
outfile = 'A2A1tof_WithData_plot.png'
figure1.savefig(outfile, dpi=250)
