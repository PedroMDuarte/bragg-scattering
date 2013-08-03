#
# Calculate Bragg scattering off of atoms in optical lattice
# Ted Corcovilos 20090129
# This version 20090226
# (using scipy module and iPython interpreter)
# all parameters in SI units unless specified
#
# Code assumes atoms are point scatterers and simply adds up the complex phases
# When changing array sizes, remember that coords are at _center_ of cell/pixel
# Changes:
# Version 7 - adding polarization
# Version 6 - added Debye-Waller factor, general cleaning
# Version 5 - corrected atomic structure factor (no exp)
#
# import libraries scipy and matplotlib
import sys, math
import numpy as np
#from W3J import * #code for Wigner 3j symbols
from scipy import * 
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# debug?
debug=0 #1 to skip loop, 0 to do it
# constants
comment="Version 7\n" #comment for csv file
filename="fm001.csv" #output filename
state=0 # 0 for ferromagnet, 1 for afm, 2 for random
plane=0 # 0 for (001), 1 for (1/2 1/2 1/2)
edge=41 # number of lattice sites along one edge of cube (needs to be odd)
# universal constants
c = 2.9979e8 # speed of light m/s
h = 6.626e-34 # planck constant Js
hbar = h/2./pi
kb = 1.3807e-23 # boltzmann constant J/K
amu = 1.66053886e-27 # atomic mass unit in kg
# Li constants
gamma = 5.92e6 # linewidth of transition in Hz (same for both states?)
Isat = 25.6 # saturation intensity of transition W/m^2 (using Foot's convention)
split = 76e6 # zeeman splitting of two magnetic sublevels in Hz (at Feshbach res.)
mass = 6 * amu # mass of Li-6
#define lattice
latticespacing = 532e-9
recoil = h*h/8/latticespacing/latticespacing/mass # recoil energy in J 
latticedepth = 10e-6 * kb # lattice depth in J, conveniant value
shofreq = sqrt(2.*latticedepth*pi**2/mass/latticespacing**2) #harmonic oscillator frequency of one site
sholength = sqrt(hbar/mass/shofreq) #harmonic oscillator length
# setup lattice coordinates
intmax = (edge-1)/2; intmin = -1*intmax
intcoords=mgrid[intmin:intmax:edge*1j,intmin:intmax:edge*1j,intmin:intmax:edge*1j] # cubic lattice centered at origin
#
#choose which sites are occupied
# define lattice as cube, for simplicity
#
#
# latticeoccupation = 0 for empty
#                     1 for state 1
#                     2 for state 2
if state == 0 :
	latticeoccupation = ones(intcoords.shape[1:4], dtype=np.int) #ferromagnetic case, all atoms in state 1
	latticestring = "polarized case\n" #label for output file
elif state == 1 :
	latticeoccupation = asarray(1 + (intcoords[0]+intcoords[1]+intcoords[2])%2, dtype=np.int) #antiferromagnetic case (alternate 1 and 2)
	latticestring = "AFM case\n" # label for output file
else:
	latticeoccupation = np.random.randint(1, 3, intcoords.shape[1:4]) #paramagnet, randomly choose state 1 or 2 for each site
	latticestring = "random case\n" # label for output file
#end if
coords = intcoords*latticespacing # scale to physical units
#
# define incoming probe beam
# From geometry
#probepower = 5e-3
#probebeamwaist = 100e-6 
#probeintensity = probepower/pi/probebeamwaist**2
# In terms of I sat
probeintensity=Isat
rabi = gamma*sqrt(probeintensity/2./Isat) # Rabi frequency in Hz
wavelength = 671e-9 #probe wavelength
detlow = split/2. # probe detuning from lower transition
detuning = array([0, detlow, detlow-split]) # indexed for "latticeoccupation"
#
# need to add stuff about probe polarization
# scattering cross section sigma will depend on
# polarization
#reset atten for debug
#atten = ones(atten.shape, dtype='complex')
#atten = asarray([0., 1., -1.])
kmag=2.*pi/wavelength #magnitude of k vectors
#define unit vectors of incoming beam:
if plane == 0 : # (0 0 1) case
	kinhat=array([sqrt(1.-(wavelength/2./latticespacing)**2),0.,-1*wavelength/2./latticespacing]) #unit vector along probe (looking for (100) peak)
else: # (1/2 1/2 1/2) case
	alpha = arcsin(sqrt(3)/4*wavelength/latticespacing) #half angle between kin and kout
	thetaQ = arctan(sqrt(2)) #polar angle of kout-kin
	thetain = pi/2. - thetaQ + alpha
	phiin = 1.25*pi
	kinhat= asarray([-1.*sin(thetain)/sqrt(2),-1.*sin(thetain)/sqrt(2),cos(thetain)])
#end if
kin = kmag * kinhat
#
# polarizations
# for this kinhat, find polarizations of interacting and spectator light for the two transitions
sigmaplus=sqrt(0.5)*array([1+0.j,0.+1.j,0.+0.j])
sigmaminus=sqrt(0.5)*array([1+0.j,0.-1.j,0.+0.j])
specsigp=cross(kinhat,sigmaplus)/sqrt(abs(dot(cross(kinhat,sigmaplus),cross(kinhat,sigmaplus))))
specsigm=cross(kinhat, sigmaminus)/sqrt(abs(dot(cross(kinhat, sigmaminus),cross(kinhat, sigmaminus))))
intsigp=cross(specsigp, kinhat)
intsigm=cross(specsigm, kinhat)
kincrossz=cross(kinhat, asarray([0,0,1]))
polin=kincrossz # Incoming polarization.  This is an arbitrary pick for now.
sigmamax = 3.*wavelength**2/2/pi # maximal scattering cross section for j=1/2 to j=3/2 transitions
sigma=sigmamax*asarray([0,abs(dot(polin,sigmaplus))**2,abs(dot(polin,sigmaminus))**2])
lineshape = rabi*(detuning + 0.5j*gamma)/(detuning*detuning+rabi**2/2+gamma**2/4)# !!Double check prefactor.  rabi or gamma?
atten = lineshape #folded cross-section into here also.  I guess atten is a bad name

#
# define ccd pixel locations
# Use rectangular array centered at kouthat
exposuretime = 1. #exposure in seconds (used to calculate photon counts...)
ccdr = 0.5 #distance from origin to center of ccd in meters
if plane == 0 : #(1 0 0) case
	ccdphi = 0 #azimuthal angle of ccd center
	ccdtheta = arccos(wavelength/2./latticespacing) #polar angle of ccd
else: # (1/2 1/2 1/2) case
	ccdtheta=pi/2. - thetaQ - alpha #see defn of thetaQ and alpha above
	ccdphi = 1.25*pi
#end if
# Define outgoing k vector
kouthat = array([cos(ccdphi)*sin(ccdtheta), sin(ccdphi)*sin(ccdtheta), cos(ccdtheta)]) #unit vector from origin to center of ccd
kout = kouthat*kmag
deltak = kout-kin
DW = exp(-1.*dot(deltak,deltak)*sholength**2/6) #Debye-Waller factor, assuming harmonic oscillator approximation
Q = deltak*latticespacing/2/pi # scattering momentum Miller indices
#
# Define the ccd
ccdphysicalsize=array([0.03,0.03]) #physical size of ccd array
ccdpixels = array([25,25]) #number of pixels on ccd
ccdpixelsize = ccdphysicalsize/ccdpixels
ccdpixelarea = ccdpixelsize[0]*ccdpixelsize[1]
ccdcenter = ccdr*kouthat
ccdu = array([-cos(ccdtheta)*cos(ccdphi),-cos(ccdtheta)*sin(ccdphi),sin(ccdtheta)]) #unit vector in plane of ccd (normal to kout)
ccdv = cross(kouthat,ccdu) # other unit vector in plane of ccd (normal to kout)
ccdstep = asarray([ccdu*ccdpixelsize[0],ccdv*ccdpixelsize[1],[0,0,0]])
ccdpixeladdress = mgrid[-ccdpixels[0]/2.:ccdpixels[0]/2.:ccdpixels[0]*1j,-ccdpixels[1]/2.:ccdpixels[1]/2.:ccdpixels[1]*1j]
ccdpixeladdress = concatenate((ccdpixeladdress,zeros((1,ccdpixeladdress.shape[1],ccdpixeladdress.shape[2]))), axis=0)
ccdcoords = tensordot(ccdstep,ccdpixeladdress,axes=[0,0])
ccdcoords += ccdr*kouthat[:,newaxis,newaxis]
ccdshape=ccdcoords.shape[1:3]
#initialize image (storing as a floating point array for now)
ccdimage=zeros(ccdshape)
#
# Begin some physics
#
# calculate phase of probe at each lattice site
# TODO: add polarization (should go in norm as sin**2 term), although see also sigma
localphase=exp(1.j*(kin[0]*coords[0]+kin[1]*coords[1]+kin[2]*coords[2]))*atten[latticeoccupation]
#define function to sum over lattice sites and calculate phase of light at position a
def phasesum(a, b, localphase, kmag): #for the vector a, sum over 3d array b to get complex phases
	shape=b.shape[1:4]
	phase = 1.+0.j
	for x in range(0,shape[0]):
		for y in range(0,shape[1]):
			for z in range(0,shape[2]):
				phase += localphase[x,y,z]*exp(1.j*kmag*sqrt((a[0]-b[0,x,y,z])**2 + (a[1]-b[1,x,y,z])**2 + (a[2]-b[2,x,y,z])**2))
	return phase

# Prefactor to normalize ccd pixel values
ccdnorm = probeintensity/(4*pi*ccdr*ccdr)*ccdpixelarea*exposuretime*wavelength/c/h*DW*DW*sigmamax #account for geometry, etc.
if debug == 0 : # The big loop over the lattice
	for x in range(0,ccdshape[0]):
		for y in range(0,ccdshape[1]):
			ccdimage[x,y]=abs(phasesum(ccdcoords[:,x,y], coords, localphase, kmag))**2.*ccdnorm

# plot
image=plt.imshow(ccdimage,  origin='lower', extent=[-ccdphysicalsize[0]/2,ccdphysicalsize[0]/2,-ccdphysicalsize[1]/2,ccdphysicalsize[1]/2], interpolation='nearest')
#
# Write output file
import csv
csvfile=open(filename, "wb")
csvfile.write(comment)
csvfile.write(latticestring)
writer=csv.writer(csvfile, dialect='excel')
writer.writerow(['kinhat', kinhat[0], kinhat[1], kinhat[2]])
writer.writerow(['kouthat', kouthat[0], kouthat[1], kouthat[2]])
writer.writerow(['Miller indices', Q[0], Q[1], Q[2]])
writer.writerow(['crystal edge length', edge])
writer.writerow(['ccdcenter (xyz)', ccdcenter[0], ccdcenter[1], ccdcenter[2]])
writer.writerow(['ccdcenter (r phi theta)', ccdr, ccdphi, ccdtheta])
writer.writerow(['ccddims', ccdphysicalsize[0], ccdphysicalsize[1]])
writer.writerow(['probe intensity (mW/cm2)', probeintensity/10.])
writer.writerow(['lattice depth (recoil units)', latticedepth/recoil])
writer.writerows(ccdimage)
csvfile.close()
