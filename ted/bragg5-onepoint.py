#log# Automatic Logger file. *** THIS MUST BE THE FIRST LINE ***
# #log# DO NOT CHANGE THIS LINE OR THE TWO BELOW
# #log# opts = Struct({'__allownew': True, 'logfile': 'bragg2.py', 'profile': 'scipy', 'pylab': 1})
# #log# args = []
# #log# It is safe to make manual edits below here.
# #log#-----------------------------------------------------------------------
#
# Calculate Bragg scattering off of atoms in optical lattice
# Ted Corcovilos 20090129
# (using scipy module and iPython interpreter)
# all parameters in SI units unless specified
#
# Code assumes atoms are point scatterers and simply adds up the complex phases
# When changing array sizes, remember that coords are at _center_ of cell/pixel
#
# import libraries scipy and matplotlib
#import psyco ; psyco.log() ; psyco.profile()
#from psyco.classes import *
import sys, math
import numpy as np
from scipy import * 
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# debug?
debug=0 #1 to skip loop, 0 to do it
# constants
comment="afm case\n" #comment for csv file
filename="detuningafmhhh.csv" #output filename
gamma = 5.92e6 # linewidth of transition in Hz
Isat = 25.6 # saturation intensity of transition W/m^2 (using Foot's convention)
split = 76e6 # zeeman splitting of two magnetic sublevels in Hz (typical value at high field)
c = 2.9979e8 # speed of light m/s
h = 6.626e-34 # planck constant Js
hbar = h/2./pi
kb = 1.3807e-23 # boltzmann constant J/K
#define lattice
#print('defining things...')
latticespacing = 532e-9
intcoords=mgrid[-25.:25.:51j,-25.:25.:51j,-25.:25.:51j] # 50 x 50 x 50 lattice centered at origin
#choose which sites are occupied
#(in future may use to define which state the atoms at each site are in)
#beamradius = 60e-6 #assuming all sites inside a sphere of this radius are occupied
#latticeoccupation = less(coords[0]*coords[0]+coords[1]*coords[1]+coords[2]*coords[2],beamradius**2)
# define lattice as cube, for simplicity
#
#
# latticeoccupation = 0 for empty
#                     1 for state 1
#                     2 for state 2
#latticeoccupation = ones(intcoords.shape[1:4], dtype=np.int) #ferromagnetic case, crystal is cube
latticeoccupation = asarray(1 + (intcoords[0]+intcoords[1]+intcoords[2])%2, dtype=np.int) #antiferromagnetic case (alternate 1 and 2)
#latticeoccupation = np.random.randint(1, 3, intcoords.shape[1:4]) #paramagnet
coords = intcoords*latticespacing # scale to physical units
# define incoming probe beam
probepower = 5e-3
probebeamwaist = 100e-6 #Yikes!
probeintensity = probepower/pi/probebeamwaist**2
rabi = gamma*sqrt(probeintensity/2./Isat) # Rabi frequency in Hz
wavelength = 671e-9 #probe wavelength
detlow = split/2. # probe detuning from lower transition
detuning = array([0, detlow, detlow-split]) # indexed for "latticeoccupation"
# stuff about probe polarization
sigma = 3.*wavelength**2/2/pi # TODO: need to multiply by reduced matirix element, depends on polarization and state
lineshape = rabi*(detuning + 0.5j*gamma)/(detuning*detuning+rabi**2/2+gamma**2/4)
#atten = exp(1.0j*lineshape)
atten = lineshape
atten[0] = 0 # reset empty case
#reset atten for debug
#atten = ones(atten.shape, dtype='complex')
#atten = asarray([0., 1., -1.])
kmag=2.*pi/wavelength
#pick unit vectors
# next line is for (001)
#kinhat=array([sqrt(1.-(wavelength/2./latticespacing)**2),0.,-1*wavelength/2./latticespacing]) #unit vector along probe (looking for (100) peak)
# next few lines are for (1/2 1/2 1/2)
alpha = arcsin(sqrt(3)/4*wavelength/latticespacing) #half angle between kin and kout
thetaQ = arctan(sqrt(2)) #polar angle of kout-kin
thetain = pi/2. - thetaQ + alpha
phiin = 1.25*pi #absorb into next line
kinhat= asarray([-1.*sin(thetain)/sqrt(2),-1.*sin(thetain)/sqrt(2),cos(thetain)])
#end 1/2 1/2 1/2
kin= kmag * kinhat
# TODO: add polarization of incoming beam
#
# define ccd pixel locations
# Use rectangular array centered at kouthat
exposuretime = 1. #exposure in seconds (used to calculate photon counts...)
ccdr = 0.5 #distance from origin to center of ccd
#next lines for (001)
#ccdphi = 0 #azimuthal angle of ccd center
#ccdtheta = arccos(wavelength/2./latticespacing) #polar angle of ccd center (this is (1 0 0) reflection)
#next lines for (1/2 1/2 1/2)
ccdtheta=pi/2. - thetaQ - alpha #see defn of thetaQ and alpha above
ccdphi = 1.25*pi
#end 1/2 1/2 1/2
kouthat = array([cos(ccdphi)*sin(ccdtheta), sin(ccdphi)*sin(ccdtheta), cos(ccdtheta)]) #unit vector from origin to center of ccd
kout=kouthat*kmag
deltak=kout-kin
Q=deltak*latticespacing/2/pi
#
# Commenting out ccdarray stuff.   Just want to look at intensity at one point.
#
ccdphysicalsize=array([0.03,0.03]) #physical size of ccd array
ccdpixels = array([25,25]) #number of pixels on ccd
ccdpixelsize = ccdphysicalsize/ccdpixels
ccdpixelarea = ccdpixelsize[0]*ccdpixelsize[1]
ccdcenter = ccdr*kouthat
#ccdu = array([-cos(ccdtheta)*cos(ccdphi),-cos(ccdtheta)*sin(ccdphi),sin(ccdtheta)]) #unit vector in plane of ccd
#ccdv = cross(kouthat,ccdu) # other unit vector in plane of ccd
#ccdstep = asarray([ccdu*ccdpixelsize[0],ccdv*ccdpixelsize[1],[0,0,0]])
#ccdpixeladdress = mgrid[-ccdpixels[0]/2.:ccdpixels[0]/2.:ccdpixels[0]*1j,-ccdpixels[1]/2.:ccdpixels[1]/2.:ccdpixels[1]*1j]
#ccdpixeladdress = concatenate((ccdpixeladdress,zeros((1,ccdpixeladdress.shape[1],ccdpixeladdress.shape[2]))), axis=0)
#ccdcoords = ccdpixeladdress*ccdustep[:,newaxis,newaxis] + ccdpixeladdress*ccdvstep[:,newaxis,newaxis]
#ccdcoords = tensordot(ccdstep,ccdpixeladdress,axes=[0,0])
#ccdcoords += ccdr*kouthat[:,newaxis,newaxis]
#ccdshape=ccdcoords.shape[1:3]
#initialize image (storing as a floating point array for now)
#ccdimage=zeros(ccdshape)
# calculate phase of probe at each lattice site
# TODO: add dispersion, abosorption, and spin state
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

ccdnorm = sigma/(4*pi*ccdr*ccdr)*ccdpixelarea*exposuretime*wavelength/c/h #account for geometry, etc.
# set up loop over detuning
detlist = asarray([linspace(-1.,2.,25)*split,zeros(25)])
#detlow = split/2. # probe detuning from lower transition
if debug == 0:
	for num in range(0,detlist.shape[1]):
		d=detlist[0,num]
		detuning = array([0, d, d-split]) # indexed for "latticeoccupation"
		lineshape = rabi*(detuning + 0.5j*gamma)/(detuning*detuning+rabi**2/2+gamma**2/4)
		#atten = exp(1.0j*lineshape)
		atten = lineshape
		atten[0] = 0 # reset empty case
		localphase=exp(1.j*(kin[0]*coords[0]+kin[1]*coords[1]+kin[2]*coords[2]))*atten[latticeoccupation]
		detlist[1,num]=abs(phasesum(ccdcenter, coords, localphase, kmag))*ccdnorm
	
	
#if debug == 0:#
#	for x in range(0,ccdshape[0]):
#		for y in range(0,ccdshape[1]):
#			# print('starting',x,y)
#			ccdimage[x,y]=abs(phasesum(ccdcoords[:,x,y], coords, localphase, kmag))*ccdnorm



value=abs(phasesum(ccdcenter, coords, localphase, kmag))*ccdnorm

# plot
#image=plt.imshow(ccdimage,  origin='lower', extent=[-ccdphysicalsize[0]/2,ccdphysicalsize[0]/2,-ccdphysicalsize[1]/2,ccdphysicalsize[1]/2], interpolation='nearest')
import csv
csvfile=open(filename, "wb")
csvfile.write(comment)
writer=csv.writer(csvfile, dialect='excel')
writer.writerow(['kinhat', kinhat[0], kinhat[1], kinhat[2]])
writer.writerow(['kouthat', kouthat[0], kouthat[1], kouthat[2]])
writer.writerow(['Q', Q[0], Q[1], Q[2]])
writer.writerow(['ccdcenter (xyz)', ccdcenter[0], ccdcenter[1], ccdcenter[2]])
writer.writerow(['ccdcenter (r phi theta)', ccdr, ccdphi, ccdtheta])
#writer.writerow(['ccddims', ccdphysicalsize[0], ccdphysicalsize[1]])
writer.writerow(['detuning (Hz)', 'signal'])
writer.writerows(detlist.swapaxes(0,1))
csvfile.close()
