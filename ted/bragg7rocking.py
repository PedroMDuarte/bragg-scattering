#
# Calculate Bragg scattering off of atoms in optical lattice
# Ted Corcovilos 20090129
# This version 20090312
# (using scipy module and iPython interpreter)
# all parameters in SI units unless specified
#
# Code assumes atoms are point scatterers and simply adds up the complex phases
# When changing array sizes, remember that coords are at _center_ of cell/pixel
# Changes:
# Version 7 - adding polarization foundations, fix some unit errors
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
debug=0 #1 to skip loops, 0 to do it

def mainloop(state, plane, fileprefix="", debug=1):
	"""mainloop(state, plane, fileprefix, debug)
	Calculates a Bragg scattering image and saves to a csv
	state=0 for ferromagnet, 1 for afm, 2 for random
	plane=0 for (001), 1 for (1/2 1/2 1/2)
	fileprefix is a prefix on the output filename, rest of name based on args
	debug=1 to skip loops, 0 to do it
	"""
	#state=0 # 0 for ferromagnet, 1 for afm, 2 for random
	#plane=0 # 0 for (001), 1 for (1/2 1/2 1/2)
	#debug=0 #1 to skip loops, 0 to do it
	# constants
	#fileprefix="test-"
	comment="Rocking curve - Version 7\n" #comment for csv file
	filenamelist=["fm001","afm001","pm001","fmhhh","afmhhh","pmhhh"]
	filename=fileprefix+filenamelist[state+3*plane]
	print filename
	# universal constants
	c = 2.9979e8 # speed of light m/s
	h = 6.626e-34 # planck constant Js
	hbar = h/2./pi
	kb = 1.3807e-23 # boltzmann constant J/K
	amu = 1.66053886e-27 # atomic mass unit in kg
	re = 2.8179e-15 # classical electron radius in m
	# Li constants
	gamma = 5.92e6 # linewidth of transition in Hz (same for both states?)
	Isat = 25.6 # saturation intensity of transition W/m^2 (using Foot's convention)
	split = 76e6 # zeeman splitting of two magnetic sublevels in Hz (at Feshbach res.)
	mass = 6. * amu # mass of Li-6
	#define lattice
	edge=41 # number of lattice sites along one edge of cube (needs to be odd)
	latticespacing = 532e-9
	recoil = h*h/8./latticespacing/latticespacing/mass # recoil energy in J 
	latticedepth = 10e-6 * kb # lattice depth in J, conveniant value
	shofreq = sqrt(2.*latticedepth*pi**2./mass/latticespacing**2.) #harmonic oscillator frequency of one site
	sholength = sqrt(hbar/mass/shofreq) #harmonic oscillator length
	# setup lattice coordinates
	intmax = (edge-1)/2.; intmin = -1*intmax
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
	probeintensity=0.2*Isat
	rabi = gamma*sqrt(probeintensity/2./Isat) # Rabi frequency in Hz
	wavelength = 671e-9 #probe wavelength
	flux = probeintensity*wavelength/c/h #photons per area per time
	#print "flux: %g" % flux
	detlist = zeros((41,3)) #initialize
	for increment in range(0,41):
		dtheta = pi/180.*(-5.+increment/4.) #scan thetain over +/- 5 deg
		#detlow = (-10.+increment*100./40.)*1.e6
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
			thetain = arccos(wavelength/2./latticespacing) + dtheta
			kinhat=array([sin(thetain),0,cos(thetain)]) #unit vector along probe (looking for (100) peak)
		else: # (1/2 1/2 1/2) case
			alpha = arcsin(sqrt(3)/4*wavelength/latticespacing) #half angle between kin and kout
			thetaQ = arctan(sqrt(2)) #polar angle of kout-kin
			thetain = pi/2. - thetaQ + alpha + dtheta
			phiin = 1.25*pi
			kinhat= asarray([-1.*sin(thetain)/sqrt(2),-1.*sin(thetain)/sqrt(2.),cos(thetain)])
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
		#
		# define ccd pixel locations
		# Use rectangular array centered at kouthat
		exposuretime = 0.000001 #exposure in seconds (used to calculate photon counts...)
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
		DW = exp(-1.*dot(deltak,deltak)*sholength**2./6.) #Debye-Waller factor, assuming harmonic oscillator approximation
		#print "DW: %g" % DW
		Q = deltak*latticespacing/2./pi # scattering momentum Miller indices
		#
		# Define the ccd
		ccdphysicalsize=array([0.05,0.05]) #physical size of ccd array
		ccdpixels = array([25,25]) #number of pixels on ccd
		#ccdpixels = array([3,3]) #small array for debugging
		ccdpixelsize = ccdphysicalsize/ccdpixels
		ccdpixelarea = ccdpixelsize[0]*ccdpixelsize[1]
		#print "ccdpixelarea: %g" % ccdpixelarea
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
		# Next line from Als-Nielsen with mods from Loudon
		lineshape = (3.*wavelength/4./pi)*gamma/2.*(detuning - 0.5j*gamma)/(detuning*detuning+gamma*gamma/4.)
		#atten=lineshape * DW * asarray([0,abs(dot(polin,sigmaplus)), abs(dot(polin,sigmaminus))]) # multiply by some constants and the polarization factor
		atten=lineshape * asarray([0,1,1]) # Structure factor ignores DW and the polarization factor
		# calculate phase of probe at each lattice site
		localphase=exp(1.j*(kin[0]*coords[0]+kin[1]*coords[1]+kin[2]*coords[2]))*atten[latticeoccupation]
		# !!rewrite phasesum to be a sum over both electric field components.  include polarization
		# Prefactor to normalize ccd pixel values:
		#ccdnorm = flux*ccdpixelarea*exposuretime/ccdr/ccdr #normalize to photons/pixel.  account for geometry, etc.
		ccdnorm = 1./edge**6 #normalize to structure factor units (remove N dependence, distance)
		#define function to sum over lattice sites and calculate phase of light at position a
		def phasesum(a, b, localphase, kmag): #for the vector a, sum over 3d array b to get complex phases
			shape=b.shape[1:4]
			phase = 0.+0.j
			for x in range(0,shape[0]):
				for y in range(0,shape[1]):
					for z in range(0,shape[2]):
						dr = sqrt((a[0]-b[0,x,y,z])**2. + (a[1]-b[1,x,y,z])**2. + (a[2]-b[2,x,y,z])**2.)
						phase += localphase[x,y,z]*exp(1.j*kmag*dr)
			return phase

		#for y in range(0,ccdshape[1]):
			#any polarization optics before detector would go in the next line
		detlist[increment,0]=dtheta
		detlist[increment,1]=abs(phasesum(ccdcoords[:,13,13], coords, localphase, kmag))**2.*ccdnorm
		#end loop
	# plot (only works if run as iPython shell script)
	#image=plt.imshow(ccdimage,  origin='lower', extent=[-ccdphysicalsize[0]/2,ccdphysicalsize[0]/2,-ccdphysicalsize[1]/2,ccdphysicalsize[1]/2], interpolation='nearest')
	#
	# Write normalized intensities in col 3 of detlist
	detlist[:,2]=detlist[:,1]/max(detlist[:,1])
	# Write output file
	csvfilename=filename+".csv" #comma separated value file
	datfilename=filename+".dat" #gnuplot ascii matrix datafile
	import csv
	#csv file
	csvfile=open(csvfilename, "wb") #output file in csv format
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
	if debug==1 : #extra info for debug mode
		writer.writerow(['atten',atten[0],atten[1],atten[2]])
		writer.writerow(['localph(0)', localphase[21,21,21]])
		writer.writerow(['phasesum(0)',phasesum(ccdcoords[:,13,13], coords, localphase, kmag)])
		writer.writerow(['phasesum(corner)',phasesum(ccdcoords[:,0,0], coords, localphase, kmag)])
	writer.writerow(['dtheta','signal','normsig'])
	writer.writerows(detlist)
	csvfile.close()
	#datfile
	datfile=open(datfilename, "wb") #output file in gnuplot-ready format (image data only)
	writer=csv.writer(datfile, delimiter=' ', quoting=csv.QUOTE_NONE) #what dialect?
	writer.writerows(detlist)
	datfile.close()
	

def main(argv=None):
	#debug=0
	if argv is None: #standard wrapper stuff
		argv = sys.argv
	for state in range(0,3): #skipping rand for now
		for plane in range(0,2):
			mainloop(state, plane, fileprefix="rocking-", debug=0)
	#mainloop(0,0,"fix-",0)

if __name__ == "__main__":
	sys.exit(main())
