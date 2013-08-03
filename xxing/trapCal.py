#!/usr/bin/python
# Simulation for traps used in EMT1

import scipy.constants as C
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
from matplotlib.ticker import OldScalarFormatter


class point:
    '''Define the coordinate of a point in our system.
    
    x, y and z directions are specified in figure 2.9 of Yean's Ph.D thesis. '''
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.rho = np.sqrt(x**2+y**2)   #radial coordinate
        self.r =np.sqrt(x**2+y**2+z**2) #vector length

    def __add__(self, other):
        ''' Overloading the + sign of point object.'''
        return point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        '''overloading the - sign of point object.'''
        return point(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def getdir(self):
        '''Return the direction of vector point'''
        if self.r == 0:
            print "Vector length is zero!!!"
            return direction(0,0)
        else:
            theta = np.arccos(self.z/self.r)
            if self.rho ==0:
                phi = 0
            else:
                phi = np.arccos(self.x/self.rho)
                if self.y < 0:
                    phi = -phi
            return direction(theta,phi)


class direction:
    '''Define the direction by polar angle theta and azimuth angle phi'''
    def __init__(self, theta, phi):
        self.theta = theta      #polar angle
        self.phi = phi          #azimuth angle

    def getposition(self,l):
        '''Convert from Spherical coordinate system (l,theta, phi) to Cartesian coordinate system (x,y,z).'''
        x = l*np.sin(self.theta)*np.cos(self.phi)
        y = l*np.sin(self.theta)*np.sin(self.phi)
        z = l*np.cos(self.theta)
        return point(x, y, z)
    
    def getprojection(self, p):
        '''Return the projection length of vector p on startpoint.'''
        e = self.getposition(1) # unit vector of direction
        e_array = np.array([e.x,e.y,e.z])
        p_array = np.array([p.x,p.y,p.z])
        if np.size(p_array.shape) == 3:
            permu = (1,2,0)
            p_array = np.transpose(p_array, permu)
        elif np.size(p_array.shape) == 2:
            p_array = np.transpose(p_array)
        
        return np.inner(e_array, p_array)

class line:
    '''Define a line in our system from a point startpoint and the positive direction pdir.'''
    def __init__(self,startpoint, pdir):
        self.startpoint = startpoint
        self.pdir = pdir

    def getposition(self, l):
        '''Return the position of a point on the line whose distance from startpoint is l.'''
        return self.startpoint + self.pdir.getposition(l)

class crosssection:
    '''Define a crosssection from a point startpoint and two directions hdir and vdir.'''
    def __init__(self,startpoint, hdir, vdir):
        self.startpoint = startpoint
        self.hdir = hdir
        self.vdir = vdir

    def getposition(self, h, v):
        '''Return the position of a point on the crosssection whose coordinate is specified as (h,v).'''
        return self.startpoint + self.hdir.getposition(h) + self.vdir.getposition(v)

    def getz(self):
        vech = self.hdir.getposition(1)
        vecv = self.vdir.getposition(1)
        vecz_array = np.cross(np.array([vech.x,vech.y,vech.z]),np.array([vecv.x,vecv.y,vecv.z]))
        vecz = point(vecz_array[0],vecz_array[1],vecz_array[2])
        return vecz.getdir()

class atom:
    '''Define the properties of the atom in use.'''

    def __init__(self, name, mass, drate, la0):
        self.name = name        #atom name
        self.mass = mass        #atom mass
        self.drate = drate      #damping rate of the upper level
        self.la0 = la0          #resonance wavelength
        self.f0 = C.c/self.la0  #resonance frequency
        self.omg = self.f0*2*C.pi       #resonance anglar frequency

class beam:
    ''' Beams used for the optical dipole traps.

    Contains all the parameters for an specific beam. The getI(p) function gives the intensity at point p. '''
    def __init__(self, P, crs, la, wx, wy):
        self.P = P      # power
        self.startpoint = crs.startpoint  # beam center point
        self.dirx = crs.hdir    # x,y,z direction in beam setup
        self.diry = crs.vdir
        self.dirz = crs.getz()
        self.la = la    # wavelength lambda
        self.wx = wx    # waist x
        self.wy = wy    # waist y
        self.I0 = self.P*2/C.pi/self.wx/self.wy      #intensity at waist center
        self.f = C.c/self.la    # frequency
        self.omg = self.f*2*C.pi        # angular frequency
        self.k = 2*C.pi/self.la/50         # wave number k
        self.Rzx = C.pi*self.wx**2/self.la       # Rayleigh range x
        self.Rzy = C.pi*self.wy**2/self.la       # Rayleigh range y

    def getI(self, p):
        '''Return the intensity at point p.'''
        px = self.dirx.getprojection(p - self.startpoint)     # coordinate of p in beam setup
        py = self.diry.getprojection(p - self.startpoint)
        pz = self.dirz.getprojection(p - self.startpoint)

        wzx = self.wx*np.sqrt(1+(pz/self.Rzx)**2)        #beam crosssection size x at z
        wzy = self.wy*np.sqrt(1+(pz/self.Rzy)**2)        #beam crosssection size y at z

        return self.P*2/C.pi/wzx/wzy*np.exp(-2*px**2/wzx**2-2*py**2/wzy**2)        #intensity at p

class couplebeams(beam):
    ''' Beams used for the optical dipole traps.

    Contains all the parameters for an specific beam. The getI(p) function gives the intensity at point p. '''
    def __init__(self, P, crs, la, wx, wy, refl =0, lattice ='off', dis = 0):
        beam.__init__(self, P, crs, la, wx, wy)
        self.refl = refl        #reflection rate
        self.lattice = lattice  #form lattice or not
        self.dis = dis          #distance from mirror to startpoint

    def getI(self, p):
        '''Return the intensity at point p.'''
        px = self.dirx.getprojection(p - self.startpoint)     # coordinate of p in beam setup
        py = self.diry.getprojection(p - self.startpoint)
        pz = self.dirz.getprojection(p - self.startpoint)
        wzx = self.wx*np.sqrt(1+(pz/self.Rzx)**2)        #beam crosssection size x at z
        wzy = self.wy*np.sqrt(1+(pz/self.Rzy)**2)        #beam crosssection size y at z

#        print px, '\n', py, '\n', pz,'\n'

        Iin = self.P*2/C.pi/wzx/wzy*np.exp(-2*px**2/wzx**2-2*py**2/wzy**2)        #intensity at p
        if self.lattice == 'off':
            return Iin * (1 + self.refl)
        elif self.lattice == 'on':
            invRzx = pz/(pz**2+self.Rzx**2)     #Rz^(-1) in x direction, where Rz is Radius of curvature
            invRzy = pz/(pz**2+self.Rzy**2)     #Rz^(-1) in y direction
            return Iin * (2*(1-self.refl)+4*self.refl*np.cos(self.k*(self.dis - pz- px**2*invRzx/2- py**2*invRzy/2))**2)    # neglected zeta = arctan(z/Rz)
        else:
            print 'lattice could be either \'on\' or \'off\'.'
            return 0

class potential:
    '''Define all kinds of potentials the atom feels.'''

    def __init__(self, atom):
        self.atom = atom

class dtrap(potential):
    '''Define the optial dipole trap which is induced by beam.'''
    
    def __init__(self, atom, beam):
        potential.__init__(self, atom)
        self.beam = beam        #the beam used for generating dipole trap

    def getU(self, p):
        '''Return the potential depth at point p due to dtrap.'''
        return -3*C.pi*C.c**2/2/self.atom.omg**3*self.beam.getI(p)*self.atom.drate*(1/(self.atom.omg-self.beam.omg)+1/(self.atom.omg+self.beam.omg))


class gravity(potential):
    '''Define the gravity potential.'''

    def __init__(self, atom, direction):
        potential.__init__(self, atom)
        self.direction = direction      #direction of gravity

    def getU(self, p):
        '''Return the gravity potential at point p.'''
        return -self.atom.mass*C.g*(p.x*np.cos(self.direction.phi)+p.y*np.sin(self.direction.phi))

class cbinpot(potential):
    '''Define the combined potential of several different potentials specified in potlist.'''

    def __init__(self, atom, potlist):
        potential.__init__(self, atom)
        self.potlist = potlist  #potentials taken into account
        #-----------------------parameters used to show the 3 2D view of potential
        self.hlength = 0.000003   #horizontal length, the direction is horizontal and perpendicular to z
        self.vlength = 0.000003   #vertical length
        self.zlength = 0.00002   #z direction length
        self.prsn = 50          #points number in each direction
        #-----------------------parameters used for fitting and frequency calculation in a specific direction
        self.llength_fit = 0.000001      #line length in fitting
        self.prsn_fit = 20              #number of point in fitting

        for pot in self.potlist:        #be sure to use the same atoms in different potentials
            if self.atom.name != pot.atom.name:
                pot.atom = self.atom
                print "Different atoms used in calculation!!!"

    def getU(self, p):
        '''Return the combined potential depth at point p.''' 
        U=0;
        for pot in self.potlist:
            U = U + pot.getU(p)
        return U
   
    def calcline(self, ln, llength, prsn):
        '''Calculate the combined potential depth along a specific line ln.
        
        llength is the line length calculated, and prsn is the number of points calculated.'''
        l = np.linspace(-llength/2,llength/2,prsn+1)
        p = ln.getposition(l)
        U = self.getU(p)/C.k*10**6
        return l, U

    def showline(self, ax, ln, llength, prsn):
        '''Show the potential depth along line ln in axel ax.'''
        l, U = self.calcline(ln, llength, prsn)
        PT = ax.plot(l,U,'-')
        return l, U, PT
    
    def getfreq(self, ln):
        '''Return the trap frequency along line ln.'''
        l, U = self.calcline(ln,self.llength_fit,self.prsn_fit)
        Ufit = np.polyfit(l,U,2)
        return (Ufit[0]*2/self.atom.mass*C.k/10**6)**(1/2.0)/2.0/C.pi

    def getTF(self, N):
        '''Return the Fermi temperature of this trap with N atoms in it.'''
#        startpoint = point(0,0,0)
#        vdir = direction(C.pi/2,C.pi*1/3)
#        hdir = direction(C.pi/2,-C.pi/6)
#        zdir = direction(0,0)
        lnv = line(startpoint, vdir)
        lnh = line(startpoint, hdir)
        lnz = line(startpoint, zdir)
        fv = self.getfreq(lnv)
        fh = self.getfreq(lnh)
        fz = self.getfreq(lnz)
#        print fv,fh,fz
        omega = (fv*fh*fz)**(1/3.0)*2.0*C.pi    #mean anglar frequency of the trap
        return C.hbar*omega/C.k*(6*N)**(1/3.0)*10**6

    def showcrosssection(self, ax, crs, hlength, vlength, prsn, style):
        '''Show the potential distribution along a specific crosssection crs in axel ax.

        hlength and vlength define the size of the plan. prsn is the number calculated in each direction.'''
        h = np.linspace(-hlength/2,hlength/2,prsn)
        v = np.linspace(-vlength/2,vlength/2,prsn)
        H,V = np.meshgrid(h,v)
#        print H *10**6, '\n',V*10**6
        p = crs.getposition(H,V)
#        print p.x*10**6,'\n', p.y*10**6,'\n', p.z*10**6
        Z = self.getU(p)/C.k*10**6

        if style == 'contour':
            levels = np.linspace(np.rint(np.min(Z)),np.rint(np.max(Z)),6)
            CS = ax.contour(H,V,Z,levels)
            ax.clabel(CS,inline=1,fontsize =10,fmt = '%1.1f')
        elif style == 'image':
            ax.plot_surface(H, V, Z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased= True)
            ax.contour(H, V, Z, zdir='x', offset = -hlength*2/3,levels = np.array([0]) )
            ax.contour(H, V, Z, zdir='y', offset = vlength*2/3, levels = np.array([0]) )
            ax.contour(H, V, Z, zdir='z', offset = -300, cmap=cm.jet)

        else:
            print "style should be either \'contour\' or \'image\'."

    def fig3view(self, style):
        '''Show the 3 2D views of the potential.'''
        
        if style == 'contour':
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
            col = 4
            figsz = (14,6)
        elif style == 'image':
            col = 5
            figsz = (14,4)
        else:
            print "style should be either \'contour\' or \'image\'."
            return
        
        plt.figure(figsize = figsz)
        
        levels = np.linspace(-50,-0,6)
        ax1 = plt.subplot2grid((2,col),(0,0),colspan= col-1)
        crs = crosssection(startpoint,zdir,vdir)
#        print crs.getz().theta/C.pi, crs.getz().phi/C.pi
        cpot.showcrosssection(ax1, crs, self.zlength, self.vlength, self.prsn, style)
        ax1.set_xlabel('z direction')
        ax1.set_ylabel('vertical direction')
        ax1.xaxis.set_major_formatter(OldScalarFormatter())
        ax1.yaxis.set_major_formatter(OldScalarFormatter())

        ax2 = plt.subplot2grid((2,col),(1,0),colspan= col-1)
        crs = crosssection(startpoint,zdir,hdir)
        cpot.showcrosssection(ax2, crs, self.zlength, self.hlength, self.prsn, style)
        ax2.set_xlabel('z direction')
        ax2.set_ylabel('horizontal direction')
        ax2.xaxis.set_major_formatter(OldScalarFormatter())
        ax2.yaxis.set_major_formatter(OldScalarFormatter())

        ax3 = plt.subplot2grid((2,col),(0,col-1),colspan= 1)
        crs = crosssection(startpoint,hdir,vdir)
        cpot.showcrosssection(ax3, crs, self.hlength, self.vlength, self.prsn, style)
        ax3.set_xlabel('horizontal direction')
        ax3.xaxis.set_major_formatter(OldScalarFormatter())
        ax3.yaxis.set_major_formatter(OldScalarFormatter())
    
        plt.show()
        
        return

    def showtubedevi(self, spoint, direction, rlength = 0.00001064):
        llength = 0.0001
        delta = 0.000000532
        r = np.arange(-rlength, rlength+delta, delta)
        plt.figure(figsize = (12,8))

        ax = plt.subplot(1,1,1)
        dU = np.zeros(r.shape[0])
        i = 0
        for ri in r:
            ln = line(spoint + direction.getposition(ri), zdir)
            l,U = self.calcline(ln,llength, 50)
            dU[i] = np.max(U)-np.min(U)
            i = i+1
        plt.plot(r/delta,dU,'bo')
        plt.axis([np.min(r/delta)-1, np.max(r/delta)+1,0,np.ceil(np.max(dU)/0.01)*0.01])
#    ln = line(startpoint +gdir.getposition(7*delta), zdir)
#    cpot.showline(ax,ln,llength,50)
        plt.show()

#------------------------------------------------------    
if __name__ == "__main__":

    li6_name = 'lithium6'
    li6_mass = 6.015121*1.66054*10**(-27)   #kg
    li6_drate = 2*C.pi*5.9*10**6    #5.9MHz
    li6_la =671*10**(-9)   #671nm
    li6 = atom(li6_name, li6_mass, li6_drate,li6_la)

    startpoint = point(0,0,0)
    hdir = direction(C.pi/2,-C.pi/6)    #horizontal direction(perpendicular to z)
#    vdir = direction(C.pi/2*89/90,C.pi*1/3)   #vertical direction(opposite to gravity)
    vdir = direction(C.pi/2,C.pi*1/3)   #vertical direction(opposite to gravity)
    xdir = direction(C.pi/2,0)          #x direction, 30 degree from hdir
    ydir = direction(C.pi/2,C.pi/2)     #y direction, 30 degree from vdir
#    zdir = direction(C.pi/2/90,C.pi*4/3)               #z direction
    zdir = direction(0,0)               #z direction
    xpdir = direction(C.pi/2,C.pi/6)    #x' direction, 30 degree from x
    ypdir = direction(C.pi/2,C.pi*2/3)  #y' direction, 30 degree from y

    sb_P = 0*1       # Single Beam power, in Watt
    sb_la = 1080*10**(-9)  # lambda, 1.08um
    sb_w = 26.3*10**(-6)  # waist, 26.3um
    sb_crs = crosssection(startpoint, hdir, vdir)
    sb = couplebeams(sb_P, sb_crs, sb_la, sb_w, sb_w)
    sbtrap = dtrap(li6, sb)        #single beam trap

    reflection = 1 
    latticestatus = 'off'
    
    north_P = 10 #North beam power, in Watt
    north_la = 1064*10**(-9)
    north_wx = 54*10**(-6)
    north_wz = 236*10**(-6)
    north_crs = crosssection(startpoint, xpdir, zdir)
    north = couplebeams(north_P, north_crs, north_la, north_wx, north_wz, refl = reflection, lattice = latticestatus)
    ntrap = dtrap(li6, north)   #north beam trap
    
    south_P = north_P #South beam power, in Watt
    south_la = 1064*10**(-9)
    south_wx = 54*10**(-6)
    south_wz = 236*10**(-6)
    south_crs = crosssection(startpoint, ypdir, zdir)
    south = couplebeams(south_P, south_crs, south_la, south_wx, south_wz, refl = reflection, lattice = latticestatus)
    strap = dtrap(li6, south)   #south beam trap
    
    green_P = 3.85*north_P #Green beam power, in Watt
    green_la = 532*10**(-9)
    green_wx = 54*10**(-6)
    green_wz = 236*10**(-6)
    green_crs = crosssection(startpoint, ydir, zdir)
    green = beam(green_P, green_crs, green_la, green_wx, green_wz)
    gtrap = dtrap(li6, green)   #green beam trap

    gdir = direction(C.pi/2,-C.pi*2/3)  #gravity direction
    grvty = gravity(li6,gdir)    #gravity
    
    cpot = cbinpot(li6,[strap,ntrap])   #combined potential
    
#    N = 1.7*10**6       #atom number
#    print cpot.getTF(N)
#    cpot.fig3view('image')
#    cpot.showtubedevi(startpoint, gdir, 0.00000532)
#    lnz = line(startpoint, zdir)
#    lnh = line(startpoint, hdir)
#    print cpot.getfreq(lnz), cpot.getfreq(lnh)

    from mpl_toolkits.mplot3d import Axes3D       
    from matplotlib import cm
    
#    fig = plt.figure(figsize = (10,10))
        
#    ax = fig.gca(projection='3d')
#    crs = crosssection(startpoint, xpdir, ypdir)
#    cpot.showcrosssection(ax, crs, 0.0002 ,0.0002 , 10, 'image')
    
#    plt.show()
    fig = plt.figure(figsize = (8,6))
        
    ax = fig.gca(projection='3d')
    crs = crosssection(startpoint, xpdir, ypdir)

    xplength = 0.00015
    yplength = xplength
    prsn = 100
    cpot.showcrosssection(ax, crs, xplength , yplength , prsn, 'image')
    
    ax.auto_scale_xyz([-xplength*2/3, xplength*2/3], [-yplength*2/3, yplength*2/3], [-300,0])
#    ax.get_xaxis().set_ticks([])
#    ax.get_yaxis().set_ticks([])

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.zaxis.set_ticklabels([])
#    ax.axis('off')
    ax.view_init(30,-30)
    plt.show()

