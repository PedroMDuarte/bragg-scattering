
import sys
import numpy as np
import random

def walk( n = 40**3, p=0.5, step=1.):
    x=0.; 
    for t in range(n):
        u = random.random()
        if ( u<=p) : x= x+1.
        else : x= x-1.
    #return x  
    return x**2./float(n)


import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif'})

s=0.8
figure = plt.figure(figsize=(16.*s,12.*s))
gs = matplotlib.gridspec.GridSpec( 4,5, hspace=0.25, wspace=0.25,
                                   top=0.98,bottom=0.07,left=0.05,right=0.98) 

axL = [[ figure.add_subplot( gs[j,i] ) for i in range(5)] for j in range(4) ] 
         
def plotHIST( ax, x, dat, xmax, **kwargs ): 
    nbins = 50
    if xmax is None:
        n, bins, patches = ax.hist(dat[:,0], nbins, normed=False, \
            alpha=0.75, **kwargs)
    else:
        n, bins, patches = ax.hist(dat[:,0], nbins, (-1., xmax), normed=False, \
            alpha=0.75, **kwargs)

    
    mu = np.mean(dat[:,0]) 
    sig = np.std(dat[:,0])
    ax.axvline( mu, color='black' )  
    ax.axvspan( mu-sig, mu+sig,  facecolor='lightgray', alpha=0.4)
    textstr = '$N^{1/3}=%d$'%N+'\n'+'$p=%0.3f$'%p + '\n' + \
               '$\mu=%.3g$'% mu + '\n' + '$\sigma=%.2g$'%sig + '\n'\
               + '$\sigma/\mu=%.3g$'%(sig/mu)
    ax.text( 0.95,0.95, textstr, \
        va='top', ha='right', fontsize=10,transform=ax.transAxes,\
        bbox = dict(boxstyle='round', facecolor='lightgrey', alpha=1.) )
    return n.max()


def plotDIST( ax, N, p , norm=None):
    n = N**3.
    q = 1.-p
    y = np.linspace( 0.1, 320., 4200)
    P = ( np.sqrt(n/y)/2. ) * 1./ np.sqrt( 2.*np.pi ) / np.sqrt( 4*n*p*q)  \
        * ( np.exp( -0.5 *  ( np.sqrt(y*n) - n*(p-q))**2. / (4*n*p*q) )
          + np.exp( -0.5 *  ( -np.sqrt(y*n) - n*(p-q))**2. / (4*n*p*q) ) )
    if norm is not None:
        P = norm * P / P.max()
    ax.plot( y, P, color='blue', lw=2.5)


 

lc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', \
       'purple', 'gray','orange','green']
fc = [ 'black', 'brown', 'red', 'gold', 'limegreen', 'blue', \
       'purple', 'gray','orange','green']

xmax = [ [9.]*4,  [12.]*4, [12.]*4, [80., 55., 40., 25.], [300., 190., 130., 75.] ] 

Ns = [ 20, 25, 30, 40]
ps = [ 
       [0.5, 0.5015, 0.51, 0.534, 0.58] ,
       [0.5, 0.5015, 0.506, 0.518, 0.545] ,
       [0.5, 0.5015, 0.504, 0.510, 0.525] ,
       [0.5, 0.501, 0.502, 0.506, 0.512],
    ] 
for i,N in enumerate(Ns):
    print "\n",N
    Nstr = '$N=%d$' % N
    Nr = 520
    dirname='noisedat/'
    for j,p in enumerate(ps[i]):
        fname = '_N%02d_%0.3f_%03d.dat'%(N,p,Nr)
        print p,",",
        sys.stdout.flush()
        try: 
            walks = np.loadtxt( dirname + 'walks'+fname)
        except: 
            walks = [] 
            for r in range(Nr):
                walks.append( walk( N**3, p ) )
            #print walks
            walks = np.column_stack((walks,walks))
            np.savetxt( dirname + 'walks'+fname, walks)

        
        maxhist = plotHIST( axL[i][j], N, walks, None, color=lc[i], \
                            facecolor=fc[i], label=Nstr)
        plotDIST( axL[i][j], N, p, norm=maxhist )

    
xlims = [ 
          [(0., 5.), (0.,10), (0.,18), (0.,80), (0.,320.) ], 
          [(0., 5.), (0.,10), (0.,18), (0.,60), (0.,200.) ], 
          [(0., 5.), (0.,10), (0.,14), (0.,50), (0.,140.) ], 
          [(0., 5.), (0.,10), (0.,10), (0.,40), (0.,100.) ],
        ] 

for i,axR in enumerate(axL):
    for j,ax in enumerate(axR):
        ax.set_xlim( xlims[i][j] ) 
        ax.grid()
        if j==0: 
            ax.set_ylabel('counts/bin', fontsize=14)
        if i==len(axL)-1:
            ax.set_xlabel('$S_{\mathbf{\pi}}$',fontsize=18)
    
  #ax.set_xlabel('$L_{\mathrm{AFM}}$',fontsize=16)

 
outfile = 'Noise_hist.png'
figure.savefig(outfile, dpi=250)
