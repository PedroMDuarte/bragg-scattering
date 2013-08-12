import numpy as np
import vec3
import pylab


import braggvectors as bv
import afm

 
N = 40 
A1 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA1))
A2 = afm.crystal(N, nafm, bv.l1064/2, (bv.kin,bv.kA2))

nafms=[2,4,6,7,8,9,10]
for nafm in anfms: 
    inangles = np.linspace(-1600./nafm**1.5,1600./nafm**1.5,20)
    outangles = np.linspace(-1600./nafm**1.5,1600./nafm**1.5,20)

    A1in =[]
    A2in =[]
    for angle in inangles:
        A1.set_kvectors( bv.kinput(angle), bv.kA1, bv.kipol ) 
        A2.set_kvectors( bv.kinput(angle), bv.kA2, bv.kipol ) 
        Nr = 20
        A1in.append( A1.sigma_coh_det( Nr, 0., 0.))
        A2in.append( A2.sigma_coh_det( Nr, 0., 0.))

    A1out =[]
    A2out =[]
    for angle in inangles:
        A1.set_kvectors( bv.kin, bv.koutput(angle), bv.kipol ) 
        A2.set_kvectors( bv.kin, bv.koutput(agnle), bv.kipol ) 
        Nr = 20
        A1in.append( A1.sigma_coh_det( Nr, 0., 0.))
        A2in.append( A2.sigma_coh_det( Nr, 0., 0.))
        







