
import numpy as np
import braggvectors


class crystal():
    """Contains information about crystal. 
       Used to calculate scattering of light"""
    def __init__(self, latsize, lcorr, spacing, kin_out=None):
        # crystal settings
        self.latsize = latsize # number of lattice sites on the side
        self.lcorr = lcorr # correlation length in units of lattice spacing

        self.a = spacing # lattice spacing
        self.mass = 6.
        #self.x, self.y, self.z = np.ogrid[ 0:size, 0:size, 0:size ]
        self.x, self.y, self.z = np.mgrid[ 0:self.latsize, \
                                           0:self.latsize, \
                                           0:self.latsize ]

        # default values for scattering parameters
        self.set_detuning( 0. )
        self.set_pbragg( 250. )
        self.lBragg = braggvectors.l671
        self.sunits( ) 

        if kin_out == None:
            self.set_kvectors( braggvectors.kin, braggvectors.kout ,\
                               braggvectors.kipol)
        else:
            self.set_kvectors( kin_out[0], kin_out[1],\
                               braggvectors.kipol)

        self.set_v0( [20.,20.,20.] ) 
