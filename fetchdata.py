
import numpy as np
from statarray import statdat

dkeys = { 'ai' : 0 ,
          'af' : 1 , 
          'latsize' : 2, 
          'afmsize' : 3,
          'det' : 4,
          'v0[0]' : 5, 
          'v0[1]' : 6, 
          'v0[2]' : 7, 
          'kipol[0]' : 8,
          'kipol[1]' : 9,
          'sigma' : 10, 
          'sigma_inelastic' : 11,
          'polsum' : 12,
          'debye-waller' : 13, 
          'alpha' : 14,
          'beta' : 15,
          'C' : 16, 
          'S' : 17 ,
} 

         
def fetch_data_A1A2( filters, xkey, datfile, ykey='sigma'):
    dat = np.loadtxt(datfile)  
    # First get A1 and A2 entries
    a1 = dat[ dat[:, dkeys['af']] == 1486.41 ]  # ANDOR1
    a2 = dat[ dat[:, dkeys['af']] == 0. ] # ANDOR2
  
    # Get the entries that have the correct values for all the filters 
    for k in filters.keys():
      a1 = a1[ a1[:,dkeys[k]] == filters[k] ]
      a2 = a2[ a2[:,dkeys[k]] == filters[k] ]
    

    # Return the elastic cross section data as a function of xkey
    x = dkeys[xkey] 
    y = dkeys[ykey] 
    sort1 = np.argsort(a1[:,x])
    sort2 = np.argsort(a2[:,x])
   
    a1 = a1[sort1,np.array([[x],[ y ]])].T
    a2 = a2[sort2,np.array([[x],[ y ]])].T
 
    a1 = statdat( a1, 0, 1 )
    a2 = statdat( a2, 0, 1 )

    return a1, a2



