import numpy as np
from scipy import stats

def statdat( dataset , Xcol, Ycol):
  out =[]
  while dataset.shape[0] != 0: 
    #print "# of rows in array = %d" % dataset.shape[0]
     
    Xval = dataset [0, Xcol]
    Yval = []

    to_delete = []
    if np.isnan(Xval):
      to_delete.append(0)
      dataset = np.delete( dataset, to_delete , axis=0)
      continue

    if Xcol >= dataset.shape[1] :
      print "Column index for stat dat is larger than data columns available!"
      exit(1)

    for i in range( dataset.shape[0] ) :
      row = dataset[i,:]
      if row[Xcol] == Xval:
        to_delete.append(i)
        Yval.append( row[Ycol]) 
    dataset = np.delete( dataset, to_delete , axis=0)


    #print "# of rows in array = %d" % dataset.shape[0]
    #print set 
    Yarray = np.array( Yval)
    mean = np.mean(Yarray)
    stddev = np.std(Yarray)
    serror = stats.sem(Yarray)
    pkpk = np.max( Yarray) - np.min( Yarray )
    #print Yval
    out.append( [Xval, mean, stddev, serror, pkpk] ) 
  return np.array( out )

def bindat( dataset , Xcol, Ycol, X0, XF, nbins):
  out =[]
  
  bins,step = np.linspace(X0,XF,nbins, endpoint=False, retstep=True)
  for b in bins:

    if Xcol >= dataset.shape[1]:
      print "Column index for stat dat is larger than data columns available!"
      exit(1)

    Xval = dataset [0, Xcol]
    Yval = []

    to_delete = []
    if np.isnan(Xval):
      to_delete.append(0)
      dataset = np.delete( dataset, to_delete, axis=0) 
      continue

    for i in range( dataset.shape[0] ) :
      row = dataset[i,:]
      if row[Xcol] >= b and row[Xcol]<b+step:
        to_delete.append(i)
        Yval.append( row[Ycol]) 
    dataset = np.delete( dataset, to_delete , axis=0)
 
    #print "# of rows in array = %d" % dataset.shape[0]
    #print set
    Yarray = np.array( Yval)
    if Yarray.size == 1 : 
      mean = np.mean(Yarray)
      stddev = 0.
      serror = 0.
      pkpk = 0.
      out.append( [b+step/2., mean, stddev, serror, pkpk] ) 
    elif Yarray.size > 1 : 
      mean = np.mean(Yarray)
      stddev = np.std(Yarray)
      serror = stats.sem(Yarray)
      pkpk = np.max( Yarray) - np.min( Yarray )
      out.append( [b+step/2., mean, stddev, serror, pkpk] ) 
    else:
      mean = np.nan
      stddev = np.nan
      serror = np.nan
      pkpk = 0.
    #print Yval
  return np.array( out )
      
