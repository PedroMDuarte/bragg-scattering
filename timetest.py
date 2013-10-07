
import timeit
import numpy as np

def crossterms( phase ):
    cfactor = 0.
    for i in range(phase.size):
        for j in range(phase.size):
            if i != j:
                cfactor += phase.flat[i]*np.conj(phase.flat[j])
    return cfactor

def squareterms( phase):
    cfactor = 0.
    for i in range(phase.size):  
        cfactor += phase.flat[i]*np.conj(phase.flat[i])
    return cfactor



def sqsum(phase):
    cfactor = np.abs( np.sum( phase ) )**2 - np.sum( np.abs( phase)**2 )
    return cfactor
                
phase = np.exp( 1j * 2.*np.pi*np.random.random(20))

print phase
print

print crossterms(phase)
print squareterms(phase)
print crossterms(phase) + squareterms(phase)

print
print np.abs( np.sum(phase) ) **2
print sqsum(phase)

print

t = timeit.Timer('crossterms(phase)','from __main__ import crossterms,phase')
print round(t.timeit(100),6)

t = timeit.Timer('sqsum(phase)','from __main__ import sqsum,phase')
print round(t.timeit(100),6)
