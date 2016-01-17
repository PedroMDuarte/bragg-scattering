# bragg-scattering
Calculation of bragg-scattering in a 3D lattice of ultracold atoms

The brute-force numerical calculation for light scattered by an array of spins (at some Q vector and with some degree of spin ordering) is defined in `afm.py`.

The code in `afm.py` is used to evaluate some specific situations, such as 

* crystal size: dependence of signal on number of lattice sites forming the crystal
* DW_100: debye-waller factor calculation (signal vs. tof) at Q=100
* DW_HHH: debye-waller factor at Q=HHH
* detuning_100: signal vs. probe detuning at Q=100
* detuning_HHH: signal vs. probe detuning at Q=HHH
* pbragg_100: signal vs. probe power at Q=100
* pbragg_HHH: signal vs. probe power at Q=HHH
* random_noise: signal as a function of ordering in the crystal.
* rocking: signal as a function of probe beam angle

A write up explaining the calculation can be found in https://github.com/PedroMDuarte/bragg-scattering/blob/master/writeup/bragg.pdf 
