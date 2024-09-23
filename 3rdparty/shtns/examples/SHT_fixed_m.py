#
#  Copyright (c) 2010-2018 Centre National de la Recherche Scientifique.
#  written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
#
#  nathanael.schaeffer@univ-grenoble-alpes.fr
#
#  This software is governed by the CeCILL license under French law and
#  abiding by the rules of distribution of free software. You can use,
#  modify and/or redistribute the software under the terms of the CeCILL
#  license as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info".
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL license and that you accept its terms.
#

###################################################################
# SHTns Python interface example using fixed m (no fft) functions #
###################################################################

import numpy        # numpy for arrays
import shtns        # shtns module compiled and installed using
                    #   ./configure --enable-python && make && make install

lmax = 31           # maximum degree of spherical harmonic representation.
m = 10              # we work with m=10 only here

sh = shtns.sht(lmax)        # create sht object with given lmax and mmax (orthonormalized)
# sh = shtns.sht(lmax, mmax, mres=2, norm=shtns.sht_schmidt | shtns.SHT_NO_CS_PHASE)    # use schmidt semi-normalized harmonics

nlat, nphi = sh.set_grid()  # build default grid (gauss grid, phi-contiguous)
print(sh.nlat, sh.nphi)     # displays the latitudinal and longitudinal grid sizes.

cost = sh.cos_theta         # latitudinal coordinates of the grid as cos(theta)

ylm10 = sh.spec_array(m)    # a spherical harmonic spectral array for m=10, same as numpy.zeros(sh.lmax +1 - sh.mres*10, dtype=complex)
zlm = sh.spec_array()       # a spectral array for all m, same as numpy.zeros(sh.nlm, dtype=complex)
ym = numpy.zeros(sh.nlat, dtype=complex)    # a spatial array at fixed m

ylm10[1] = 1.0              # set sh coefficient l=11, m=10 to value 1

sh.SH_to_spat_m(ylm10, ym, m)   # Legendre transform at given m

sh.spat_to_SH_m(ym, zlm[sh.idx(m, m):sh.idx(m+1, m+1)], m)  # slicing works as long as it is contiguous

print(zlm[sh.idx(m+1, m)])       # should be one
