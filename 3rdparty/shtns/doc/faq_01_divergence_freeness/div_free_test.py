#! /usr/bin/env python3

import numpy as np
import shtns
import sys
import stat
import os

class Spharmt(object):
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models.  Provides an interface to shtns compatible
    with pyspharm (pyspharm.googlecode.com).
    """
    def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian'):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,0)
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,0)
        self.lats = np.arcsin(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)
        self.nlons = nlons
        self.nlats = nlats
        self.ntrunc = ntrunc
        self.nlm = self._shtns.nlm
        self.degree = self._shtns.l
        self.lap = -self.degree*(self.degree+1.0)#.astype(np.complex)
        self.invlap = np.copy(self.lap)
        self.invlap[0] = 1
        self.invlap = 1.0/self.invlap
        self.invlap[0] = 0
        self.rsphere = rsphere
        
        if not 'nooutput' in sys.argv:
            print("N: "+str(self.nlons)+", "+str(self.nlats))
            print("Mtrunc: "+str(self.ntrunc))
            print("Nlm: "+str(self.nlm))

    def grdtospec(self,data):
        """compute spectral coefficients from gridded data"""
        return self._shtns.analys(data)
    def spectogrd(self,dataspec):
        """compute gridded data from spectral coefficients"""
        return self._shtns.synth(dataspec)
    def getuv(self,vrtspec,divspec):
        """compute wind vector from spectral coeffs of vorticity and divergence"""
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)
    def getvrtdivspec(self,u,v):
        """compute spectral coeffs of vorticity and divergence from wind vector"""
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec
    def getuv_from_stream(self,strmspec):
        """compute wind vector from spectral coeffs of stream function"""
        return self._shtns.synth_grad(strmspec)

file_ext = ""

output_file_name="{:s}.csv"

def savefile(data, name, t=0):
    d = x.spectogrd(data)
    d = np.flip(d, 0)
    np.savetxt(output_file_name.format(name, t/(60*60)), d, delimiter="\t")



def savefileg(data, name, t=0):
    d = np.flip(data, 0)
    np.savetxt(output_file_name.format(name, t/(60*60)), d, delimiter="\t")



nlats = 64

if len(sys.argv) > 1:
    nlats = int(sys.argv[1])

nlons = nlats*2

if len(sys.argv) > 2:
    nlons = int(sys.argv[2])

ntrunc = int(nlons/3)  # spectral truncation (for alias-free computations)

rsphere = 1 #2314123


# setup up spherical harmonic instance, set lats/lons of grid
x = Spharmt(nlons, nlats, ntrunc, rsphere, gridtype='gaussian')
lons, lats = np.meshgrid(x.lons, x.lats)

#
# SETUP U/V here from the paper
#
# Case-1: page 5
#

# u, v, psi
ug = np.zeros((nlats, nlons))
vg = np.zeros((nlats, nlons))
psi = np.zeros((nlats, nlons))
for j in range(nlats):
    for i in range(nlons):
        # Theta
        the = x.lats[j]
        # Lambda
        lam = x.lons[i]

        # We drop the time dependent terms as well as 'k'
        ug[j,i] = np.sin(lam/2.0)**2 * np.sin(2.0*the)
        vg[j,i] = 1.0/2.0 * np.sin(lam) * np.cos(the)
        psi[j,i] = np.sin(lam/2.0)**2 * np.cos(the)**2





#
# Using velocity forward/backward transformation
#
vrtspec, divspec =  x.getvrtdivspec(ug,vg)
ug_vrtdiv, vg_vrtdiv = x.getuv(vrtspec, divspec)

lmax_ug_vrtdiv = np.max(np.abs(ug-ug_vrtdiv))
lmax_vg_vrtdiv = np.max(np.abs(vg-vg_vrtdiv))
if not 'nooutput' in sys.argv:
    print('Lmax error on ug (based on vrtdiv transf.):', lmax_ug_vrtdiv)
    print('Lmax error on vg (based on vrtdiv transf.):', lmax_vg_vrtdiv)



#
# Using stream function
#
strm = x.grdtospec(psi)
ug_psi,vg_psi = x.getuv_from_stream(strm)     # divergence-free built-in

lmax_ug_psi = np.max(np.abs(ug-ug_psi))
lmax_vg_psi = np.max(np.abs(vg-vg_psi))



if not 'nooutput' in sys.argv:
    print('Lmax  error on ug (based on psi):', lmax_ug_psi)
    print('Lmax error on vg (based on psi):', lmax_vg_psi)


if not 'nooutput' in sys.argv:
    #
    # Output data to .csv files
    #
    savefileg(ug_psi, "output_ug_psi", 0)
    savefileg(vg_psi, "output_vg_psi", 0)


#
# Test for divergence freeness
#

vrtspec_, divspec_ =  x.getvrtdivspec(ug_vrtdiv,vg_vrtdiv)
lmax_div_free_vrtdiv = np.max(np.abs(x.spectogrd(divspec_)))

vrtspec_, divspec_ =  x.getvrtdivspec(ug_psi,vg_psi)
lmax_div_free_psi = np.max(np.abs(x.spectogrd(divspec_)))


if not 'noheader' in sys.argv:
    print("NLATS\tNLONS\t|\tlmax_ug_vrtdiv\tlmax_vg_vrtdiv\tlmax_ug_psi\t|\tlmax_vg_psi\tlmax_div_free_psi\tlmax_div_free_vrtdiv")



print(str(nlats)+"\t"+str(nlons)+"\t"+\
        "|\t"+
        str(lmax_ug_vrtdiv)+"\t"+
        str(lmax_vg_vrtdiv)+"\t"+
        str(lmax_div_free_vrtdiv)+"\t"
        "|\t"+
        str(lmax_ug_psi)+"\t"+
        str(lmax_vg_psi)+"\t"+
        str(lmax_div_free_psi)+"\t"+
        ""
    )


if not 'nooutput' in sys.argv:
    savefileg(ug_psi, "output_ug_psi")
    savefileg(vg_psi, "output_vg_psi")
    savefileg(ug_psi, "output_ug_vrtdiv")
    savefileg(vg_psi, "output_vg_vrtdiv")

