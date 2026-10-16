#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:44:35 2026

@author: u242031
"""
import numpy as np
import pyarts3.arts as paa3


freq=np.linspace(10e9,120e9,4)
refl_data=np.random.rand(len(freq))
refl1d=paa3.SortedGriddedField1(data=refl_data, grids=[freq])


lat=np.linspace(-90,90,3)
lon=np.linspace(-180,180,4)
freq2d=np.linspace(10e9,100e9,5)
refl2d_data=np.random.rand(len(lat),len(lon),len(freq2d))
refl2d=paa3.SortedGriddedField3(data=refl2d_data, grids=[lat,lon,freq2d])


Lambertian=paa3.LambertianSurfaceScatterer(refl1d)
LambertianField=paa3.LambertianSurfaceScattererField(refl2d)

surface_models=paa3.MapOfSurfaceScatteringModel()
surface_models.add('type1',Lambertian)
surface_models['types2']=LambertianField


print(surface_models)