#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Manfred Brath

"""
import numpy as np

import pyarts

ws=pyarts.Workspace()

#Number of frequencies
N_f = 4

#set test frequency grid
ws.f_grid = np.linspace(0,2,N_f) # Hz

#set quadrature test weights
#Here we use a simple sin function to make dummy quadrature weights
#The quadrature weights are normalized to the range of the frequency grid
quadrature_weights= np.linspace(1,np.pi,N_f)
quadrature_weights /= np.sum(quadrature_weights)
quadrature_weights *= (ws.f_grid.value[-1]-ws.f_grid.value[0])

#create  test spectral irradiance field
spectral_irradiance_field=np.zeros((N_f,1,1,1,2))
spectral_irradiance_field[:,0,0,0,0]=-np.linspace(0,1,N_f)
spectral_irradiance_field[:,0,0,0,1]=np.linspace(0,1,N_f)
ws.spectral_irradiance_field=spectral_irradiance_field

#create test spectral radiance field
spectral_radiance_field=np.zeros((N_f,1,1,1,2,1,1))
spectral_radiance_field[:,0,0,0,0,0,0]=np.linspace(0,1,N_f)
spectral_radiance_field[:,0,0,0,1,0,0]=np.linspace(0,1,N_f)
ws.spectral_radiance_field=spectral_radiance_field

#reference calculation
irradiance_field_ref=np.sum(spectral_irradiance_field*quadrature_weights[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis],0)

#for radiance field reference, we have to remove the last dimension of dummy_field_rad,
#because in the ARTS function the last dimension (stokes dimension) is also removed.
radiance_field_ref=np.sum(spectral_radiance_field[:,:,:,:,:,:,0]*quadrature_weights[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis,np.newaxis],0)

#do the integration for the spectral irradiance field
ws.RadiationFieldSpectralIntegrate(ws.irradiance_field, ws.f_grid, ws.spectral_irradiance_field, quadrature_weights)

#do the integration for the spectral radiance field
ws.RadiationFieldSpectralIntegrate(ws.radiance_field, ws.f_grid, ws.spectral_radiance_field, quadrature_weights)


#compare
ws.Compare(ws.irradiance_field, irradiance_field_ref, 1e-6)
ws.Compare(ws.radiance_field, radiance_field_ref, 1e-6)
