################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for deriving Jovian (atmospheric) data from the      #
# arts-xml-data package and convert it to the common spatial grids             #
# (p_grid), such that they can be applied in radiative transfer calculations.  #
# It is for a 1D atmosphere (for 3D use DemoJupiterAtmo3D.arts instead).       #
#                                                                              #
# It provides following output:                                                #
#   atmosphere_dim    as the WSV                                               #
#   p_grid            as the WSV                                               #
#   lat_grid          as the WSV                                               #
#   lon_grid          as the WSV                                               #
#   z_field           as the WSV                                               #
#   t_field           as the WSV                                               #
#   vmr_field         as the WSV                                               #
#   wind_u/v/w_field  as the WSV                                               #
#   mag_u/v/w_field   as the WSV                                               #
#   abs_species       as the WSV                                               #
#                                                                              #
# The user is supposed to select (atmospheric case, species to include) from   #
# lists. Details of setting rules are given at the place of the settings.      #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/common/createvars.arts                                            #
#                                                                              #
# This template makes use of the following include files                       #
#   includes/mars/atmo_mars.arts                                               #
#   includes/mars/getatmo_mars.arts                                            #
#   includes/common/getgrids_3Dspaced5deg.arts                                 #
#   includes/common/makeatmo3D.arts                                            #
#   includes/mars/getwind_mars.arts                                            #
#   includes/common/makefield3D.arts                                           #
#   includes/jupiter/getmagfield_jupiter.arts"                                 #
#   includes/common/makemagfield3D.arts"                                       #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
# set up name arrays and the like for selections
ws.execute_controlfile("planetary_toolbox/includes/jupiter/atmo_jupiter.arts")
# NOT modify
# prepare the variables for the atmosphere case & species selections
ws.IndexCreate("atmo")
ws.ArrayOfIndexCreate("basespecies")
ws.ArrayOfIndexCreate("h2ospecies")
ws.ArrayOfIndexCreate("nh3species")
ws.ArrayOfIndexCreate("ch4species")
ws.ArrayOfIndexCreate("h2species")
ws.ArrayOfIndexCreate("Necase")
ws.ArrayOfIndexCreate("windcase")
ws.ArrayOfIndexCreate("Bcase")
ws.NumericCreate("latmin")
ws.NumericCreate("latmax")
ws.NumericCreate("lonmin")
ws.NumericCreate("lonmax")
ws.IndexCreate("auxfield_zeropad")
ws.IndexCreate("vmr_zeropad")
ws.IndexCreate("interp_order")
# do NOT modify
# set atmospheric dimensionality to 1D (for 3D use another template!)
ws.AtmosphereSet3D()
# only MODIFY if you know, what you are doing (else the default setting should
#  be fine).
#
# interpolation order for atmospheric data
ws.IndexSet(ws.interp_order, 1)
# assume species-vmr to be zero at pressures/altitudes that are not covered by
#  the species' profile data?
ws.IndexSet(ws.vmr_zeropad, 1)
# as above, but for other data (namely: wind)
ws.IndexSet(ws.auxfield_zeropad, 1)
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Define limits for vertical grid (in terms of pressure)
# ---
# The grid itself is taken from the data (z_field).
# Setting limits to very low and high values will preserve the data grid limit
#  at the respective end.
# Jupiter does not have a surface. Altitude=0km is defined at pressure=1e5(Pa).
#  Jupiter data in the toolbox reaches from p=1e6 at z~-85km to p=1e-9 at
# z~3000km. Pressure @ 100km around 7e2, @ 300km 3e-1, @ 500km 5e-3.
ws.NumericSet(ws.pmin, 0.05)
ws.NumericSet(ws.pmax, 500000.0)
# Define limits for lat/lon grids
# ---
# allowed for latitudes: [-90,90]
# allowed for longitudes:
#    values in [-360,360]
#    lon_max - lon_min <= 360
# note: in case of global coverage (lon_max - lon_min == 360) cyclicity rules
# are applied (e.g. the measurement is allowed to switch over the 360-0 border
# meridian. equality of x and 360+x is recognised).
ws.NumericSet(ws.latmin, -20.0)
ws.NumericSet(ws.latmax, 30.0)
ws.NumericSet(ws.lonmin, 30.0)
ws.NumericSet(ws.lonmax, 50.0)
# ===========================================
# Select the atmospheric scenario to be used
# ---
# (0) mean
# (1) oval*
# *Note: only T, z, and Ne are from oval. all the rest will be from mean. Oval
#  refers to the auroral oval of Jupiter.
ws.IndexSet(ws.atmo, 1)
# ===========================================
# Select the trace gases (and possible sub-scenarios) to be used
# ---
# Basic species
# ---
# refers to species with only one version here. no sub-options/further
#  specifications required
# Select ALL species you like to take into account.
### C2H2, C2H4, C2H6, C3H8, CO, CO2, H2S, HCN, He, PH3
#      0,    1,    2,    3,  4,   5,   6,   7,  8,   9
# H2-CIA-H2-0, H2-CIA-He-0, H2-CIA-CH4-0, CH4-CIA-CH4-0
#         10 ,         11 ,          12 ,           13
ws.ArrayOfIndexSet(ws.basespecies, [4, 5, 6])
# Species with further sub-scenarios
# ---
# only set UP TO ONE for each species (else the species will be included
#  several times, which usually does not make sense).
# EMPTY selection de-selects the whole species.
### H2O: low, high
#           0,   1
# select UP TO ONE
ws.ArrayOfIndexSet(ws.h2ospecies, [0])
### NH3: low, high
#          0,    1
# select UP TO ONE
ws.ArrayOfIndexSet(ws.nh3species, [1])
### electron density: low, medium, high
#                       0,      1,    2
# select UP TO ONE
ws.ArrayOfIndexSet(ws.Necase, [])
# Select species with separate isotopologue profiles available
# ---
# Here it is ok to select more than one entry. General species (aka 'all')
# selects all (remaining, i.e., not yet selected) abs lines of the species. That
# is, general species (i.e., highest index) shall be LAST in selection.
### CH4: CDH3 (212), all (remaining)
#           0      ,   1
# select as many AS YOU WANT, but 'all' index has to be in LAST position
ws.ArrayOfIndexSet(ws.ch4species, [])
### H2: HD (12), all (remaining)
#        0     ,   1
# select as many AS YOU WANT, but 'all' index has to be in LAST position
ws.ArrayOfIndexSet(ws.h2species, [1])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# do NOT modify
# now, let the prepared include files do the actual work:
# (a) read in the raw atmosphere including all selected species
ws.execute_controlfile("planetary_toolbox/includes/jupiter/getatmo_jupiter.arts")
# (b) get the common grids for the atmosphere
ws.execute_controlfile("planetary_toolbox/includes/common/getgrids_3Dspaced5deg.arts")
# (c) do the conversion from raw data with individual grids to the common grids
ws.execute_controlfile("planetary_toolbox/includes/common/makeatmo3D.arts")
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Get non-abs_species data: wind, magnetic field
# ---
# Wind
# Only select ONE element per wind component (else the latter will overwrite the
#  earlier).
# Note, that for Venus no vertical or N-S wind is available, but only E-W wind.
### wind: Galileo,  thermal (for explanation see TN2, p.40)
#               0,        1
# select EXACTLY ONE.
ws.ArrayOfIndexSet(ws.windcase, [1])
# if you want winds (but only then), UNCOMMENT the following 5 lines
# INCLUDE "planetary_toolbox/includes/jupiter/getwind_jupiter.arts"
# Copy( rawfield, wind_u_raw )
# INCLUDE "planetary_toolbox/includes/common/makefield3D.arts"
# Copy( wind_u_field, finalfield )
# IndexSet( abs_f_interp_order, 1 )
# ===========================================
# Magnetic Field
# NOTE: Using the components separately hardly makes sense, we therefore
#  automatically read in all 3 components of the field.
# Only select ONE element (else the latter will overwrite the earlier).
### B: Khurana (for explanation see TN2, p.42)
#            0
# select EXACTLY ONE.
ws.ArrayOfIndexSet(ws.Bcase, [0])
# To SWITCH ON magnetic field, UNCOMMENT the commands below. Else no magnetic
# field will be taken into account.
#
# Copy( lat_true, lat_grid )
# Copy( lon_true, lon_grid )
# INCLUDE "planetary_toolbox/includes/jupiter/getmagfield_jupiter.arts"
# INCLUDE "planetary_toolbox/includes/common/makemagfield.arts"
# ===========================================
