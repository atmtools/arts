################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for deriving Venus (atmospheric) data from the       #
# arts-xml-data package and convert it to the common spatial grids             #
# (p_grid), such that they can be applied in radiative transfer calculations.  #
# It is for a 1D atmosphere (for 3D use DemoVenusAtmo3D.arts instead).         #
#                                                                              #
# It provides following output:                                                #
#   atmosphere_dim    as the WSV                                               #
#   p_grid            as the WSV                                               #
#   z_field           as the WSV                                               #
#   t_field           as the WSV                                               #
#   vmr_field         as the WSV                                               #
#   wind_u/v/w_field  as the WSV                                               #
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
#   includes/venus/atmo_venus.arts                                             #
#   includes/venus/getatmo_venus.arts                                          #
#   includes/common/getgrids_1D.arts                                           #
#   includes/common/makeatmo1D.arts                                            #
#   includes/venus/getwind_venus.arts                                          #
#   includes/common/makefield1D.arts                                           #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
# set up name arrays and the like for selections
ws.execute_controlfile("planetary_toolbox/includes/venus/atmo_venus.arts")
# do NOT modify
# prepare the variables for the atmosphere case & species selections
ws.IndexCreate("atmo")
ws.ArrayOfIndexCreate("basespecies")
ws.ArrayOfIndexCreate("h2ospecies")
ws.ArrayOfIndexCreate("hdospecies")
ws.ArrayOfIndexCreate("Necase")
ws.ArrayOfIndexCreate("NSwind")
ws.ArrayOfIndexCreate("EWwind")
ws.IndexCreate("auxfield_zeropad")
ws.IndexCreate("vmr_zeropad")
ws.IndexCreate("interp_order")
# do NOT modify
# set atmospheric dimensionality to 1D (for 3D use another template!)
ws.AtmosphereSet1D()
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
# For Venus, surface pressure is about 1e7, pressure @ 100km around 1e1, and
#  pressure @ 500km (highest altitude on Venus, where we have data) around
#  5e-12.
ws.NumericSet(ws.pmin, 10.0)
ws.NumericSet(ws.pmax, 1e99)
# ===========================================
# Select the atmospheric scenario to be used
# ---
# (0) night               (from SpicaV)
# (1) night, cold         (from SpicaV)
# (2) night               (from VIRA)
# (3) day                 (from VIRA)
# (4) day, high latitudes (from VIRA)
ws.IndexSet(ws.atmo, 3)
# ===========================================
# Select the trace gases (and possible sub-scenarios) to be used
# ---
# Basic species
# ---
# refers to species with only one version here. no sub-options/further
#  specifications required
#
# NOTE: if you select CO2-CO2 CIA here, you have to use an appropriate set for
#  the abs_xsec_agenda.
#
# Select ALL species you like to take into account.
### CO, CO2, CO2-CO2 CIA, CO2-CO2 PWR, H2SO4, HCl, HF, N2, NO, NO2, O, O2,
#    0 , 1 ,      2     ,     3      ,   4  ,  5 ,  6 , 7,  8 , 9, 10, 11,
# O3, OCS, SO, SO2
# 12,  13, 14,  15
ws.ArrayOfIndexSet(ws.basespecies, [0, 1, 2, 4, 13, 14, 15])
# Species with sub-scenarios
# ---
# only set UP TO ONE for each species (else the species will be included
#  several times, which usually does not make sense).
# EMPTY selection de-selects the whole species.
### H2O: low, medium, high
#         0 ,   1   ,  2
# select UP TO ONE
ws.ArrayOfIndexSet(ws.h2ospecies, [1])
# HDO
# NOTE: ONLY make a (non-empty) selection here, if you REALLY want to take
#  HDO separately (i.e. with a separate, specific profile) into account.
# If EMPTY selection here, HDO is routinely covered by the H2O (general) species
#  selectable above using the corresponding H2O profile scaled by isotopologue
#  abundance.
# low/medium/high profiles are related to the low/medium/high profiles of H2O,
#  hence are suggested to be used consistently with the H2O profile (i.e., we
#  discourage use of, e.g. H2O-low with HDO-medium or -high, and so on). The
#  uncorrected profile originates from independent measurements of HDO. (for
#  details see TN2 p.21)
### low, medium, high, original (uncorrected)
#    0 ,   1,     2  ,    3
# select UP TO ONE
ws.ArrayOfIndexSet(ws.hdospecies, [3])
# Electron density: depending on solar zenith angle.
# NOTE that SZA<90 are only available for day data, SZA>90 for night data
#  (non-matching selection will lead to a runtime error!). That is, COORDINATE
#  selection here with daytime according to your selection of ATMO.
#       day       day       day       day       day      night       night
### SZA0-30, SZA30-50, SZA50-70, SZA70-80, SZA80-90, SZA90-100, SZA100-120
#         0,        1,        2,        3,        4,         5,          6
# select UP TO ONE
ws.ArrayOfIndexSet(ws.Necase, [1])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# do NOT modify
# now, let the prepared include files do the actual work:
# (a) read in the raw atmosphere including all selected species
ws.execute_controlfile("planetary_toolbox/includes/venus/getatmo_venus.arts")
# (b) get the common grids for the atmosphere
ws.execute_controlfile("planetary_toolbox/includes/common/getgrids_1D.arts")
# (c) do the conversion from raw data with individual grids to the common grids
ws.execute_controlfile("planetary_toolbox/includes/common/makeatmo1D.arts")
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Get non-abs_species data: wind
# ---
# NOTE: If you want N-S or E-W wind, you have to use a 3D atmosphere. Use the
#  other, 3D atmosphere, template for that!
# Note: For Venus, no vertical wind data is available (only N-S and E-W wind).
#  Hence, no considering of wind in Venus 1D case (it won't have any effect
#  here) and you HAVE to use the 3D template, if you want wind.
################################################################################
# END USER SETTINGS                                                            #
################################################################################
