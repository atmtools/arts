#
# 2014-03-17 Jana Mendrok
#
################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs. However, MAKE A COPY of this template and store in a       #
# separate folder BEFORE starting to modify.                                   #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for "full" radiative transfer simulations, i.e.,     #
# providing synthetic sensor measurements. The demo case covers following      #
# conditions/capabilities:                                                     #
#  - passive sensor                                                            #
#  - 1D calculations only                                                      #
#  - monochromatic pencilbeam simulation (no FOV/antenna pattern convolution,  #
#     no bandwidth/instrumental line shape/... convolution)                    #
#  - non-particle (i.e., "clearsky") calculations                              #
#  - no magnetic field effects (Zeeman, Faraday)                               #
#  - for certain conditions (shallow viewing angles), refraction calculations  #
#    might fail                                                                #
#                                                                              #
#  + scalar to full polarisation RT possible (limitations might apply in       #
#     combination with, e.g., certain surface property settings)               #
#  + multiple viewing directions (allowing their specification by viewing      #
#     angles and tangent altitudes in parallel) from fixed platform altitude   #
#  + selection of atmospheric scenario and abs_species from toolbox data       #
#     (handled by INCLUDE file)                                                #
#  + consideration of wind Doppler shifts (by switching on winds in            #
#     atmospheric scenario INCLUDE file)                                       #
#  + selection of absorption calculation approach (on-the-fly, lookup table    #
#     calculation, pre-calculated lookup table from file)                      #
#  + consideration of refraction selectable                                    #
#  + selection of surface properties (handled by INCLUDE file)                 #
#                                                                              #
#                                                                              #
# The user is supposed to set certain parameters or select from a choice of    #
# settings. Details of setting rules are given at the place of the settings.   #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
#                                                                              #
# The file provides following OUTPUT (written to file):                        #
#   iy         as the WSV                                                      #
#               radiance; units selectable                                     #
#               one file per pencilbeam (i.e., per pointing direction =        #
#                 per viewing angle/tangent altitude)                          #
#               name of output file is DEMONAME.iy.LOS-ID.xml, where DEMONAME  #
#                 is the name of this controlfile (ripped off the .arts        #
#                 extension; i.e., DEMONAME changes/adapts when you save the   #
#                 template under a different name) and LOS-ID is the running   #
#                 index of the pointing direction. Will be located in the      #
#                 directory, from where the controlfile is executed.           #
#   iy_aux     as the WSV                                                      #
#               auxiliary output parameters (particularly of along-the-path    #
#               type), selectable variety                                      #
#               one file per pencilbeam (i.e., viewing angle)                  #
#               all different output parameters per pencilbeam in one file     #
#               name of output file is DEMONAME.iy_aux.LOS-ID.xml              #
#                                                                              #
#   Further variables can be output (actually, ANY workspace variable) to file #
#   using WriteXML( in=VARIABLENAME ). The demo is setup to also provide:      #
#     f_grid         as the WSV                                                #
#     abs_species    as the WSV                                                #
#                                                                              #
#   Files tmp*.xml* are intermediate data needed within the calculation. They  #
#   remain in your execution directory afterwards, but can safely be deleted.  #
#                                                                              #
#                                                                              #
# This template is set up to make use of the following include files           #
# (internally they might use further includes. See their respective headers    #
# for details):                                                                #
#   demos/common/DemoVenusAtmo1D.arts                                          #
#   demos/common/DemoVenusSurface1D.arts                                       #
#   includes/common/createvars.arts                                            #
#   includes/common/makegeometry1D_unrefracted.arts                            #
#   includes/common/makegeometry1D_refracted-air_effectivetanh.arts            #
#   includes/common/makegeometry1D_refracted-air_geometrictanh.arts            #
#   includes/common/makegeometry1D_refracted-air+electrons_effectivetanh.arts  #
#   includes/common/makegeometry1D_refracted-air+electrons_geometrictanh.arts  #
#   includes/common/use-absOnTheFly.arts                                       #
#   includes/common/make-and-use-absLUT.arts                                   #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_venus.arts")
# do NOT modify
# prepare some variable containers we will need for later INCLUDE file executions.
# this INCLUDE has to appear at the beginning. DO NOT include it more than ONCE
# in a full ARTS run.
ws.execute_controlfile("planetary_toolbox/includes/common/createvars.arts")
# do NOT modify
ws.NumericCreate("obsh")
ws.VectorCreate("tanh")
ws.VectorCreate("viewang")
# No Jacobian Matrix in the basic test
ws.jacobianOff()
ws.Touch(ws.nlte_vibrational_energies)
################################################################################
# START USER SETTINGS - Modify settings according to you wishes                #
################################################################################
# unit of intensity output( "PlanckBT" or "RJBT" or "1" )
ws.StringSet(ws.iy_unit, "PlanckBT")
# monochromatic frequency grid
ws.VectorLinSpace(ws.f_grid, 595000000000.0, 600000000000.0, 10000000.0)
ws.WriteXML(ws.output_file_format, ws.f_grid, "", 0)
# number of Stokes parameters to consider ( 1: scalarRT, 4: fully polarised )
ws.IndexSet(ws.stokes_dim, 1)
# ok to use extrapolation of partition functions (i.e. fitting range of
#  partition function paameterization does not cover the actual atmospheric
#  temperatures)?
# ( 0: not ok, 1: ok )
ws.IndexSet(ws.bad_partition_functions_ok, 1)
# ---
# Definition of atmosphere
# ---
# which scenario, which absorption species, ...
#  - use a COPY of DemoVenusAtmo1D.arts adapted to your needs/wishes
#  - provides as output: atmosphere_dim (1D), abs_species , p_grid, z/t/vmr_field,
#     wind_w_field (if switched on in adaptation of DemoVenusAtmo1D.arts)
#####
ws.execute_controlfile("planetary_toolbox/demos/common/DemoVenusAtmo1D.arts")
ws.WriteXML(ws.output_file_format, ws.abs_species, "", 0)
# ---
# Sensor parameters: viewing geometry
# ---
# platform altitude (will be identical for all viewing directions)
ws.NumericSet(ws.obsh, 450000.0)
# pointing direction
# note: You can specify line of sight (LOS) in terms of zenith angle AND tangent
#  altitude. Both can be used within one ARTS run. For that, tangent altitude
#  LOS, after being converted to LOS zenith angles, are appended to zenith angle
#  LOS, and both together are stored in WSV allzang (viewang directions are
#  first! if you want allzang as output, use WriteXML to write it to file). That
#  is, RT calculations for n_LOS = n_viewang + n_tanh pencilbeam LOS directions
#  will be performed.
#
# (a) zenith viewing angles (e.g., for slantlooking, but working for ALL geometries)
# - can be empty
ws.VectorSet(ws.viewang, np.array([180.0, 150.0, 120.0]))
# (b) tangent altitudes (in m) for limb views
# - can be empty
# - depending on your choice of type of geometry/refraction below, specified
#    tangent altitudes will either be "theoretical" tangent altitudes (as would
#    be without refraction) or "true" tangent altitudes (the true tangent
#    altitude of a refracted ray)
# - specified tanh are allowed to be negative. In case of theoretical tangent
#    altitudes (but not for true ones!), they can even be "false" tangent
#    altitudes, i.e., pointing BELOW the actual surface level.
ws.VectorSet(
    ws.tanh, np.array([0.0, 4000.0, 10000.0, 12000.0, 20000.0, 40000.0, 60000.0])
)
# ---
# Type of geometry/refraction
# ---
#  Prepares the (1D) viewing angles for the RT calculation from the above
#   specified viewang and tanh. Select proper INCLUDE for the type/level of
#   refraction that shall be considered.
#  NOTE: uncomment the ONE(!) you want. If more than one uncommented, the run
#   will crash.
# ---
# 1) no refraction at all
# INCLUDE "planetary_toolbox/includes/common/makegeometry1D_unrefracted.arts"
# 2) refraction (by "air" only)
# 2-1) tanh => true tangent altitudes of refracted rays
ws.execute_controlfile(
    "planetary_toolbox/includes/common/makegeometry1D_refracted-air_truetanh.arts"
)
# 2-2) tanh => theoretical tangent altitudes of unrefracted rays
#  (true tangent altitudes will be lower)
# INCLUDE "planetary_toolbox/includes/common/makegeometry1D_refracted-air_geometrictanh.arts"
# 3) refraction (by "air" and free electrons)
#  These require abs=species to contain "free_electrons".
# 3-1) tanh => true tangent altitudes of refracted rays
# INCLUDE "planetary_toolbox/includes/common/makegeometry1D_refracted-air+electrons_truetanh.arts"
# 3-2) tangent altitudes = theoretical tangent altitudes of unrefracted rays
#  (true ones will be lower)
# INCLUDE "planetary_toolbox/includes/common/makegeometry1D_refracted-air+electrons_geometrictanh.arts"
# ---
# Absorption calculation approach
# ---
#  NOTE: uncomment the ONE(!) you want. If more than one uncommented, the LATEST
#   RULES, but you might waste time with in-the-end-not-applied calculations.
# ---
# 1) use on-the-fly absorption
#  (time consuming! but recommended when doing calculations with winds or other
#  Doppler shifts.)
# INCLUDE "planetary_toolbox/includes/common/use-absOnTheFly.arts"
# 2) calculate and use absorption lookup table (generally recommended)
ws.execute_controlfile("planetary_toolbox/includes/common/make-and-use-absLUT.arts")
# 2-1) if you want to keep the absLUT data for later use or analysis, uncomment:
# WriteXML( in=abs_lookup )
# 3) use a previously calculated absorption lookup table
#  if you want to use an absLUT created earlier, uncomment the following 3 lines
#  (you might also have to specify the name. check online-doc of ReadXML)
# Copy( propmat_clearsky_agenda, propmat_clearsky_agenda__LookUpTable )
# ReadXML( out=abs_lookup )
# abs_lookupAdapt
# ---
# Surface settings
# ---
# definition of surface (altitude, reflection model, temperature)
#  - use a copy of DemoVenusSurface1D.arts adapted to your needs/wishes
#  - provides as output: z_surface, surface_rtprop_agenda (additionally:
#     surface_skin_t (cases B-1: manually set (scalar) surface temperature) or
#     t_surface (cases B-3: surface temperature field from file))
#####
ws.execute_controlfile("planetary_toolbox/demos/common/DemoVenusSurface1D.arts")
# ---
# Define (auxiliary) data output
# ---
# Uncomment all parameters you want as auxiliary output (i.e., in addition to
#  total radiance/brigthness temperature). For meaning of each paramters see
#  online-doc of the WSM selected for iy_main_agenda (here: iyEmissionStandard).
# NOTE: Last element NOT to be followed by comma.
# NOTE: Only use "Absorption, species X" up to the number of entries in
#  abs_species (clearsky calculations in Venus have at maximum 19 abs_species
#  entries, i.e. highest valid index is 18).
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Radiative background", "Optical depth"])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# setting agendas needed for RT calc (there are alternative settings, though)
#####
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# no scattering, no jacobian
#####
ws.jacobianOff()
ws.cloudboxOff()
# the checks necessary for full RT calc
#####
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
# and the (clearsky) RT calc
#####
ws.NumericCreate("za")
ws.AgendaCreate("forloop_agenda_angles")


@arts_agenda
def forloop_agenda_angles(ws):
    ws.Extract(ws.za, ws.allzang, ws.forloop_index)
    ws.rte_losSet(za=ws.za, aa=ws.za)
    ws.Print(ws.rte_los, 0)
    ws.iyCalc()
    ws.WriteXMLIndexed(ws.output_file_format, ws.forloop_index, ws.iy, "", 0)
    ws.WriteXMLIndexed(ws.output_file_format, ws.forloop_index, ws.iy_aux, "", 0)


ws.forloop_agenda_angles = forloop_agenda_angles

ws.IndexCreate("nangles")
ws.nelemGet(ws.nangles, ws.allzang)
ws.IndexStepDown(ws.nangles, ws.nangles)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_angles)
ws.ForLoop(ws.forloop_agenda, 0, ws.nangles, 1)
