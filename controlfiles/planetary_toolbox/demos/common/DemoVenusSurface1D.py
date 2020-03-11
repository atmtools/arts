################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template INCLUDE file for setting Martian surface properties for   #
# following conditions:                                                        #
#  - 1D calculations only                                                      #
#  - non-scattering calculations (?)                                           #
#                                                                              #
# It offers setups for the following parameters and settings                   #
#  - PART A: Surface Altitude (Topography)                                     #
#    - CASE 1: Surface at a given altitude                                     #
#    - CASE 2: Surface at altitude of a given atmospheric level                #
#   (- CASE 3: Surface altitude from read-in topography file (NOT available))  #
#  - PART B: Surface temperature + Surface reflectivity                        #
#    - CASE 1: Surface of a given temperature                                  #
#    - CASE 2: Surface with temperature from atmospheric temperature at        #
#               surface altitude                                               #
#    - CASE 3: Surface with (atmospheric scenario dependent) surface           #
#               temperature (global/regional) field from file                  #
#    where each has the subCASEs                                               #
#      - subCASE a: Blackbody surface                                          #
#      - subCASE b: Lambertian surface of given reflectivity                   #
#      - subCASE c: Specular surface of given reflectivity                     #
#      - subCASE d: Specular surface of given surface refractive index         #
#      - subCASE e: Specular surface with surface refractive index             #
#                    (global/regional) field from file                         #
#                                                                              #
# Files to be included before this file:                                       #
#   ../controlfiles/general/agendas.arts  for agenda predefinitions            #
#   includes/common/createvars.arts       for creating required variable       #
#                                          containers                          #
#   demos/common/DemoVenusAtmo1D.arts     for atmospheric scenario selection   #
#                                                                              #
# It requires the following input:                                             #
#   p/lat/lon_grid         as the WSVs                                         #
#   t_field                as the WSV                                          #
#                                                                              #
# It provides following output (in parentheses: restricted to certain CASEs):  #
#   z_surface              as the WSV                                          #
#   surface_rtprop_agenda  as the WSA                                          #
#   (t_surface)            as the WSV                                          #
#   (surface_skin_t)       as the WSV                                          #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
ws.IndexCreate("surface_level")
ws.NumericCreate("surface_altitude")
ws.NumericCreate("reflectivity")
ws.NumericCreate("surface_refractive_index_real")
ws.NumericCreate("surface_refractive_index_imag")
ws.StringCreate("surface_refractive_index_file")
ws.StringCreate("surface_temperature_file")
ws.StringCreate("surface_altitude_file")
ws.GriddedField2Create("surface_temperature_field")
ws.GriddedField2Create("surface_altitude_field")
ws.GriddedField5Create("surface_refractive_index_field")
ws.MatrixSetConstant(ws.t_surface, 1, 1, -1.0)
# only MODIFY if you know, what you are doing (else the default setting should
#  be fine).
# File containing surface altitude (global/regional) field
# Used for cases: A-3
#
# There is NO altitude data available in the toolbox for Venus. The code can be
#  used as template, though, if you DO have data from other sources.
# StringSet( surface_altitude_file, "YOUR-DATA-LOCATION/z_surface.xml.gz" )
# File containing surface refractive index (global/regional) field
# Used for cases: B-1e, B-2e, B-3e
#
# The only available Venus surface refractive index data in toolbox. So, change
#  only if you really want to use (and have available) other data.
ws.StringSet(
    ws.surface_refractive_index_file,
    "planets/Venus/MPS/Venus.surface_complex_refr_index_field.xml",
)
# File containing surface temperature (global/regional) field
# Used for cases: B-3
#
# The only available Venus surface temperature data in toolbox. So, change only
#  if you really want to use (and have available) other data.
#
ws.StringSet(ws.surface_temperature_file, "planets/Venus/MPS/Venus.t_surface.xml")
################################################################################
# START USER SETTINGS - Modify settings according to you wishes                #
################################################################################
# NOTE: For each case (start marked by 'CASE') you are supposed to UNCOMMENT all
#  command lines till the start of the next (sub)CASE or PART description/
#  definition.
#  Furthermore, you are also supposed to SET the values of certain parameters,
#  marked by 'SET' in the comment line before the respective command lines.
# Define geographic position of the 1D case
# ---
# NOTE: only required for some cases (when data is extracted from global
#  fields), namely for: A-3, B-1e, B-2e, B-3
#  If you do not select any of these cases, the values set here do not matter.
# SET one(!) number (-90...90)
ws.VectorSet(ws.lat_true, array([20.0]))
# SET one(!) number (-360...360)
ws.VectorSet(ws.lon_true, array([-50.0]))
# -------
# PART A: Surface Altitude (Topography)
# -------
# CASE 1) Surface at a given altitude (in m)
# ---
# NOTE: Can NOT be lower than the lowest atmospheric level of the scenario.
#  For Venus, this is at max 0m, but might be higher depending on choice of pmax.
#
# SET number (0m <= surface_altitude!):
ws.NumericSet(ws.surface_altitude, 0.0)
ws.MatrixSetConstant(ws.z_surface, 1, 1, ws.surface_altitude)
# CASE 2) Surface at altitude of a given atmospheric level
# ---
# We extract the respective level from the atmospheric z_field by its index
#  surface_level. surface_level=0 indicates the lowest atmospheric level given
#  by the atmospheric scenario.
# NOTE: surface_level=0 is the intended application of this case (A-2).
#
# SET index (0 <= surface_level!):
# IndexSet( surface_level, 0 )
# Extract( z_surface, z_field, surface_level )
# CASE 3) Surface altitude from read-in topography file
# ---
# NOTE: NOT available for Venus as no surface altitude data available in
#  toolbox. You might use the below code as template, though in case you DO have
#  data from other sources.
#
# ReadXML( surface_altitude_field, surface_altitude_file )
# GriddedFieldLatLonRegrid( out=surface_altitude_field,
#                          in=surface_altitude_field )
# FieldFromGriddedField( out=z_surface, in=surface_altitude_field )
# -------
# PART B: Surface temperature + Surface reflectivity
# -------
# Note: As these two are governed by the surface_rtprop_agenda, it is not
#  possible to prepare their settings independently of each other. Hence, they
#  are handled together in PART B.
# CASE 1) Surface of a given temperature
# ---
# NOTE: This setting of surface_skin_t below is for all subCASEs 1a-1e below,
#  i.e., UNCOMMENT and set value if any of the CASE 1 subcases is used.
#
# SET number:
# NumericSet( surface_skin_t, 180. )
# subCASE 1a) Blackbody surface
# ---
# Copy( surface_rtprop_agenda, surface_rtprop_agenda__Blackbody_SurfTFix )
# subCASE 1b) Lambertian surface of given reflectivity
# ---
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# SET index (number of rays to consider in reflection calulation, for details
#  see WSV description):
# IndexSet( lambertian_nza, 18 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__lambertian_ReflFix_SurfTFix )
# subCASE 1c) Specular surface of given reflectivity
# ---
# NOTE: only for scalar case: stokes_dim=1.
#  For stokes_dim>1 specify surface_reflectivity (Tensor3; see WSV decription)
#  and use surface_rtprop_agenda__Specular_Pol_ReflFix_SurfTFix.
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFix )
# subCASE 1d) Specular surface of given surface refractive index
# ---
#
# SET number (1 <= surface_refractive_index_real):
# NumericSet( surface_refractive_index_real, 2.0 )
# SET number:
# NumericSet( surface_refractive_index_imag, 0.0 )
# complex_refr_indexConstant( surface_complex_refr_index,
#                            surface_refractive_index_real,
#                            surface_refractive_index_imag )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__Specular_Pol_RefracFix_SurfTFix )
# subCASE 1e) Specular surface with surface refractive index (global/regional)
#             field from file
# ---
# NOTE: Requires lat_true and lon_true to be set. Requires
#  surface_refractive_index_file to be set.
#
# ReadXML( out=surface_refractive_index_field,
#         filename=surface_refractive_index_file )
# AgendaSet( surface_rtprop_agenda ){
#  specular_losCalc
#  surface_complex_refr_indexFromGriddedField5(
#      complex_refr_index_field=surface_refractive_index_field )
#  surfaceFlatRefractiveIndex
# }
# CASE 2) Surface with temperature from atmospheric temperature at surface
#         altitude
# ---
# NOTE: Temperature is derived in-situ everytime surface_rtprop_agenda is
#  executed. Hence, we do not need a global setting for this case, but include
#  the appropriate method within surface_rtprop_agenda definition of each
#  subCASE.
# subCASE 2a) Blackbody surface
# ---
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__Blackbody_SurfTFromt_field )
# subCASE 2b) Lambertian surface of given reflectivity
# ---
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# SET index (number of rays to consider in reflection calulation, for details
#  see WSV description):
# IndexSet( lambertian_nza, 18 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field )
# subCASE 2c) Specular surface of given reflectivity
# ---
# NOTE: only for scalar case: stokes_dim=1.
#  For stokes_dim>1 specify surface_reflectivity (Tensor3; see WSV decription)
#  and use surface_rtprop_agenda__Specular_Pol_ReflFix_SurfTFix.
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field )
# subCASE 2d) Specular surface of given surface refractive index
# ---
#
# SET number (1 <= surface_refractive_index_real):
# NumericSet( surface_refractive_index_real, 2.0 )
# SET number:
# NumericSet( surface_refractive_index_imag, 0.0 )
# complex_refr_indexConstant( surface_complex_refr_index,
#                            surface_refractive_index_real,
#                            surface_refractive_index_imag )
# AgendaSet( surface_rtprop_agenda ){
#  specular_losCalc
#  InterpAtmFieldToPosition( out=surface_skin_t, field=t_field )
#  surfaceFlatRefractiveIndex
# }
# subCASE 2e) Specular surface with surface refractive index (global/regional)
#             field from file
# ---
# NOTE: Requires lat_true and lon_true to be set. Requires
#  surface_refractive_index_file to be set.
#
# ReadXML( out=surface_refractive_index_field,
#         filename=surface_refractive_index_file )
# AgendaSet( surface_rtprop_agenda ){
#  specular_losCalc
#  surface_complex_refr_indexFromGriddedField5(
#      complex_refr_index_field=surface_refractive_index_field )
#  InterpAtmFieldToPosition( out=surface_skin_t, field=t_field )
#  surfaceFlatRefractiveIndex
# }
# CASE 3) Surface with (atmospheric scenario dependent) surface temperature
#         (global/regional) field from file
# ---
# NOTE: The field reading and preparation settings below is for all
#  subCASEs 3a-3e, i.e., UNCOMMENT if any of the CASE 3 subcases is used.
# NOTE: Requires lat_true and lon_true to be set. Requires
#  surface_temperature_file to be set.
#
ws.StringSet(ws.infostr, "Surface temperature data taken from: ")
ws.Append(ws.infostr, ws.surface_temperature_file)
ws.Print(ws.infostr)
ws.ReadXML(ws.surface_temperature_field, ws.surface_temperature_file)
ws.GriddedFieldLatLonRegrid(
    ws.surface_temperature_field,
    ws.lat_true,
    ws.lon_true,
    ws.surface_temperature_field,
    1,
)
ws.FieldFromGriddedField(
    ws.t_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.surface_temperature_field
)
# subCASE 3a) Blackbody surface
# ---
# AgendaSet( surface_rtprop_agenda ){
#  Ignore( rtp_los )
#  InterpSurfaceFieldToPosition( out=surface_skin_t, field=t_surface )
#  surfaceBlackbody
# }
# subCASE 3b) Lambertian surface of given reflectivity
# ---
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# SET index (number of rays to consider in reflection calulation, for details
#  see WSV description):
# IndexSet( lambertian_nza, 18 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# AgendaSet( surface_rtprop_agenda ){
#  InterpSurfaceFieldToPosition( out=surface_skin_t, field=t_surface )
#  surfaceLambertianSimple( za_pos=0.5 )
# }
# subCASE 3c) Specular surface of given reflectivity
# ---
# NOTE: only for scalar case: stokes_dim=1.
#  For stokes_dim>1 specify surface_reflectivity (Tensor3; see WSV decription)
#  and use surface_rtprop_agenda__Specular_Pol_ReflFix_SurfTFix.
#
# SET number (0 <= reflectivity <= 1; emissivity = 1-reflectivity):
# NumericSet( reflectivity, 0.15 )
# VectorSetConstant( surface_scalar_reflectivity, 1, reflectivity )
# Copy( surface_rtprop_agenda,
#      surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface )
# subCASE 3d) Specular surface of given surface refractive index
# ---
#
# SET number (1 <= surface_refractive_index_real):
# NumericSet( surface_refractive_index_real, 2.0 )
# SET number:
# NumericSet( surface_refractive_index_imag, 0.0 )
# complex_refr_indexConstant( surface_complex_refr_index,
#                            surface_refractive_index_real,
#                            surface_refractive_index_imag )
# AgendaSet( surface_rtprop_agenda ){
#  specular_losCalc
#  InterpSurfaceFieldToPosition( out=surface_skin_t, field=t_surface )
#  surfaceFlatRefractiveIndex
# }
# subCASE 3e) Specular surface with surface refractive index (global/regional)
#             field from file
# ---
# NOTE: Requires lat_true and lon_true to be set. Requires
#  surface_refractive_index_file to be set.
#
ws.ReadXML(
    out=ws.surface_refractive_index_field, filename=ws.surface_refractive_index_file
)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.surface_complex_refr_indexFromGriddedField5(
        complex_refr_index_field=ws.surface_refractive_index_field
    )
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# end of Arts2 tag
