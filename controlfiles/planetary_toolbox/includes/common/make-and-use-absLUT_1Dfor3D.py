################################################################################
#                                                                              #
# This is a (plug&play-type) include file. The USER is NOT supposed to MODIFY  #
# it, but choose another include file to suit his/her needs.                   #
#                                                                              #
################################################################################
#                                                                              #
# This INCLUDE file is for                                                     #
#  - calculating an absorption lookup table                                    #
#  - including CIA continua (works also, if no CIA species present)            #
#  - NOTE: temporary workaround for (homogeneous) 3D atmosphere as abs_lookup  #
#     from true 3D currently has some issues                                   #
#                                                                              #
# It performs the following actions:                                           #
#  - sets propmat_clearsky_agenda: use absorption lookup table                 #
#  - sets abs_xsec_agenda: include CIA and allow temperature extrapolation of  #
#     CIA data (works fine even in absence of CIA species)                     #
#  - read CIA input data from toolbox data package                             #
#  - read required spectroscopic line files from toolbox data package          #
#  - (re)calculates 1D atmospheric fields                                      #
#  - calculate the lookup table                                                #
#                                                                              #
# It requires the following input:                                             #
#   abs_species               as the WSV                                       #
#   p_grid                    as the WSV                                       #
#   t/z/vmr_field_raw         as the WSV                                       #
#   abs_p/t/nls_interp_order  as the WSVs                                      #
#   bad_partition_functions_ok  (Index)   Partition function handling flag     #
#                                                                              #
# It provides following output (in parentheses: side products):                #
#   abs_lookup                as the WSV                                       #
#   propmat_clearsky_agenda   as the WSA                                       #
#   abs_xsec_agenda           as the WSA                                       #
#   (abs_lines)               as the WSV                                       #
#   (abs_lines_per_species)   as the WSV                                       #
#   (abs_cia_data)            as the WSV                                       #
#   (abs_p/t/t_pert/vmrs/nls/nls_pert) as the WSVs                             #
#   (abs_xsec_agenda_checked) as the WSV                                       #
#                                                                              #
# The file shall NOT be modified by the USER.                                  #
# However, it is ok to adapt the catalogue type & location to use (e.g., use   #
#  your own copy of HITRAN). You might also adapt the low and high frequency   #
#  limits considered IF you KNOW what you are doing.                           #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# for abs.coeffs prepare & use a lookup table
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# for CIA allow temperature extrapolation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__withCIAextraT)
# for the selected abs_species and atmospheric scenario, calculate the
#  absorption lookup table
#####
ws.ReadXML(ws.abs_cia_data, "spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml.gz")
ws.abs_linesReadFromSplitArtscat(
    basename="spectroscopy/Perrin/", fmin=0.0, fmax=4000000000000.0
)
# abs_linesReadFromHitran( filename="your-HITRAN-file", fmin=0, fmax=4e12 )
ws.abs_lines_per_speciesCreateFromLines()
ws.VectorCreate("lat_1D")
ws.VectorCreate("lon_1D")
ws.IndexCreate("atmdim_1D")
ws.Tensor3Create("t_1D")
ws.Tensor3Create("z_1D")
ws.Tensor4Create("vmr_1D")
ws.Tensor3Create("fdummy")
ws.Touch(ws.fdummy)
ws.AtmosphereSet1D(ws.atmdim_1D, ws.lat_1D, ws.lon_1D)
ws.AtmFieldsCalc(
    t_field=ws.t_1D,
    z_field=ws.z_1D,
    vmr_field=ws.vmr_1D,
    lat_grid=ws.lat_1D,
    lon_grid=ws.lon_1D,
    atmosphere_dim=ws.atmdim_1D,
    interp_order=ws.interp_order,
    vmr_zeropadding=ws.vmr_zeropad,
)
ws.atmfields_checkedCalc(
    atmosphere_dim=ws.atmdim_1D,
    lat_grid=ws.lat_1D,
    lon_grid=ws.lon_1D,
    t_field=ws.t_1D,
    vmr_field=ws.vmr_1D,
    wind_u_field=ws.fdummy,
    wind_v_field=ws.fdummy,
    wind_w_field=ws.fdummy,
    mag_u_field=ws.fdummy,
    mag_v_field=ws.fdummy,
    mag_w_field=ws.fdummy,
    bad_partition_functions_ok=ws.bad_partition_functions_ok,
)
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupSetup(atmosphere_dim=ws.atmdim_1D, t_field=ws.t_1D, vmr_field=ws.vmr_1D)
ws.abs_lookupCalc()
# abs_xsec_agenda_checkedCalc
# atmfields_checkedCalc( bad_partition_functions_ok=bad_partition_functions_ok )
# abs_lookupSetup
# abs_lookupCalc
