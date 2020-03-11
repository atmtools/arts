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
#                                                                              #
# It performs the following actions:                                           #
#  - sets propmat_clearsky_agenda: use absorption lookup table                 #
#  - sets abs_xsec_agenda: include CIA and allow temperature extrapolation of  #
#     CIA data (works fine even in absence of CIA species)                     #
#  - read CIA input data from toolbox data package                             #
#  - read required spectroscopic line files from toolbox data package          #
#  - calculate the lookup table                                                #
#                                                                              #
# It requires the following input:                                             #
#   abs_species               as the WSV                                       #
#   atmosphere_dim            as the WSV                                       #
#   p/lat/lon_grid            as the WSV                                       #
#   t/z/vmr_field             as the WSV                                       #
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
ws.ReadSplitARTSCAT(basename="spectroscopy/Perrin/", fmin=0.0, fmax=4000000000000.0)
# abs_linesReadFromHitran( filename="your-HITRAN-file", fmin=0, fmax=4e12 )
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.abs_lookupSetup()
ws.abs_lookupCalc()
