################################################################################
#                                                                              #
# This is a (plug&play-type) include file. The USER is NOT supposed to MODIFY  #
# it, but choose another include file to suit his/her needs.                   #
#                                                                              #
################################################################################
#                                                                              #
# This INCLUDE file is for                                                     #
#  - using on-the-fly calculated absorption coefficients                       #
#  - including CIA continua (works also, if no CIA species present)            #
#                                                                              #
# It performs the following actions:                                           #
#  - sets propmat_clearsky_agenda: use on-the-fly absorption e                 #
#  - sets abs_xsec_agenda: include CIA and allow temperature extrapolation of  #
#     CIA data (works fine even in absence of CIA species)                     #
#  - read CIA input data from toolbox data package                             #
#  - read required spectroscopic line files from toolbox data package          #
#                                                                              #
# It requires the following input:                                             #
#   abs_species               as the WSV                                       #
#                                                                              #
# It provides following output (in parentheses: side products):                #
#   abs_lookup                as the WSV                                       #
#   propmat_clearsky_agenda   as the WSA                                       #
#   abs_xsec_agenda           as the WSA                                       #
#   (abs_lines)               as the WSV                                       #
#   (abs_lines_per_species)   as the WSV                                       #
#   (abs_cia_data)            as the WSV                                       #
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
# use on-the-fly calculated abs.coeffs
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# for CIA allow temperature extrapolation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__withCIAextraT)
# for the selected abs_species get the line and CIA data
#####
ws.ReadXML(ws.abs_cia_data, "spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml.gz")
ws.abs_linesReadFromSplitArtscat(
    basename="spectroscopy/Perrin/", fmin=0.0, fmax=4000000000000.0
)
# abs_linesReadFromHitran( filename="your-HITRAN-file", fmin=0, fmax=4e12 )
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
