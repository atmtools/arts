################################################################################
#                                                                              #
# This is a (plug&play-type) include file. The USER is NOT supposed to MODIFY  #
# it, but choose another include file to suit his/her needs.                   #
#                                                                              #
################################################################################
#                                                                              #
# This INCLUDE file is for                                                     #
#  - considering refraction of "air"                                           #
#  - 1D calculations only                                                      #
#  - several viewing angles (incl. tangent altitudes) from constant observer   #
#     position                                                                 #
#  - for use with iyCalc (not yCalc)                                           #
#  - for receiver-only setups (no receiver-transmitter paths!)                 #
#                                                                              #
# It performs the following actions:                                           #
#  - sets ppath_agenda: receiver-viewingangle determined path (no transmitter) #
#  - sets ppath_step_agenda: refracted ppath calculation                       #
#  - sets refr_index_air_agenda: refr_index_airMicrowavesGeneral for air               #
#  - calculates viewing angles from given vector of tangent altitudes          #
#     (where the given tangent altitudes refer to the ones the receiver is     #
#     pointing at. actual refracted rays will have their (effective) tangent   #
#     point at lower altitudes than set by the user)                           #                                               #
#  - creates a common vector from (given) viewing angle vector and the viewing #
#     angles associated with the tangent altitudes                             #
#  - sets sensor positions: constant receiver position, empty transmitter      #
#                                                                              #
# It requires the following input:                                             #
#   viewzang        Vector; the viewing angles (1D, i.e. zenith only)          #
#   tanh            Vector; the tangent altitudes                              #
#   obsh            Numeric; the receiver altitude                             #
#   refellipsoid    as the WSV                                                 #
#   atmosphere_dim  as the WSV                                                 #
# It also uses (OVERWRITES!) sensor_pos                                        #
#                                                                              #
# It provides following output:                                                #
#   allzang         Vector; merged (1D) viewing angles resulting from viewang  #
#                    and calculated viewing angles resulting from tanh         #
#                                                                              #
# The file shall NOT be modified by the USER.                                  #
#                                                                              #
# This template creates (and makes internal use of) the following non-WSV:     #
#  (These are created here, i.e., they can not be used by earlier parts of the #
#  script or created again (it also implies, this include can only be included #
#  once in an ARTS run!). They can be used in later parts of the script,       #
#  though.)                                                                    #
#   zang            Vector                                                     #
#   ntanh           Index                                                      #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# receiver-viewingangle-path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# refracted path
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
# refraction from "air"
# (using general planet applicable method refr_index_airMicrowavesGeneral)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesGeneral)
# DO NOT MODIFY
# preprocessing of viewing geometry parameters
ws.VectorCreate("allzang")
ws.Copy(ws.allzang, ws.viewang)
ws.VectorCreate("zang")
ws.IndexCreate("ntanh")
ws.nelemGet(ws.ntanh, ws.tanh)
ws.MatrixSetConstant(ws.sensor_pos, ws.ntanh, 1, ws.obsh)
ws.VectorZtanToZa1D(ws.zang, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.tanh)
ws.Append(ws.allzang, ws.zang)
# WriteXML( in=allzang )
# Print( allzang )
# for use with yCalc
# Matrix1ColFromVector( sensor_los, allzang )
# nrowsGet( itmp, sensor_los )
# MatrixSetConstant( sensor_pos, itmp, 1, obsh )
# for use with looped iyCalc, i.e., we have to set rte_pos, not sensor_pos
ws.VectorSetConstant(ws.rte_pos, 1, ws.obsh)
ws.VectorSet(ws.rte_pos2, array([], dtype=float64))
