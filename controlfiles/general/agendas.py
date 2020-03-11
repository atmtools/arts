# DEFINITIONS:  -*-sh-*-
#
############
# Agenda predefinitions
#
# There shall be NO pre-set (default) agendas in general.arts anymore, hiding
# them away from the user. Instead, the user is forced to explicitly set all
# the agendas needed. Here, we predefine a range of typical agenda settings,
# which the user can than copy into his/her setup using:
#
# Copy( agenda-name, predefined-agenda-name )
#
############
#
# Authors: Jana Mendrok
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
##################################
#    iy_main_agenda
##################################
###
ws.AgendaCreate("iy_main_agenda__Emission")


@arts_agenda
def iy_main_agenda__Emission(ws):
    ws.ppathCalc()
    ws.iyEmissionStandard()


ws.iy_main_agenda__Emission = iy_main_agenda__Emission

###
ws.AgendaCreate("iy_main_agenda__Transmission")


@arts_agenda
def iy_main_agenda__Transmission(ws):
    ws.Ignore(ws.iy_unit)
    ws.Ignore(ws.iy_id)
    ws.ppathCalc(cloudbox_on=0)
    ws.iyTransmissionStandard()


ws.iy_main_agenda__Transmission = iy_main_agenda__Transmission

###
ws.AgendaCreate("iy_main_agenda__Freqloop")


@arts_agenda
def iy_main_agenda__Freqloop(ws):
    ws.Ignore(ws.diy_dx)
    ws.Ignore(ws.iy_id)
    ws.Ignore(ws.iy_unit)
    ws.Ignore(ws.nlte_field)
    ws.Ignore(ws.cloudbox_on)
    ws.Ignore(ws.jacobian_do)
    ws.iyLoopFrequencies()
    ws.Touch(ws.ppath)


ws.iy_main_agenda__Freqloop = iy_main_agenda__Freqloop

###
ws.AgendaCreate("iy_main_agenda__ScattMC")


@arts_agenda
def iy_main_agenda__ScattMC(ws):
    ws.Ignore(ws.rte_pos2)
    ws.Ignore(ws.diy_dx)
    ws.Ignore(ws.iy_id)
    ws.Ignore(ws.nlte_field)
    ws.iyMC()
    ws.Touch(ws.ppath)


ws.iy_main_agenda__ScattMC = iy_main_agenda__ScattMC

##################################
#    iy_loop_freqs_agenda
##################################
###
ws.AgendaCreate("iy_loop_freqs_agenda__Emission")


@arts_agenda
def iy_loop_freqs_agenda__Emission(ws):
    ws.ppathCalc()
    ws.iyEmissionStandard()


ws.iy_loop_freqs_agenda__Emission = iy_loop_freqs_agenda__Emission

###
ws.AgendaCreate("iy_loop_freqs_agenda__Transmission")


@arts_agenda
def iy_loop_freqs_agenda__Transmission(ws):
    ws.Ignore(ws.iy_id)
    ws.ppathCalc()
    ws.iyTransmissionStandard()


ws.iy_loop_freqs_agenda__Transmission = iy_loop_freqs_agenda__Transmission

##################################
#    iy_space_agenda
##################################
###
# cosmic background (standard)
#
ws.AgendaCreate("iy_space_agenda__CosmicBackground")


@arts_agenda
def iy_space_agenda__CosmicBackground(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.MatrixCBR(ws.iy, ws.stokes_dim, ws.f_grid)


ws.iy_space_agenda__CosmicBackground = iy_space_agenda__CosmicBackground

##################################
#    iy_transmitter_agenda
##################################
###
# Unit (unpolarised) intensity
#
ws.AgendaCreate("iy_transmitter_agenda__UnitUnpolIntensity")


@arts_agenda
def iy_transmitter_agenda__UnitUnpolIntensity(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.MatrixUnitIntensity(ws.iy, ws.stokes_dim, ws.f_grid)


ws.iy_transmitter_agenda__UnitUnpolIntensity = iy_transmitter_agenda__UnitUnpolIntensity

###
# Unit polarised intensity
#
ws.AgendaCreate("iy_transmitter_agenda__UnitPolIntensity")


@arts_agenda
def iy_transmitter_agenda__UnitPolIntensity(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.iy_transmitterSinglePol()


ws.iy_transmitter_agenda__UnitPolIntensity = iy_transmitter_agenda__UnitPolIntensity

##################################
#    iy_surface_agenda
##################################
###
# use surface_rtprop_agenda through WSM iySurfaceRtpropAgenda (standard)
#
ws.AgendaCreate("iy_surface_agenda__UseSurfaceRtprop")


@arts_agenda
def iy_surface_agenda__UseSurfaceRtprop(ws):
    ws.SurfaceDummy()
    ws.iySurfaceRtpropAgenda()


ws.iy_surface_agenda__UseSurfaceRtprop = iy_surface_agenda__UseSurfaceRtprop

##################################
#    iy_cloudbox_agenda
##################################
# Note: iy_cloudbox_agenda only needs to be defined when a iy_main_agenda uses
#  iyEmissionStandard while a scattering calculation shall be performed. In
#  other words, when a scattering solver is applied that provides cloudbox_field
#  as output (these are, so far, the solvers DOIT, Disort, RT4).
#  Scattering solvers that have their dedicated iy* methods to be applied in the
#  iy_main_agenda (so far, iyMC and iyFOS) do not require iy_cloudbox_agenda to
#  be set.
###
# linear interpolation of the intensity field of the cloud box,
#
ws.AgendaCreate("iy_cloudbox_agenda__LinInterpField")


@arts_agenda
def iy_cloudbox_agenda__LinInterpField(ws):
    ws.iyInterpCloudboxField()


ws.iy_cloudbox_agenda__LinInterpField = iy_cloudbox_agenda__LinInterpField

###
# quadratic interpolation of the intensity field of the cloud box
#
ws.AgendaCreate("iy_cloudbox_agenda__QuadInterpField1D")


@arts_agenda
def iy_cloudbox_agenda__QuadInterpField1D(ws):
    ws.iyInterpCloudboxField(za_interp_order=4)


ws.iy_cloudbox_agenda__QuadInterpField1D = iy_cloudbox_agenda__QuadInterpField1D

##################################
#    ppath_agenda
##################################
###
# sensor-only path
#
ws.AgendaCreate("ppath_agenda__FollowSensorLosPath")


@arts_agenda
def ppath_agenda__FollowSensorLosPath(ws):
    ws.Ignore(ws.rte_pos2)
    ws.ppathStepByStep()


ws.ppath_agenda__FollowSensorLosPath = ppath_agenda__FollowSensorLosPath

###
# plane parellel version
#
ws.AgendaCreate("ppath_agenda__PlaneParallel")


@arts_agenda
def ppath_agenda__PlaneParallel(ws):
    ws.Ignore(ws.ppath_lraytrace)
    ws.Ignore(ws.rte_pos2)
    ws.Ignore(ws.t_field)
    ws.Ignore(ws.vmr_field)
    ws.Ignore(ws.f_grid)
    ws.ppathPlaneParallel()


ws.ppath_agenda__PlaneParallel = ppath_agenda__PlaneParallel

###
# transmitter-receiver (link) path
#
ws.AgendaCreate("ppath_agenda__TransmitterReceiverPath")


@arts_agenda
def ppath_agenda__TransmitterReceiverPath(ws):
    ws.Ignore(ws.cloudbox_on)
    ws.Ignore(ws.ppath_inside_cloudbox_do)
    ws.rte_losGeometricFromRtePosToRtePos2()
    ws.ppathFromRtePos2()


ws.ppath_agenda__TransmitterReceiverPath = ppath_agenda__TransmitterReceiverPath

##################################
#    ppath_step_agenda
##################################
###
# Geometrical path calculation (i.e., refraction neglected)
#
ws.AgendaCreate("ppath_step_agenda__GeometricPath")


@arts_agenda
def ppath_step_agenda__GeometricPath(ws):
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.ppath_lraytrace)
    ws.ppath_stepGeometric()


ws.ppath_step_agenda__GeometricPath = ppath_step_agenda__GeometricPath

###
# Refracted path calculation
#
ws.AgendaCreate("ppath_step_agenda__RefractedPath")


@arts_agenda
def ppath_step_agenda__RefractedPath(ws):
    ws.ppath_stepRefractionBasic()


ws.ppath_step_agenda__RefractedPath = ppath_step_agenda__RefractedPath

##################################
#    propmat_clearsky_agenda
##################################
ws.FlagOff(ws.propmat_clearsky_agenda_checked)
###
# "on-the-fly" absorption
#
ws.AgendaCreate("propmat_clearsky_agenda__OnTheFly")


@arts_agenda
def propmat_clearsky_agenda__OnTheFly(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()


ws.propmat_clearsky_agenda__OnTheFly = propmat_clearsky_agenda__OnTheFly

###
# "on-the-fly" absorption with Zeeman; requires precalculations
#
ws.AgendaCreate("propmat_clearsky_agenda__OnTheFly_ZeemanPreCalc")


@arts_agenda
def propmat_clearsky_agenda__OnTheFly_ZeemanPreCalc(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()
    ws.propmat_clearskyAddZeeman()


ws.propmat_clearsky_agenda__OnTheFly_ZeemanPreCalc = (
    propmat_clearsky_agenda__OnTheFly_ZeemanPreCalc
)

###
# "on-the-fly" absorption with Faraday rotation
#
ws.AgendaCreate("propmat_clearsky_agenda__OnTheFly_Faraday")


@arts_agenda
def propmat_clearsky_agenda__OnTheFly_Faraday(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()
    ws.propmat_clearskyAddFaraday()


ws.propmat_clearsky_agenda__OnTheFly_Faraday = propmat_clearsky_agenda__OnTheFly_Faraday

###
# absorption lookup table
# this could save (considerable) time, especially for batch calculations.
#
ws.AgendaCreate("propmat_clearsky_agenda__LookUpTable")


@arts_agenda
def propmat_clearsky_agenda__LookUpTable(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.rtp_nlte)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddFromLookup()


ws.propmat_clearsky_agenda__LookUpTable = propmat_clearsky_agenda__LookUpTable

###
# absorption lookup table with Zeeman
# this could save (considerable) time, especially for batch calculations.
#
ws.AgendaCreate("propmat_clearsky_agenda__LookUpTable_Zeeman")


@arts_agenda
def propmat_clearsky_agenda__LookUpTable_Zeeman(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddFromLookup()
    ws.propmat_clearskyAddZeeman()


ws.propmat_clearsky_agenda__LookUpTable_Zeeman = (
    propmat_clearsky_agenda__LookUpTable_Zeeman
)

##################################
#    abs_xsec_agenda
##################################
ws.FlagOff(ws.abs_xsec_agenda_checked)
###
# Without CIA (only explicit lines and
# "classical" continua / full abosorption models)
#
ws.AgendaCreate("abs_xsec_agenda__noCIA")


@arts_agenda
def abs_xsec_agenda__noCIA(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()
    ws.abs_xsec_per_speciesAddConts()


ws.abs_xsec_agenda__noCIA = abs_xsec_agenda__noCIA

###
# With CIA
#
ws.AgendaCreate("abs_xsec_agenda__withCIA")


@arts_agenda
def abs_xsec_agenda__withCIA(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()
    ws.abs_xsec_per_speciesAddConts()
    ws.abs_xsec_per_speciesAddCIA()


ws.abs_xsec_agenda__withCIA = abs_xsec_agenda__withCIA

###
# With CIA and sufficient temperature extrapolation
#
ws.AgendaCreate("abs_xsec_agenda__withCIAextraT")


@arts_agenda
def abs_xsec_agenda__withCIAextraT(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()
    ws.abs_xsec_per_speciesAddConts()
    ws.abs_xsec_per_speciesAddCIA(T_extrapolfac=6.0)


ws.abs_xsec_agenda__withCIAextraT = abs_xsec_agenda__withCIAextraT

##################################
#    refr_index_air_agenda
##################################
###
# no refraction (n==1.0)
#
ws.AgendaCreate("refr_index_air_agenda__NoRefrac")


@arts_agenda
def refr_index_air_agenda__NoRefrac(ws):
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.rtp_pressure)
    ws.Ignore(ws.rtp_temperature)
    ws.Ignore(ws.rtp_vmr)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)


ws.refr_index_air_agenda__NoRefrac = refr_index_air_agenda__NoRefrac

###
# refraction parameterization for microwaves + Earth (p,T,H2O)
#
ws.AgendaCreate("refr_index_air_agenda__GasMicrowavesEarth")


@arts_agenda
def refr_index_air_agenda__GasMicrowavesEarth(ws):
    ws.Ignore(ws.f_grid)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airMicrowavesEarth()


ws.refr_index_air_agenda__GasMicrowavesEarth = refr_index_air_agenda__GasMicrowavesEarth

###
# refraction parameterization
# valid in infrared for Earth (dry air only)
#
ws.AgendaCreate("refr_index_air_agenda__GasInfraredEarth")


@arts_agenda
def refr_index_air_agenda__GasInfraredEarth(ws):
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.rtp_vmr)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airInfraredEarth()


ws.refr_index_air_agenda__GasInfraredEarth = refr_index_air_agenda__GasInfraredEarth

###
# refraction parameterization from Newell&Baird (p,T; N2,O2,CO2,H2,He)
# valid for arbitrary planetary atmospheres in microwave
#
ws.AgendaCreate("refr_index_air_agenda__GasMicrowavesGeneral")


@arts_agenda
def refr_index_air_agenda__GasMicrowavesGeneral(ws):
    ws.Ignore(ws.f_grid)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airMicrowavesGeneral()


ws.refr_index_air_agenda__GasMicrowavesGeneral = (
    refr_index_air_agenda__GasMicrowavesGeneral
)

###
# refraction due to free electrons
#
ws.AgendaCreate("refr_index_air_agenda__FreeElectrons")


@arts_agenda
def refr_index_air_agenda__FreeElectrons(ws):
    ws.Ignore(ws.rtp_pressure)
    ws.Ignore(ws.rtp_temperature)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airFreeElectrons()


ws.refr_index_air_agenda__FreeElectrons = refr_index_air_agenda__FreeElectrons

###
# combined refraction from gases (Newell&Baird) and free electrons
# valid for arbitrary planetary atmospheres in microwave
#
ws.AgendaCreate("refr_index_air_agenda__GasMicrowavesGeneralAndElectrons")


@arts_agenda
def refr_index_air_agenda__GasMicrowavesGeneralAndElectrons(ws):
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airMicrowavesGeneral()
    ws.refr_index_airFreeElectrons()


ws.refr_index_air_agenda__GasMicrowavesGeneralAndElectrons = (
    refr_index_air_agenda__GasMicrowavesGeneralAndElectrons
)

###
# Combination of Microwaves + free electrons. OK for Earth and microwave.
#
ws.AgendaCreate("refr_index_air_agenda__GasMicrowavesEarthAndElectrons")


@arts_agenda
def refr_index_air_agenda__GasMicrowavesEarthAndElectrons(ws):
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airMicrowavesEarth()
    ws.refr_index_airFreeElectrons()


ws.refr_index_air_agenda__GasMicrowavesEarthAndElectrons = (
    refr_index_air_agenda__GasMicrowavesEarthAndElectrons
)

###
# all surface agendas are placed in s special file
#
ws.execute_controlfile("agendas_surface.arts")
