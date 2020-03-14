# DEFINITIONS:  -*-sh-*-
#
# ARTS control file for testing 1D propagation path calculations
# with different refractive index calculation methods (None, MicrowavesEarth, MicrowavesGeneral)
# for the different planets (particularly dependence on atmospheric composition)
#
# Jana Mendrok 2013-02-25

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# Use sensor_pos/los to store each case
ws.VectorCreate("za")
ws.VectorCreate("ztan")
ws.VectorSet(ws.ztan, np.array([50000.0, 46000.0, 34000.0, 10000.0, 5000.0, 2500.0]))
ws.nelemGet(ws.nelem, ws.ztan)
ws.MatrixSetConstant(ws.sensor_pos, ws.nelem, 1, 600000.0)
ws.VectorSet(ws.rte_pos2, [])
ws.IndexCreate("ilast")
ws.nrowsGet(ws.ilast, ws.sensor_pos)
ws.IndexStepDown(ws.ilast, ws.ilast)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.IndexSet(ws.stokes_dim, 1)
ws.AtmosphereSet1D()
ws.jacobianOff()
ws.cloudboxOff()
ws.NumericSet(ws.ppath_lmax, 5000.0)
# A dummy frequency grid
ws.VectorSet(ws.f_grid, np.array([1.0e10]))
ws.IndexCreate("tpoutlev")
ws.IndexSet(ws.tpoutlev, 1)
ws.IndexCreate("outlev")
ws.IndexSet(ws.outlev, 0)
ws.NumericCreate("tanh1")
ws.NumericCreate("tanlat1")
ws.VectorCreate("tp1")
ws.NumericCreate("tanh2")
ws.NumericCreate("tanlat2")
ws.VectorCreate("tp2")
ws.NumericCreate("diff")
ws.StringCreate("infostring")


@arts_agenda
def forloop_agenda(ws):
    ws.VectorExtractFromMatrix(ws.rte_pos, ws.sensor_pos, ws.forloop_index, "row")
    ws.VectorExtractFromMatrix(ws.rte_los, ws.sensor_los, ws.forloop_index, "row")
    ###
    # CASE 0-0: geometric path (no refraction)
    ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
    ws.ppathCalc()
    ws.TangentPointExtract(tan_pos=ws.tp1)
    ws.TangentPointPrint(level=ws.tpoutlev)
    ###
    # CASE 0-1: refracted path with n==1 (no refraction)
    ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
    ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__NoRefrac)
    ws.ppathCalc()
    ws.TangentPointExtract(tan_pos=ws.tp2)
    ws.TangentPointPrint(level=ws.tpoutlev)
    # Compare geom and n==1 results
    ws.Extract(ws.tanh1, ws.tp1, 0)
    ws.Extract(ws.tanh2, ws.tp2, 0)
    ws.Compare(
        ws.tanh1,
        ws.tanh2,
        1.0,
        "Geometric and n==1 tangent point altitude differ too much.",
    )
    ws.Extract(ws.tanlat1, ws.tp1, 1)
    ws.Extract(ws.tanlat2, ws.tp2, 1)
    ws.Compare(
        ws.tanlat1,
        ws.tanlat2,
        0.04,
        "Geometric and n==1 tangent point latitude differ too much.",
    )
    ###
    # CASE 1: refracted path with n from MicrowavesEarth
    ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
    ws.ppathCalc()
    # WriteXMLIndexed( "ascii", forloop_index, ppath, "thayer" )
    ws.TangentPointExtract(tan_pos=ws.tp1)
    ws.TangentPointPrint(level=ws.tpoutlev)
    ###
    # CASE 2: refracted path with n from MicrowavesGeneral
    ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesGeneral)
    ws.ppathCalc()
    # WriteXMLIndexed( "ascii", forloop_index, ppath, "MicrowavesGeneral" )
    ws.TangentPointExtract(tan_pos=ws.tp2)
    ws.TangentPointPrint(level=ws.tpoutlev)
    # Print diff of MicrowavesEarth and MicrowavesGeneral
    ws.Extract(ws.tanh1, ws.tp1, 0)
    ws.Extract(ws.tanh2, ws.tp2, 0)
    ws.NumericScale(ws.tanh2, ws.tanh2, -1.0)
    ws.NumericAdd(ws.diff, ws.tanh1, ws.tanh2)
    ws.StringSet(
        ws.infostring,
        "Tangent altitude diff of MicrowavesEarth and MicrowavesGeneral [m]: ",
    )
    ws.Print(ws.infostring, ws.outlev)
    ws.Print(ws.diff, ws.outlev)
    # Extract( tanlat1, tp1, 1 )
    # Extract( tanlat2, tp2, 1 )
    # NumericScale( tanlat2, tanlat2, -1 )
    # NumericAdd( diff, tanlat1, tanlat2 )
    # StringSet( infostring, "Tangent latitude diff of MicrowavesEarth and MicrowavesGeneral [deg]: " )
    # Print( infostring, outlev )
    # Print( diff, outlev )


ws.forloop_agenda = forloop_agenda

#####
# CASE A: Earth
#####
# Main atmospheric contributors
#  When comparing effects of planetary composition, we don't want further
#  confusion by H2O. so we leave it out here.
# abs_speciesSet( species=["H2O","N2","O2"] )
ws.abs_speciesSet(species=["N2", "O2", "CO2"])
ws.refellipsoidEarth(ws.refellipsoid, "Sphere")
# convert (geometric) tangent altitudes to zenith angles at observer
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.ztan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
# Pressure grid rougly matching 0 to 80 km.
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
ws.AtmRawRead(basename="planets/Earth/Fascod/tropical/tropical")
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
#####
# CASE B: Mars
#####
# Main atmospheric contributors
ws.abs_speciesSet(species=["N2", "O2", "CO2", "H2"])
ws.refellipsoidMars(ws.refellipsoid, "Sphere")
# convert (geometric) tangent altitudes to zenith angles at observer
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.ztan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
# Pressure grid rougly matching 0 to 80 km.
ws.VectorNLogSpace(ws.p_grid, 41, 765.0, 0.08)
ws.AtmRawRead(
    basename="planets/Mars/MPS/Mars.Ls0.day.dust-medium/Mars.Ls0.day.dust-medium.sol-avg/Mars.Ls0.day.dust-medium.sol-avg"
)
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
#####
# CASE C: Venus
#####
# Main atmospheric contributors
ws.abs_speciesSet(species=["N2", "O2", "CO2"])
ws.refellipsoidVenus(ws.refellipsoid, "Sphere")
# convert (geometric) tangent altitudes to zenith angles at observer
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.ztan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
# Pressure grid rougly matching 0 to 80 km.
ws.VectorNLogSpace(ws.p_grid, 41, 9200000.0, 490.0)
ws.AtmRawRead(basename="planets/Venus/MPS/Venus.vira.day/Venus.vira.day")
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
#####
# CASE D: Jupiter
#####
# Main atmospheric contributors
ws.abs_speciesSet(species=["CO2", "H2", "He"])
ws.refellipsoidJupiter(ws.refellipsoid, "Sphere")
# convert (geometric) tangent altitudes to zenith angles at observer
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.ztan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
# Pressure grid rougly matching 0 to 80 km.
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1800.0)
ws.AtmRawRead(basename="planets/Jupiter/MPS/Jupiter.mean/Jupiter.mean")
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
