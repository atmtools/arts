import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# First run basic tests of implementation
#####
ws.TMatrixTest()
# Define particle and ssd grids
#####
ws.VectorCreate("data_za_grid")
ws.VectorNLinSpace(ws.data_za_grid, 19, 0.0, 180.0)
ws.VectorCreate("data_aa_grid")
ws.VectorNLinSpace(ws.data_aa_grid, 19, 0.0, 180.0)
ws.VectorCreate("data_f_grid")
ws.VectorSet(ws.data_f_grid, array([2.3e11, 2.4e11]))
ws.VectorCreate("data_t_grid")
ws.VectorSet(ws.data_t_grid, array([220.0, 250.0, 270.0]))
# complex_refr_indexIceMatzler06(
#  data_f_grid = data_f_grid,
#  data_T_grid = data_t_grid )
ws.ReadXML(
    ws.complex_refr_index, "../refice/TestRefice.complex_refr_indexREFERENCE.xml"
)
ws.StringCreate("part_shape")
ws.StringSet(ws.part_shape, "cylindrical")
ws.NumericCreate("part_dveq")
ws.NumericSet(ws.part_dveq, 0.0001)
# [m]
ws.NumericCreate("part_ar")
ws.NumericSet(ws.part_ar, 2.0)
ws.NumericCreate("part_mass")
ws.NumericSet(ws.part_mass, 4.79983e-10)
# [kg]; m = Pi/6. * dveq^3 * density_ice
ws.SingleScatteringDataCreate("ref")
# TMatrix calculation: AZIMUTHALLY RANDOMLY ORIENTED (here: horizontally aligned) PARTICLE
#####
ws.scat_data_singleTmatrix(
    shape=ws.part_shape,
    diameter_volume_equ=ws.part_dveq,
    aspect_ratio=ws.part_ar,
    mass=ws.part_mass,
    ptype="azimuthally_random",
    data_f_grid=ws.data_f_grid,
    data_t_grid=ws.data_t_grid,
    data_za_grid=ws.data_za_grid,
    data_aa_grid=ws.data_aa_grid,
)
# Write data to file. And read from file (making sure, format is ok).
#
ws.WriteXML("ascii", ws.scat_data_single, "TestTMatrix.scat_data_single.azi-random.xml")
ws.WriteXML("ascii", ws.scat_meta_single)
ws.ReadXML(ws.scat_data_single, "TestTMatrix.scat_data_single.azi-random.xml")
ws.ReadXML(ws.scat_meta_single)
# Compare to reference data to ensure calcs provide expected results.
#
ws.ReadXML(ws.ref, "TestTMatrix.azi-random.ssdREFERENCE.xml")
ws.Compare(ws.scat_data_single, ws.ref, 2e-10)
# TMatrix calculation: TOTALLY RANDOMLY ORIENTED PARTICLE
#####
ws.scat_data_singleTmatrix(
    shape=ws.part_shape,
    diameter_volume_equ=ws.part_dveq,
    aspect_ratio=ws.part_ar,
    mass=ws.part_mass,
    ptype="totally_random",
    data_f_grid=ws.data_f_grid,
    data_t_grid=ws.data_t_grid,
    data_za_grid=ws.data_za_grid,
)
# Write data to file. And read from file (making sure, format is ok).
#
ws.WriteXML("ascii", ws.scat_data_single, "TestTMatrix.scat_data_single.tot-random.xml")
ws.WriteXML("ascii", ws.scat_meta_single)
ws.ReadXML(ws.scat_data_single, "TestTMatrix.scat_data_single.tot-random.xml")
ws.ReadXML(ws.scat_meta_single)
# Compare to reference data to ensure calcs provide expected results.
#
ws.ReadXML(ws.ref, "TestTMatrix.tot-random.ssdREFERENCE.xml")
ws.Compare(ws.scat_data_single, ws.ref, 1e-12)
# Check particle size parameter conversions
# (volume equivalent sphere diameter and maximum dimension)
#
# that's actually not limited to TMatrix calcs, but a non-spherical particles
# issue. anyway, we keep that test here.
#####
ws.NumericCreate("dveq0")
ws.NumericCreate("dveq")
ws.NumericCreate("dmax")
ws.NumericCreate("volume")
ws.NumericCreate("darea")
ws.StringCreate("shape")
ws.NumericCreate("aratio")
#
ws.NumericSet(ws.dveq0, 0.0001)
ws.StringSet(ws.shape, "cylindrical")
#
ws.NumericSet(ws.aratio, 3.45)
ws.diameter_maxFromDiameter_volume_equ(ws.dmax, ws.darea, ws.shape, ws.dveq0, ws.aratio)
ws.diameter_volume_equFromDiameter_max(ws.dveq, ws.volume, ws.shape, ws.dmax, ws.aratio)
ws.Compare(ws.dveq, ws.dveq0, 1e-12)
#
ws.NumericSet(ws.aratio, 0.22)
ws.diameter_maxFromDiameter_volume_equ(ws.dmax, ws.darea, ws.shape, ws.dveq0, ws.aratio)
ws.diameter_volume_equFromDiameter_max(ws.dveq, ws.volume, ws.shape, ws.dmax, ws.aratio)
ws.Compare(ws.dveq, ws.dveq0, 1e-12)
ws.StringSet(ws.shape, "spheroidal")
#
ws.NumericSet(ws.aratio, 3.45)
ws.diameter_maxFromDiameter_volume_equ(ws.dmax, ws.darea, ws.shape, ws.dveq0, ws.aratio)
ws.diameter_volume_equFromDiameter_max(ws.dveq, ws.volume, ws.shape, ws.dmax, ws.aratio)
ws.Compare(ws.dveq, ws.dveq0, 1e-12)
#
ws.NumericSet(ws.aratio, 0.22)
ws.diameter_maxFromDiameter_volume_equ(ws.dmax, ws.darea, ws.shape, ws.dveq0, ws.aratio)
ws.diameter_volume_equFromDiameter_max(ws.dveq, ws.volume, ws.shape, ws.dmax, ws.aratio)
ws.Compare(ws.dveq, ws.dveq0, 1e-12)
