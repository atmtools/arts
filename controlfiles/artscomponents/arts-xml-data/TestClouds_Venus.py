#
# Testing functionality (meeting format requirements, etc.) of cloud/dust fields
#  data.
#
# General test setup: Read a set of corresponding single scattering and particle
#  number density data, regrid to defined cloudbox, perform cloud-related
#  checks with cloudbox_checkedCalc, and apply cloud data in a (FOS, particle
#  absorption only) RT calculation.
#
#
# This case is for Venus and specifically tests
#
# (CASES A-E)
#  - the five Venus scenarions: Venus.spicav.night, Venus.spicav.night_cold,
#     Venus.vira.night, Venus.vira.day, Venus.vira.day_highlat
#  - per scenario use each single scattering data file (ssd) and particle number
#     density (pnd) field at least once (results in 32 combinations in total).
#     Various ssd can belong to one pnd field due to varied refractive index
#     assumptions. Various pnd fields can go with one ssd file due to varied
#     altitude distribution.
#  - regridding pnd fields to the computational grid within cloudbox
#  - perform cloud related internal check
#  - perfrom an RT calculation using FOS scattering method in zero-th scattering
#     order (particle absorption only) mode
#
# Jana Mendrok 2013-10-06

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_venus.arts")
ws.NumericCreate("pmin")
ws.NumericCreate("pmax")
ws.StringCreate("dummyfile")
ws.StringSet(ws.dummyfile, "dummy.tmp")
ws.ArrayOfGriddedField3Create("addvmr")
ws.Touch(ws.addvmr)
ws.WriteXML(ws.output_file_format, ws.addvmr, ws.dummyfile, 0)
# VectorCreate( ally )
# VectorSet( ally, [] )
# WriteXML( in=ally )
# some basic RT settings
#####
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([2.208e11]))
ws.NumericSet(ws.pmin, 0.1)
ws.NumericSet(ws.pmax, 1e99)
# and some further settings in order to be able to do an RT calc
#####
ws.jacobianOff()
ws.sensorOff()
# and agenda settings needed for RT calc
#####
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field,
)
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
ws.StringCreate("casename")
ws.IndexCreate("ncases")
ws.IndexCreate("itmp")
# set basic case folder
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Venus/MPS/")
# Array with case names
ws.ArrayOfStringCreate("atmcasearray")
ws.ArrayOfStringSet(
    ws.atmcasearray,
    [
        "Venus.spicav.night",
        "Venus.spicav.night_cold",
        "Venus.vira.night",
        "Venus.vira.day",
        "Venus.vira.day_highlat",
    ],
)
# case unspecific cloud specification (the pndfields are specific. however, we
# need them for ScatElementsPndAndScatAdd. so we use a dummy version, which we will
# overwrite later on.
ws.StringCreate("cloudpath")
ws.StringCreate("ssdname")
ws.StringCreate("pndname")
ws.StringCreate("psdname")
ws.StringCreate("psdprofname")
ws.ArrayOfStringCreate("pndcasearray")
ws.ArrayOfStringCreate("ssdcasearray")
ws.ArrayOfStringCreate("assd")
ws.StringSet(ws.cloudpath, "planets/Venus/SAT/")
ws.ArrayOfStringSet(
    ws.pndcasearray,
    [
        "pnd_field__H2SO4__LowerCloud-Mode1-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__LowerCloud-Mode1-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__LowerCloud-Mode2-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__LowerCloud-Mode2-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__LowerCloud-Mode3-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__LowerCloud-Mode3-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__LowerHaze-Mode1-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__LowerHaze-Mode1-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode1-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode1-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode2-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode2-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode3-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__MiddleCloud-Mode3-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__UpperCloud-Mode1-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__UpperCloud-Mode1-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__UpperCloud-Mode2-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__UpperCloud-Mode2-bulk__KH80-Nstd-box-profile.xml",
        "pnd_field__H2SO4__UpperHaze-Mode1-bulk__KH80-Nalt-box-profile.xml",
        "pnd_field__H2SO4__UpperHaze-Mode1-bulk__KH80-Nstd-box-profile.xml",
    ],
)
ws.ArrayOfStringSet(
    ws.ssdcasearray,
    [
        "Venus.scat_data__H2SO4__LowerCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerCloud-Mode3-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerCloud-Mode3-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerHaze-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__LowerHaze-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode3-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__MiddleCloud-Mode3-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperCloud-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperCloud-Mode2-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperHaze-Mode1-bulk__RI-std.xml",
        "Venus.scat_data__H2SO4__UpperHaze-Mode1-bulk__RI-std.xml",
    ],
)
ws.ArrayOfStringCreate("allspecies")
ws.ArrayOfStringCreate("basespecies")
ws.ArrayOfStringCreate("addspecies")
ws.ArrayOfStringCreate("addspeciesnames")
ws.ArrayOfStringSet(ws.basespecies, ["CO2", "CO"])
ws.ArrayOfStringSet(ws.addspecies, ["H2O"])
ws.ArrayOfStringSet(ws.addspeciesnames, [".H2O_mid"])
ws.Append(ws.allspecies, ws.basespecies)
ws.Append(ws.allspecies, ws.addspecies)
ws.abs_speciesSet(species=ws.allspecies)
ws.abs_linesReadFromSplitArtscat(
    ws.abs_lines, ws.abs_species, "spectroscopy/Perrin/", 0.0, 1000000000000.0
)
ws.abs_lines_per_speciesCreateFromLines()
# create LUT valid for all cases
ws.abs_lookupSetupWide(
    p_min=ws.pmin,
    p_max=10000000.0,
    t_min=140.0,
    t_max=740.0,
    h2o_min=0.0,
    h2o_max=4e-05,
)
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupCalc()
ws.abs_lookupAdapt()
# sensor specification (that's case independent...): LOS zenith angle and altitude
# VectorCreate( ztan )
# VectorSet( ztan, [30e3] )
# nelemGet( itmp, ztan )
# MatrixSetConstant( sensor_pos, itmp, 1, 600e3 )
# VectorZtanToZa1D( ztan, sensor_pos, refellipsoid, atmosphere_dim, ztan )
# Print( ztan )
# Matrix1ColFromVector( sensor_los, ztan )
ws.MatrixSet(ws.sensor_los, array([[113.9]]))
# tanh~30km
ws.nrowsGet(ws.itmp, ws.sensor_los)
ws.MatrixSetConstant(ws.sensor_pos, ws.itmp, 1, 600000.0)
# case unspecific surface settings
ws.VectorSet(ws.surface_scalar_reflectivity, array([0.4]))
# scattering solver required settings
ws.ReadXML(ws.fos_scatint_angles, "scattering/fosangles_360.xml")
ws.VectorSet(
    ws.fos_iyin_za_angles,
    array(
        [
            0.0,
            30.0,
            50.0,
            80.0,
            90.0,
            91.0,
            92.0,
            93.0,
            95.0,
            105.0,
            110.0,
            130.0,
            150.0,
            180.0,
        ]
    ),
)
#####
# CASES A-C (night) and D-E (day)
#####
# we go with a foorloop through the different cases
ws.AgendaCreate("forloop_agenda_particles")


@arts_agenda
def forloop_agenda_particles(ws):
    # cloud specification
    # we set a the cloudbox over the whole atmosphere. as we might include
    # mesospheric clouds.
    ws.cloudboxSetManually(
        p1=ws.pmax, p2=ws.pmin, lat1=0.0, lat2=0.0, lon1=0.0, lon2=0.0
    )
    ws.Extract(ws.psdprofname, ws.pndcasearray, ws.forloop_index)
    ws.Append(ws.pndname, ws.psdprofname)
    # Print( pndname, 0 )
    ws.Extract(ws.psdname, ws.ssdcasearray, ws.forloop_index)
    ws.Append(ws.ssdname, ws.psdname)
    ws.Append(ws.assd, ws.ssdname)
    # Print( assd, 0 )
    ws.ScatSpeciesInit()
    ws.ScatElementsPndAndScatAdd(scat_data_files=ws.assd, pnd_field_files=[""])
    ws.scat_dataCalc()
    # scat_dataCheck
    ws.ReadXML(ws.pnd_field_raw, ws.pndname)
    ws.pnd_fieldCalcFrompnd_field_raw(zeropadding=1)
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.scat_data_checkedCalc()
    ws.sensor_checkedCalc()
    ws.propmat_clearsky_agenda_checkedCalc()
    ws.yCalc()
    # Print( y, 0 )
    # ReadXML( out=ally )
    # Append( ally, y )
    # WriteXML( in=ally )


ws.forloop_agenda_particles = forloop_agenda_particles

ws.AgendaCreate("forloop_agenda_addspecies")


@arts_agenda
def forloop_agenda_addspecies(ws):
    ws.Extract(ws.caseext, ws.addspeciesnames, ws.forloop_index)
    ws.Append(ws.basename, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.basename)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.ReadXML(out=ws.addvmr, filename=ws.dummyfile)
    ws.Append(ws.addvmr, ws.gf3tmp)
    ws.WriteXML(ws.output_file_format, ws.addvmr, ws.dummyfile, 0)


ws.forloop_agenda_addspecies = forloop_agenda_addspecies

ws.AgendaCreate("forloop_agenda_atmcase")


@arts_agenda
def forloop_agenda_atmcase(ws):
    # construct atmcase name (Venus.CaseName)
    ws.Extract(ws.casefull, ws.atmcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # construct pnd location name (same Venus.CaseName as clear-sky atmo)
    ws.Append(ws.pndname, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.pndname, ws.caseext)
    ws.Append(ws.pndname, ws.atmcase)
    ws.StringSet(ws.caseext, ".")
    ws.Append(ws.pndname, ws.caseext)
    # Print( pndname, 0 )
    # keep the atmcasestring and make upper-level folder name
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    ws.Append(ws.basename, ws.atmcase)
    ws.abs_speciesSet(species=ws.basespecies)
    ws.AtmRawRead(basename=ws.basename)
    ws.abs_speciesAdd(species=ws.addspecies)
    ws.Delete(ws.addvmr)
    ws.Touch(ws.addvmr)
    ws.WriteXML(ws.output_file_format, ws.addvmr, ws.dummyfile, 0)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_addspecies)
    ws.nelemGet(ws.ncases, ws.addspecies)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
    ws.ReadXML(out=ws.addvmr, filename=ws.dummyfile)
    ws.Append(ws.vmr_field_raw, ws.addvmr)
    # now derive common p_grid and regrid atm fields to this
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.VectorCrop(ws.p_grid, ws.p_grid, ws.pmin, ws.pmax)
    ws.AtmFieldsCalc()
    # surface also needed for basics_checkedCalc
    ws.Extract(ws.z_surface, ws.z_field, 0)
    # and now the dusty calcs
    ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__FOSN0)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_particles)
    ws.nelemGet(ws.ncases, ws.pndcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_atmcase = forloop_agenda_atmcase

ws.StringSet(ws.iy_unit, "PlanckBT")
ws.Copy(ws.pndname, ws.cloudpath)
ws.Copy(ws.ssdname, ws.cloudpath)
ws.nelemGet(ws.ncases, ws.atmcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_atmcase)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
