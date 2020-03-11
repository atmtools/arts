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
# This case is for Mars and specifically tests
#
# (CASES 1-72)
#  - 72 Mars scenarions: 4seasons x 2daytimes x 3dustloads x 3solaractivities
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
ws.execute_controlfile("general/planet_mars.arts")
ws.NumericCreate("pmin")
ws.NumericCreate("pmax")
# VectorCreate( ally )
# VectorSet( ally, [] )
# WriteXML( in=ally )
# some basic RT settings
#####
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([2.208e11]))
ws.NumericSet(ws.pmin, 1e-05)
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
ws.StringSet(ws.basename, "planets/Mars/MPS/")
# Arrays with (sub)case names
ws.ArrayOfStringCreate("seasoncasearray")
ws.ArrayOfStringSet(
    ws.seasoncasearray, ["Mars.Ls0", "Mars.Ls90", "Mars.Ls180", "Mars.Ls270"]
)
ws.ArrayOfStringCreate("timecasearray")
ws.ArrayOfStringSet(ws.timecasearray, [".day", ".night"])
ws.ArrayOfStringCreate("dustcasearray")
ws.ArrayOfStringSet(ws.dustcasearray, [".dust-high", ".dust-low", ".dust-medium"])
ws.ArrayOfStringCreate("solarcasearray")
# regarding cloud data, all solar cases are identical. hence, we testprocess only one of them.
ws.ArrayOfStringSet(ws.solarcasearray, [".sol-avg"])
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
ws.StringSet(ws.cloudpath, "planets/Mars/SAT/")
ws.ArrayOfStringSet(
    ws.pndcasearray,
    [
        "pnd_field__Dust__small-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__small-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__small-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__large-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__large-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__large-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__verylarge-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__verylarge-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__verylarge-size-bulk__verywell-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__well-mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__mixed__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__confined__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__very-confined__tau1075cm-1_1e-01.xml",
        "pnd_field__Dust__medium-size-bulk__highly-confined__tau1075cm-1_1e-01.xml",
        "pnd_field__CO2ice__mesospheric-day-bulk__meso-gauss-profile__tau1024nm_2e-01.xml",
        "pnd_field__CO2ice__mesospheric-day-bulk__meso-gauss-profile__tau1024nm_2e-01.xml",
        "pnd_field__CO2ice__mesospheric-day-bulk__meso-gauss-profile__tau1024nm_2e-01.xml",
        "pnd_field__CO2ice__polarnight-chan1-bulk__20km-narrow-box-profile__ext1024nm_3e-01.xml",
        "pnd_field__CO2ice__polarnight-chan1-bulk__20km-narrow-box-profile__ext1024nm_3e-01.xml",
        "pnd_field__CO2ice__polarnight-chan1-bulk__20km-narrow-box-profile__ext1024nm_3e-01.xml",
        "pnd_field__CO2ice__polarnight-chan1-bulk__8km-narrow-box__ext1024nm_3e-01.xml",
        "pnd_field__CO2ice__polarnight-chan1-bulk__8km-tower-box-profile__ext1024nm_3e-01.xml",
        "pnd_field__CO2ice__polarnight-chan4-bulk__10km-narrow-box-profile__tau1024nm_1e+00.xml",
        "pnd_field__CO2ice__polarnight-chan4-bulk__10km-narrow-box-profile__tau1024nm_1e+00.xml",
        "pnd_field__CO2ice__polarnight-chan4-bulk__10km-narrow-box-profile__tau1024nm_1e+00.xml",
        "pnd_field__CO2ice__polarnight-chan4-bulk__4km-narrow-box__tau1024nm_1e+00.xml",
        "pnd_field__H2Oice__Type1-bulk__high-wide-profile__tau825cm-1_5e-02.xml",
        "pnd_field__H2Oice__Type1-bulk__low-narrow-profile__tau825cm-1_5e-02.xml",
        "pnd_field__H2Oice__Type2-bulk__high-wide-profile__tau825cm-1_2e-01.xml",
    ],
)
ws.ArrayOfStringSet(
    ws.ssdcasearray,
    [
        "Mars.scat_data__Dust__small-size-bulk__RI-minabs.xml",
        "Mars.scat_data__Dust__small-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__small-size-bulk__RI-maxabs.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-minabs.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-maxabs.xml",
        "Mars.scat_data__Dust__large-size-bulk__RI-minabs.xml",
        "Mars.scat_data__Dust__large-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__large-size-bulk__RI-maxabs.xml",
        "Mars.scat_data__Dust__verylarge-size-bulk__RI-minabs.xml",
        "Mars.scat_data__Dust__verylarge-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__verylarge-size-bulk__RI-maxabs.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__Dust__medium-size-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__mesospheric-day-bulk__RI-maxabs.xml",
        "Mars.scat_data__CO2ice__mesospheric-day-bulk__RI-minabs.xml",
        "Mars.scat_data__CO2ice__mesospheric-day-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__polarnight-chan1-bulk__RI-maxabs.xml",
        "Mars.scat_data__CO2ice__polarnight-chan1-bulk__RI-minabs.xml",
        "Mars.scat_data__CO2ice__polarnight-chan1-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__polarnight-chan1-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__polarnight-chan1-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__polarnight-chan4-bulk__RI-maxabs.xml",
        "Mars.scat_data__CO2ice__polarnight-chan4-bulk__RI-minabs.xml",
        "Mars.scat_data__CO2ice__polarnight-chan4-bulk__RI-std.xml",
        "Mars.scat_data__CO2ice__polarnight-chan4-bulk__RI-std.xml",
        "Mars.scat_data__H2Oice__Type1-bulk__RI-maetzler06.xml",
        "Mars.scat_data__H2Oice__Type1-bulk__RI-maetzler06.xml",
        "Mars.scat_data__H2Oice__Type2-bulk__RI-maetzler06.xml",
    ],
)
ws.abs_speciesSet(species=["CO2", "H2O", "CO"])
ws.abs_linesReadFromSplitArtscat(
    ws.abs_lines, ws.abs_species, "spectroscopy/Perrin/", 0.0, 1000000000000.0
)
ws.abs_lines_per_speciesCreateFromLines()
# create LUT valid for all cases
ws.abs_lookupSetupWide(
    p_min=ws.pmin, p_max=1000.0, t_min=105.0, t_max=405.0, h2o_min=0.0, h2o_max=0.0004
)
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupCalc()
ws.abs_lookupAdapt()
# sensor specification (that's case independent...): LOS zenith angle and altitude
ws.MatrixSet(ws.sensor_los, array([[121.8]]))
# tanh~1km
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
# now we go with several nested foorloop through the different cases.
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

ws.AgendaCreate("forloop_agenda_dust")


@arts_agenda
def forloop_agenda_dust(ws):
    # construct atmcase name III (Mars.LsXX.YY.dust-ZZ)
    ws.Extract(ws.casefull, ws.dustcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # construct pnd location name (same Mars.LsXX.YY.dust-ZZ as clear-sky atmo)
    ws.Append(ws.pndname, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.pndname, ws.caseext)
    ws.Append(ws.pndname, ws.atmcase)
    ws.StringSet(ws.caseext, ".")
    ws.Append(ws.pndname, ws.caseext)
    # Print( pndname, 0 )
    # keep the casestring till dust and make upper-level folder name
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    # construct atmcase name IV (Mars.LsXX.YY.dust-ZZ.sol-WW)
    ws.Extract(ws.casefull, ws.solarcasearray, 0)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    ws.Append(ws.basename, ws.atmcase)
    # Print( basename, 0 )
    ws.AtmRawRead(basename=ws.basename)
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


ws.forloop_agenda_dust = forloop_agenda_dust

ws.AgendaCreate("forloop_agenda_time")


@arts_agenda
def forloop_agenda_time(ws):
    # construct atmcase name II (Mars.LsXX.d/n)
    ws.Extract(ws.casefull, ws.timecasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_dust)
    ws.nelemGet(ws.ncases, ws.dustcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_time = forloop_agenda_time

ws.AgendaCreate("forloop_agenda_season")


@arts_agenda
def forloop_agenda_season(ws):
    # construct atmcase name I (Mars.LsXX)
    ws.Extract(ws.casefull, ws.seasoncasearray, ws.forloop_index)
    ws.Copy(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_time)
    ws.nelemGet(ws.ncases, ws.timecasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_season = forloop_agenda_season

ws.StringSet(ws.iy_unit, "PlanckBT")
ws.Copy(ws.pndname, ws.cloudpath)
ws.Copy(ws.ssdname, ws.cloudpath)
ws.nelemGet(ws.ncases, ws.seasoncasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_season)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
