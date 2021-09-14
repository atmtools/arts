#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Perform a radiative transfer simulation for microwave sensors.

 Index met_mm_accuracy selects the number of frequencies to be used
 from the available accuracies. This is set by the user in the main controlfile
 Must be documented here:
 N  Accuracy   Speed           Number of frequencies/channel
 0: fast       very-very fast  1
 1: normal     fast            varying
 2: high       slow            varying
 3: reference  extremely slow  varying

More information on met_mm_accuracy can be found in the ARTS-Documentation
web-page in the section Technical Reports.
http://arts.mi.uni-hamburg.de/docs/met_mm_setup_documentation.pdf

Parameters:

Returns:
"""


def testmetmm(sensor_name="sensor_amsub"):
    import pyarts as py
    import sensor_descriptions
    #import importlib
    #sensor_d = importlib.import_module("sensor_descriptions." + "sensor_amsub")

    verbosity = 2
    ws = py.workspace.Workspace(verbosity)
    ws.execute_controlfile("general/general.arts")
    ws.execute_controlfile("general/continua.arts")
    ws.execute_controlfile("general/agendas.arts")
    ws.execute_controlfile("general/planet_earth.arts")

    ws.verbositySetScreen(ws.verbosity, verbosity)

    ws.AtmosphereSet1D()
    ws.IndexSet(ws.stokes_dim, 1)
    ws.StringSet(ws.iy_unit, "PlanckBT")

    ws.ArrayOfIndexCreate("channels")
    ws.ArrayOfIndexCreate("viewing_angles")
    ws.ArrayOfIndexSet(ws.channels, [-1])
    ws.ArrayOfIndexSet(ws.viewing_angles, [0])

    ws.MatrixSetConstant(ws.sensor_pos, int(1), int(1), 850e3)
    ws.MatrixSetConstant(ws.sensor_los, int(1), int(1), float(180))

    ws.IndexCreate("met_mm_accuracy")
    ws.IndexSet(ws.met_mm_accuracy, int(1))

    # sensor_descriptions/prepare_metmm.arts
    ws.ArrayOfArrayOfIndexCreate("met_mm_available_accuracies")
    ws.ArrayOfIndexCreate("freq_number_tmp")

    ws.VectorCreate("freq_spacing_tmp")

    ws.ArrayOfIndexCreate("met_mm_freq_number")
    ws.NumericCreate("current_spacing")

    ws.VectorCreate("met_mm_freq_spacing")
    ws.IndexCreate("met_mm_nchannels")

    # sensor_descriptions/sensor_amsub.arts
    getattr(sensor_descriptions, sensor_name)(ws)

    ########    apply_metmm.arts
    ws.Select(ws.antenna_dlos, ws.antenna_dlos, ws.viewing_angles)
    ws.Select(ws.met_mm_backend, ws.met_mm_backend, ws.channels)
    ws.Select(ws.met_mm_polarisation, ws.met_mm_polarisation, ws.channels)
    ws.Select(ws.met_mm_freq_number, ws.met_mm_freq_number, ws.channels)
    ws.Select(ws.met_mm_freq_spacing, ws.met_mm_freq_spacing, ws.channels)

    ws.f_gridMetMM(
        freq_spacing=ws.met_mm_freq_spacing, freq_number=ws.met_mm_freq_number
    )

    ws.sensor_responseMetMM()

    #    common_metmm.arts
    ws.output_file_formatSetZippedAscii()
    ws.NumericSet(ws.ppath_lmax, float(250))

    # Agenda for scalar gas absorption calculation
    ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)

    # Surface
    ws.Copy(
        ws.surface_rtprop_agenda,
        ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
    )

    # (standard) emission calculation
    ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)

    # cosmic background radiation
    ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)

    # standard surface agenda (i.e., make use of surface_rtprop_agenda)
    ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)

    # sensor-only path
    ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)

    # no refraction
    ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)

    # Set propmat_clearsky_agenda to use lookup table
    ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
    #ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)

    # Spectroscopy
    species = [
        "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        "O2-66, O2-CIAfunCKDMT100",
        "N2,  N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
        "O3",
    ]
    # ws.abs_speciesSet( species=species )

    ws.abs_speciesSet(
        species=species,
    )

    ws.ReadARTSCAT(
        filename="instruments/metmm/abs_lines_metmm.xml.gz",
        fmin=0.0,
        fmax=float(1e99),
        globalquantumnumbers="",
        localquantumnumbers="",
        normalization_option="None",
        mirroring_option="None",
        population_option="LTE",
        lineshapetype_option="VP",
        cutoff_option="None",
        cutoff_value=750e9,
        linemixinglimit_value=-1.0,
    )

    ws.abs_linesSetCutoff(option="ByLine", value=750e9)
    ws.abs_linesSetNormalization(option="VVH")
    ws.abs_lines_per_speciesCreateFromLines()
    ws.abs_lines_per_speciesSetCutoffForSpecies(
        option="ByLine", value=5e9, species_tag="O3"
    )

    ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.4)

    ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")

    ws.batch_atm_fields_compactAddConstant(
        name="abs_species-O2", value=0.2095, condensibles=[]
    )
    ws.batch_atm_fields_compactAddConstant(
        name="abs_species-N2", value=0.7808, condensibles=[]
    )

    ws.abs_lookupSetupBatch()
    ws.abs_xsec_agenda_checkedCalc()
    ws.lbl_checkedCalc()
    ws.abs_lookupCalc()

    @py.workspace.arts_agenda
    def ybatch_calc_agenda(ws):
        ws.Extract(
            ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index
        )
        # Split up *atm_fields_compact* to
        # generate p_grid, t_field, z_field, vmr_field:
        ws.AtmFieldsAndParticleBulkPropFieldFromCompact()

        # Optionally set Jacobian parameters.
        ws.jacobianOff()

        # No scattering
        ws.cloudboxOff()

        # get some surface properties from corresponding atmospheric fields
        ws.Extract(ws.z_surface, ws.z_field, int(0))
        ws.Extract(ws.t_surface, ws.t_field, int(0))

        # Checks
        #ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
        ws.atmfields_checkedCalc()
        ws.atmgeom_checkedCalc()
        ws.cloudbox_checkedCalc()
        ws.sensor_checkedCalc()
        # Perform RT calculations
        ws.yCalc()

    ws.Copy(ws.ybatch_calc_agenda, ybatch_calc_agenda)

    ws.nelemGet(v=ws.batch_atm_fields_compact)
    ws.IndexSet(ws.ybatch_n, ws.nelem)
    # ====================================================================

    # Execute the batch calculations:
    # First check, then execute the batch RT calculations
    ws.propmat_clearsky_agenda_checkedCalc()
    ws.ybatchCalc(robust=1)
    # ====================================================================

    # Store results
    ws.WriteXML("ascii", ws.ybatch, "TestMetMM-pyarts.ybatch.xml")
    ws.WriteXML("ascii", ws.f_grid, "TestMetMM-pyarts.f_grid.xml")

    # Verify results
    ws.ArrayOfVectorCreate("ybatch_ref")
    ws.ReadXML( ws.ybatch_ref, "ybatchREFERENCE.xml")

    ws.Compare( ws.ybatch, ws.ybatch_ref, 0.01,
         "Total BT should be close to the reference values")

    print("Success! We reached the finish!")


if __name__ == "__main__":
    # normal way to run
    testmetmm(sensor_name="sensor_amsub")

    # test all sensors
    if 1==0:
        list = [
            "sensor_amsua",
            "sensor_amsub",
            "sensor_atms",
            "sensor_deimos",
            "sensor_hatpro",
            "sensor_ismar_downward",
            "sensor_upward",
            "sensor_marss",
            "sensor_mhs",
            "sensor_mwhs2",
            "sensor_saphir"
        ]
        for name in nlist:
            testmetmm(sensor_name=name)




