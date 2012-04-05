
add_custom_command(
  OUTPUT clean-tests-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-tests DEPENDS clean-tests-cmd)

arts_test_run_ctlfile(fast.abs            Abs/TestAbs.arts)
arts_test_run_ctlfile(fast.ppath1d        Ppath/TestPpath1D.arts)
arts_test_run_ctlfile(fast.ppath2d        Ppath/TestPpath2D.arts)
arts_test_run_ctlfile(fast.ppath3d        Ppath/TestPpath3D.arts)
arts_test_run_ctlfile(fast.pencilbeam     PencilBeam/TestPencilBeam.arts)
arts_test_run_ctlfile(fast.clearsky       ClearSky/TestClearSky.arts)
arts_test_run_ctlfile(fast.groundbased    GroundBased/TestGbased.arts)
arts_test_run_ctlfile(fast.transmission   Transmission/TestTransmission.arts)
arts_test_run_ctlfile(fast.radiolink      RadioLink/TestRadioLink.arts)
arts_test_run_ctlfile(fast.amsub          AMSU/TestAMSUB_fast.arts)
arts_test_run_ctlfile(fast.doit           DOIT/TestDOIT.arts)
arts_test_run_ctlfile(fast.odinsmr        OdinSMR/TestOdinSMR.arts)

arts_test_run_ctlfile(fast.montecarlo.dataprepare
                      MonteCarlo/TestMonteCarloDataPrepare.arts)

arts_test_run_ctlfile(slow.montecarlo.general
                      MonteCarlo/TestMonteCarloGeneral.arts)
arts_test_ctlfile_depends(slow.montecarlo.general
                          fast.montecarlo.dataprepare)

arts_test_run_ctlfile(slow.montecarlo.generalgaussian
                      MonteCarlo/TestMonteCarloGeneralGaussian.arts)
arts_test_ctlfile_depends(slow.montecarlo.generalgaussian
                          fast.montecarlo.dataprepare)

arts_test_run_ctlfile(fast.montecarlo.rtecalc MonteCarlo/TestRteCalcMC.arts)
arts_test_ctlfile_depends(fast.montecarlo.rtecalc fast.montecarlo.dataprepare)

arts_test_run_ctlfile(slow.amsub AMSU/TestAMSUB.arts)

arts_test_run_ctlfile(slow.doitbatch DOITBatch/TestDOITBatch.arts)

arts_test_run_ctlfile(fast.mhs MHS/TestMHS_fast.arts)
arts_test_run_ctlfile(slow.mhs MHS/TestMHS.arts)

arts_test_run_ctlfile(slow.wfuns Wfuns/TestWfuns.arts)

