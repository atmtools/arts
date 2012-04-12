
add_custom_command(
  OUTPUT clean-cfiles-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-cfiles DEPENDS clean-cfiles-cmd)

#######################
### ARTS Components ###
arts_test_run_ctlfile(fast.ppath1d artscomponents/ppath/TestPpath1D.arts)
arts_test_run_ctlfile(fast.ppath2d artscomponents/ppath/TestPpath2D.arts)
arts_test_run_ctlfile(fast.ppath3d artscomponents/ppath/TestPpath3D.arts)

arts_test_run_ctlfile(fast.pencilbeam
                      artscomponents/pencilbeam/TestPencilBeam.arts)

arts_test_run_ctlfile(fast.clearsky artscomponents/clearsky/TestClearSky.arts)

arts_test_run_ctlfile(fast.groundbased
                      artscomponents/groundbased/TestGbased.arts)

arts_test_run_ctlfile(fast.transmission
                      artscomponents/transmission/TestTransmission.arts)

arts_test_run_ctlfile(fast.absorption
                      artscomponents/absorption/TestAbs.arts)

arts_test_run_ctlfile(fast.doit artscomponents/doit/TestDOIT.arts)

arts_test_run_ctlfile(slow.doitbatch
                      artscomponents/doitbatch/TestDOITBatch.arts)

arts_test_run_ctlfile(fast.montecarlo.dataprepare
                      artscomponents/montecarlo/TestMonteCarloDataPrepare.arts)

arts_test_run_ctlfile(slow.montecarlo.general
                      artscomponents/montecarlo/TestMonteCarloGeneral.arts)
arts_test_ctlfile_depends(slow.montecarlo.general fast.montecarlo.dataprepare)

arts_test_run_ctlfile(slow.montecarlo.generalgaussian
                      artscomponents/montecarlo/TestMonteCarloGeneralGaussian.arts)
arts_test_ctlfile_depends(slow.montecarlo.generalgaussian
                          fast.montecarlo.dataprepare)

arts_test_run_ctlfile(fast.montecarlo.rtecalc
                      artscomponents/montecarlo/TestRteCalcMC.arts)
arts_test_ctlfile_depends(fast.montecarlo.rtecalc fast.montecarlo.dataprepare)


