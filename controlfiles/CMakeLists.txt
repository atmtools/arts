
add_custom_command(
  OUTPUT clean-cfiles-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-cfiles DEPENDS clean-cfiles-cmd)

set (ARTS_CONTROLFILE_DIR share/arts/controlfiles)

install (DIRECTORY artscomponents
                   testdata
                   general
                   instruments
         DESTINATION ${ARTS_CONTROLFILE_DIR}
         PATTERN ".svn" EXCLUDE)


#######################
### ARTS Components ###
arts_test_run_ctlfile(fast.absorption artscomponents/absorption/TestAbs.arts)
arts_test_run_ctlfile(fast.absorption.doppler
                      artscomponents/absorption/TestAbsDoppler.arts)
arts_test_run_ctlfile(fast.absorption.compareArtscat3and4 
                      artscomponents/absorption/TestCompareArtscat3And4.arts)
arts_test_run_ctlfile(slow.isoratios artscomponents/absorption/TestIsoRatios.arts)

arts_test_run_ctlfile(fast.surface artscomponents/surface/TestSurface.arts)

arts_test_run_ctlfile(fast.ppath1d artscomponents/ppath/TestPpath1D.arts)
arts_test_run_ctlfile(fast.ppath2d artscomponents/ppath/TestPpath2D.arts)
arts_test_run_ctlfile(fast.ppath3d artscomponents/ppath/TestPpath3D.arts)

arts_test_run_ctlfile(fast.pencilbeam
                      artscomponents/pencilbeam/TestPencilBeam.arts)

arts_test_run_ctlfile(fast.clearsky artscomponents/clearsky/TestClearSky.arts)
arts_test_run_ctlfile(slow.clearsky2 artscomponents/clearsky/TestClearSky2.arts)

arts_test_run_ctlfile(fast.fgrid
                      artscomponents/groundbased/TestFgrid.arts)
arts_test_run_ctlfile(fast.groundbased
                      artscomponents/groundbased/TestGbased.arts)
arts_test_run_ctlfile(fast.groundbased.winds
                      artscomponents/groundbased/TestWinds.arts)

arts_test_run_ctlfile(fast.transmission
                      artscomponents/transmission/TestTransmission.arts)

arts_test_run_ctlfile(fast.faraday
                      artscomponents/faraday/TestFaradayRotation.arts)

arts_test_run_ctlfile(fast.radiolink.sat2sat
                      artscomponents/radiolink/TestRadioLink.arts)
arts_test_run_ctlfile(fast.radiolink.sat2surface
                      artscomponents/radiolink/TestRadioLink2.arts)
arts_test_run_ctlfile(slow.radiolink.occultation
                      artscomponents/radiolink/TestRadioOccultation.arts)

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

arts_test_run_ctlfile(slow.wfuns artscomponents/wfuns/TestWfuns.arts)

#arts_test_run_ctlfile(slow.zeeman artscomponents/zeeman/TestZeeman.arts)

#arts_test_run_ctlfile(slow.disort artscomponents/disort/TestDISORT.arts)

###################
### Instruments ###
arts_test_run_ctlfile(fast.odinsmr instruments/odinsmr/TestOdinSMR.arts)
arts_test_run_ctlfile(slow.odinsmr.1D instruments/odinsmr/TestOdinSMR_1D.arts)

arts_test_run_ctlfile(fast.amsub instruments/amsu/TestAMSUB_fast.arts)
arts_test_run_ctlfile(slow.amsub instruments/amsu/TestAMSUB.arts)

arts_test_run_ctlfile(fast.mhs instruments/mhs/TestMHS_fast.arts)
arts_test_run_ctlfile(slow.mhs instruments/mhs/TestMHS.arts)

arts_test_run_ctlfile(slow.hirs instruments/hirs/TestHIRS.arts)

#arts_test_run_ctlfile(slow.avhrr instruments/avhrr/TestAVHRR.arts)

