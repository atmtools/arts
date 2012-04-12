
add_custom_command(
  OUTPUT clean-tests-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-tests DEPENDS clean-tests-cmd)

arts_test_run_ctlfile(fast.radiolink      RadioLink/TestRadioLink.arts)
arts_test_run_ctlfile(fast.amsub          AMSU/TestAMSUB_fast.arts)
arts_test_run_ctlfile(fast.odinsmr        OdinSMR/TestOdinSMR.arts)

arts_test_run_ctlfile(slow.amsub AMSU/TestAMSUB.arts)

arts_test_run_ctlfile(fast.mhs MHS/TestMHS_fast.arts)
arts_test_run_ctlfile(slow.mhs MHS/TestMHS.arts)

arts_test_run_ctlfile(slow.wfuns Wfuns/TestWfuns.arts)

