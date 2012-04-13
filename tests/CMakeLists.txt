
add_custom_command(
  OUTPUT clean-tests-cmd
  COMMAND rm -f "Test*.xml" "Test*.rep" "Test*.bin"
  )
add_custom_target (clean-tests DEPENDS clean-tests-cmd)

