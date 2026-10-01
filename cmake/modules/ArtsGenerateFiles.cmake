# Code generators that only touch the files they actually change
#
# Most automatic source files of ARTS are written by small programs that are
# built and run as part of the build.  Such a generator has to be re-run
# whenever one of the libraries it links to changes, and rewriting an output
# with the very same text is enough to make the build system recompile
# everything that includes it.
#
# arts_generate() therefore lets the generator write into a scratch directory of
# its own and then copies the result next to the other build files with
# "copy_if_different".  A file that comes out exactly as it was keeps its time
# stamp and nothing downstream is recompiled.
#
#   arts_generate(<generator>
#                 OUTPUT   <file>...  # generated files, named relative to the
#                                     # current binary directory
#                 [ARGS    <arg>...]  # arguments passed to the generator
#                 [DEPENDS <dep>...]  # further inputs of the generator
#                 [COMMENT <text>])
#
# The whole scratch directory is copied rather than the listed files alone,
# because some generators also write files whose names CMake does not know
# (make_enums writes one header per enumeration).
function(arts_generate generator)
  cmake_parse_arguments(PARSE_ARGV 1 GEN "" "COMMENT" "OUTPUT;ARGS;DEPENDS")

  # One scratch directory per generator: several of them write into the same
  # binary directory in parallel and must not see each other's unfinished work
  set(scratch ${CMAKE_CURRENT_BINARY_DIR}/tmp/${generator})

  set(scratch_output "")
  foreach(file ${GEN_OUTPUT})
    list(APPEND scratch_output ${scratch}/${file})
  endforeach()

  add_custom_command(
    OUTPUT ${scratch_output}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${scratch}
    COMMAND ${CMAKE_COMMAND} -E chdir ${scratch} $<TARGET_FILE:${generator}> ${GEN_ARGS}
    DEPENDS ${generator} ${GEN_DEPENDS}
    COMMENT "${GEN_COMMENT}"
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${GEN_OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${scratch} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${scratch_output}
    COMMENT "${GEN_COMMENT} (copying the files that changed)"
    VERBATIM
  )
endfunction()
