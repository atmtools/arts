# Evans' PolRadTran benchmark scripts (3rdparty/polradtran), unchanged, with
# his original programs rt3-evans, rt4-evans and scatcnv-evans
if(NOT ENABLE_RT3 AND NOT ENABLE_RT4)
  return()
endif()

add_executable(polradtran-scripts polradtran-scripts.cpp)
add_dependencies(check-deps polradtran-scripts)

set(polradtran_dir ${ARTS_SOURCE_DIR}/3rdparty/polradtran)
macro(polradtran_script script)
  add_test(NAME "cpp.fast.polradtran-${script}"
           COMMAND polradtran-scripts ${polradtran_dir}/${script} ${CMAKE_CURRENT_BINARY_DIR}/${script} ${ARGN})
endmacro()
if(ENABLE_RT3)
  polradtran_script(runmietest rt3=$<TARGET_FILE:rt3-evans>)
  polradtran_script(runtesta rt3=$<TARGET_FILE:rt3-evans>)
  add_dependencies(check-deps rt3-evans)
endif()
if(ENABLE_RT4)
  polradtran_script(runtestr scatcnv=$<TARGET_FILE:scatcnv-evans> rt4=$<TARGET_FILE:rt4-evans>)
  polradtran_script(runtestc rt4=$<TARGET_FILE:rt4-evans>)
  add_dependencies(check-deps rt4-evans scatcnv-evans)
endif()

# ARTS's RT4 (src/core/polradtran/rt4) on Evans' RT4 scripts, with his optics (cl340d14.dda, and scatcnv's conversion)
if(ENABLE_RT4)
  add_executable(polradtran-rt4-arts polradtran-rt4-arts.cpp)
  target_link_libraries(polradtran-rt4-arts arts_rt4)
  add_test(NAME "cpp.fast.polradtran-rt4-arts"
           COMMAND polradtran-rt4-arts ${polradtran_dir} $<TARGET_FILE:scatcnv-evans> ${CMAKE_CURRENT_BINARY_DIR}/rt4-arts)
  add_dependencies(check-deps polradtran-rt4-arts)
endif()
