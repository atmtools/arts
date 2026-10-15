# ARTS's RT4 (src/core/polradtran/rt4) on Evans' RT4 script runtestc, with his
# optics (3rdparty/polradtran/cl340d14.dda), against his expected output
add_executable(polradtran-rt4-arts polradtran-rt4-arts.cpp)
target_link_libraries(polradtran-rt4-arts arts_rt4)
add_test(NAME "cpp.fast.polradtran-rt4-arts" COMMAND polradtran-rt4-arts ${ARTS_SOURCE_DIR}/3rdparty/polradtran)
add_dependencies(check-deps polradtran-rt4-arts)
