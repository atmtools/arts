if (ENABLE_ARTS_LGPL)
  add_library(disort-cpp-eigen STATIC disort-eigen.lgpl.cpp)
else (ENABLE_ARTS_LGPL)
  add_library(disort-cpp-eigen STATIC disort-eigen.gpl.cpp)
endif ()
target_link_libraries(disort-cpp-eigen PUBLIC matpack)

add_library(disort-cpp STATIC disort.cpp)

target_link_libraries(disort-cpp PUBLIC matpack legendre time_report operators disort-cpp-eigen)
target_include_directories(disort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# VDISORT CHANGE BEGIN: polarized four-Stokes sibling of disort::main_data.
add_library(vdisort-cpp STATIC vdisort.cpp)
target_link_libraries(vdisort-cpp PUBLIC matpack legendre time_report operators rtepack)
target_include_directories(vdisort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# VDISORT CHANGE END

add_subdirectory(test)
