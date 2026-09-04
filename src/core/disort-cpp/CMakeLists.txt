if (ENABLE_ARTS_LGPL)
  add_library(disort-cpp-eigen STATIC disort-eigen.lgpl.cpp)
else (ENABLE_ARTS_LGPL)
  add_library(disort-cpp-eigen STATIC disort-eigen.gpl.cpp)
endif ()
target_link_libraries(disort-cpp-eigen PUBLIC matpack)

add_library(disort-cpp-common STATIC common.cpp)
target_link_libraries(disort-cpp-common PUBLIC matpack legendre)
target_include_directories(disort-cpp-common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(disort-cpp STATIC disort.cpp disort-brdf.cpp)

target_link_libraries(disort-cpp PUBLIC matpack legendre time_report operators disort-cpp-eigen disort-cpp-common)
target_include_directories(disort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# VDISORT CHANGE BEGIN: polarized four-Stokes sibling of disort::main_data.
add_library(vdisort-cpp STATIC vdisort.cpp vdisort-brdf.cpp)
target_link_libraries(vdisort-cpp PUBLIC matpack legendre time_report operators rtepack disort-cpp-common)
target_include_directories(vdisort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# VDISORT CHANGE END

add_subdirectory(test)
