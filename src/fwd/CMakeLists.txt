add_library(fwdlbl STATIC
  fwd_lbl.cpp
  fwd_lbl_mtckd_voigt.cpp
)
target_link_libraries(fwdlbl PUBLIC matpack Faddeeva)

add_library(fwdcia STATIC
  fwd_cia.cpp
  ../cia.cc # ## FIXME: Include this elsewise, or do manual computations ###
)
target_link_libraries(fwdcia PUBLIC matpack)
target_include_directories(fwdcia PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

# add_dependencies(fwdcia make_tokval)  ### FIXME: remove when messages.h is gone? ###
add_library(fwdpredef STATIC
  fwd_predef.cpp
)
target_link_libraries(fwdpredef PUBLIC matpack predef)

add_library(fwdabs STATIC
  fwd_abs.cpp
)
target_link_libraries(fwdabs PUBLIC fwdlbl fwdcia fwdpredef)

add_library(fwdprofile STATIC
  fwd_profile.cpp
)
target_link_libraries(fwdprofile PUBLIC fwdabs)

# ## Main interface, everyone should just "include <fwd.h>" and link to fwd"
add_library(fwd INTERFACE)
target_link_libraries(fwd INTERFACE fwdabs fwdprofile)
