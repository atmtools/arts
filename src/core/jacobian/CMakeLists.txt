add_library(jacobian STATIC
    covariance_matrix.cc
    jac_log.cc
    jac_logrel.cc
    jac_polyfit.cc
    jac_rel.cc
    jac_rh.cc
    jacobian.cc
    retrieval_target.cc
)

target_link_libraries(jacobian PUBLIC atm subsurface surface sensor scattering lbl)
target_include_directories(jacobian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
