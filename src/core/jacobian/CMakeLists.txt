add_library(jacobian STATIC
    jacobian.cc
    covariance_matrix.cc
    retrieval_target.cc
    jac_polyfit.cc
    jac_rel.cc
    jac_log.cc
)

target_link_libraries(jacobian PUBLIC atm subsurface surface sensor scattering lbl)
target_include_directories(jacobian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
