add_library(fwd STATIC
  lbl.cpp
  lbl_mtckd_voigt.cpp
)

target_link_libraries(fwd PUBLIC matpack Faddeeva)
