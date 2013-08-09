if (ENABLE_FORTRAN AND NOT NO_TMATRIX)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tmd.arts.f)
    set (ENABLE_TMATRIX true)
    message (STATUS "TMATRIX enabled")
  else()
    message (STATUS "TMATRIX disabled")
  endif()

  if (ENABLE_TMATRIX)
    add_library (tmatrix STATIC
      tmd.arts.f
      ampld.arts.f
      lpd.f
      avg_tmatrix.f
      )
    set_target_properties (tmatrix PROPERTIES
      COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")

    add_executable(tmatrix_tmd
      tmd.lp.f
      lpd.f)

    add_executable(tmatrix_ampld
      ampld.lp.f
      lpd.f)
  endif()

endif()

