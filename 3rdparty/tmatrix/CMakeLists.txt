if (ENABLE_FORTRAN AND NOT NO_TMATRIX)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tmatrix.f)
    set (ENABLE_TMATRIX true)
    message (STATUS "TMATRIX enabled")
  else()
    message (STATUS "TMATRIX disabled")
  endif()

  if (ENABLE_TMATRIX)
    add_library (tmatrix STATIC
      tmatrix.f
      )
    set_target_properties (tmatrix PROPERTIES
                           COMPILE_FLAGS "-fdefault-double-8 -fdefault-integer-8 -fdefault-real-8")
  endif()

endif()

