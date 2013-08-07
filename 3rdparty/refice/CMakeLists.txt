if (ENABLE_FORTRAN AND NOT NO_REFICE)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/REFICE.f)
    set (ENABLE_REFICE true)
    message (STATUS "REFICE enabled")
  else()
    message (STATUS "REFICE disabled")
  endif()

  if (ENABLE_REFICE)
    add_library (refice STATIC
      REFICE.f
      )
    set_target_properties (refice PROPERTIES
      COMPILE_FLAGS "-fdefault-double-8 -fdefault-integer-8 -fdefault-real-8")

  endif()

endif()

