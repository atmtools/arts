if (ENABLE_FORTRAN AND NOT NO_DISORT)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/DISORT.f)
    set (ENABLE_DISORT true)
    message (STATUS "DISORT enabled")
  else()
    message (STATUS "DISORT disabled")
  endif()

  if (ENABLE_DISORT)
    add_library (disort STATIC
      D1MACH.f
      DISORT.f
      DISOTEST.f
      ErrPack.f
      LINPAK.f
      R1MACH.f
      )
    set_target_properties (disort PROPERTIES
      COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
  endif()

endif()

