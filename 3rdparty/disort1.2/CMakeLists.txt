if (ENABLE_FORTRAN AND NOT NO_DISORT)
  set (ENABLE_DISORT true)

  add_library (disort STATIC
    D1MACH.f
    DISORT.f
    ErrPack.f
    LINPAK.f
    R1MACH.f
    )
  set_target_properties (disort PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")

  add_executable(test_disort DISOTEST.f)
  target_link_libraries(test_disort disort)
  set_target_properties (test_disort PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set (ENABLE_DISORT false)
endif()

if (ENABLE_DISORT)
  message (STATUS "Disort enabled (use -DNO_DISORT=1 to disable)")
elseif(NO_DISORT)
  message (STATUS "Disort disabled")
else()
  message (STATUS "Disort disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

