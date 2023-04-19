if (ENABLE_FORTRAN AND NOT NO_FASTEM)
  set (ENABLE_FASTEM true)

  add_library (fastem STATIC
    parkind1.F90
    mod_rttov_fastem5_coef.F90
    rttov_fastem5.F90
    )
  set_target_properties (fastem PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set(ENABLE_FASTEM false)
endif()

if (ENABLE_FASTEM)
  message (STATUS "FASTEM enabled (use -DNO_FASTEM=1 to disable)")
elseif(NO_FASTEM)
  message (STATUS "FASTEM disabled")
else()
  message (STATUS "FASTEM disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

