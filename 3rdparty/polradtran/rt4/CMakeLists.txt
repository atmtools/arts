if (ENABLE_FORTRAN AND NOT NO_RT4)
  set (ENABLE_RT4 true)

  add_library (rt4 STATIC
    radtran4.f
    radmat.f
    radutil4.f
    radintg4.f
    radscat4.f
    )
  set_target_properties (rt4 PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set(ENABLE_RT4 false)
endif()

if (ENABLE_RT4)
  message (STATUS "RT4 enabled (use -DNO_RT4=1 to disable)")
elseif(NO_RT4)
  message (STATUS "RT4 disabled")
else()
  message (STATUS "RT4 disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

