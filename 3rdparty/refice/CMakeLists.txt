if (ENABLE_FORTRAN AND NOT NO_REFICE)
  set (ENABLE_REFICE true)

  add_library (refice STATIC
    REFICE.f
    )
  set_target_properties (refice PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set(ENABLE_REFICE false)
endif()

if (ENABLE_REFICE)
  message (STATUS "Refice enabled (use -DNO_REFICE=1 to disable)")
elseif(NO_REFICE)
  message (STATUS "Refice disabled")
else()
  message (STATUS "Refice disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

