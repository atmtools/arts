if (ENABLE_FORTRAN AND NOT NO_TMATRIX)
  set (ENABLE_TMATRIX true)

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

  if (FORTRAN_COMPILER MATCHES "ifort.*")
    add_executable(tmatrix_tmq
      tmq.lp.f
      lpq.f)

    add_executable(tmatrix_amplq
      amplq.lp.f
      lpq.f)
  endif()

else()
  set(ENABLE_TMATRIX false)
endif()

if (ENABLE_TMATRIX)
  message (STATUS "Tmatrix enabled (use -DNO_TMATRIX=1 to disable)")
elseif(NO_TMATRIX)
  message (STATUS "Tmatrix disabled")
else()
  message (STATUS "Tmatrix disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

