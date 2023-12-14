if (ENABLE_FORTRAN AND NOT NO_SURFEMOCEAN)
  set (ENABLE_SURFEMOCEAN true)

  add_library (surfemocean STATIC
    parkind1.F90
    rttov_const.F90
    mod_rttov_fastem5_coef.F90
    rttov_surfem_ocean_coef_mod.F90
    rttov_surfem_ocean_mod.F90
    surfem_ocean.F90
    )
  set_target_properties (surfemocean PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set(ENABLE_SURFEMOCEAN false)
endif()

if (ENABLE_SURFEMOCEAN)
  message (STATUS "SURFEMOCEAN enabled (use -DNO_SURFEMOCEAN=1 to disable)")
elseif(NO_SURFEMOCEAN)
  message (STATUS "SURFEMOCEAN disabled")
else()
  message (STATUS "SURFEMOCEAN disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

