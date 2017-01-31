if (ENABLE_FORTRAN AND NOT NO_RELMAT)
  set (ENABLE_RELMAT true)

  add_library (relmat STATIC
    module_common_var.F90
    module_error.F90
    module_maths.F90
    module_molecSp.F90
    module_read.F90
    module_phsub.F90
    module_linemixing.F90
    arts_interface.F90
    )
  target_link_libraries( relmat ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
  set_target_properties (relmat PROPERTIES
    COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
else()
  set(ENABLE_RELMAT false)
endif()

if (ENABLE_RELMAT)
  message (STATUS "RELMAT enabled (use -DNO_RELMAT=1 to disable)")
elseif(NO_RELMAT)
  message (STATUS "RELMAT disabled")
else()
  message (STATUS "RELMAT disabled (use -DENABLE_FORTRAN=1 to enable)")
endif()

