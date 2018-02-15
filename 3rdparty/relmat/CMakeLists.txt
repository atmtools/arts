if (ENABLE_FORTRAN AND NOT NO_RELMAT)
  set (ENABLE_RELMAT true)
  if(FASTWIGNER)
    add_library (relmat STATIC
      ../wigner/fastwigxj/src/ffastwigxj.f
      module_common_var.F90
      module_error.F90
      module_maths.F90
      module_molecSp.F90
      module_read.F90
      module_phsub.F90
      module_LLS.F90
      module_linemixing.F90
      arts_interface.F90)
    target_link_libraries( relmat wigner ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} wigner )
    add_definitions("-DUSE_EXTERNAL_WIGNER=1")
  else()
    add_library (relmat STATIC
      module_common_var.F90
      module_error.F90
      module_maths.F90
      module_molecSp.F90
      module_read.F90
      module_phsub.F90
      module_LLS.F90
      module_linemixing.F90
      arts_interface.F90)
    target_link_libraries( relmat wigner ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
    add_definitions("-DUSE_EXTERNAL_WIGNER=0")
  endif()
  set_target_properties (relmat PROPERTIES COMPILE_FLAGS "${FORTRAN_EXTRA_FLAGS}")
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
