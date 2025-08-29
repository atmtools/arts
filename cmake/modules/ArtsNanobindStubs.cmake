function (ARTS_ADD_CPP_STUBS)

  foreach (MODULENAME IN LISTS ARGN)
    nanobind_add_stub(
      pyarts_${MODULENAME}_cpp_stub
      MODULE arts.${MODULENAME}
      OUTPUT ${ARTS_BINARY_DIR}/python/src/pyarts3/arts/${MODULENAME}.pyi
      PYTHON_PATH ${ARTS_BINARY_DIR}/python/src/pyarts3
      DEPENDS pyarts_cpp
    )
    list(APPEND deplist "pyarts_${MODULENAME}_cpp_stub")
  endforeach()

  nanobind_add_stub(
    pyarts_cpp_stub
    MODULE arts
    OUTPUT ${ARTS_BINARY_DIR}/python/src/pyarts3/arts/__init__.pyi
    PYTHON_PATH ${ARTS_BINARY_DIR}/python/src/pyarts3
    DEPENDS pyarts_cpp "${deplist}"
  )
  set_property(
    TARGET pyarts_cpp_stub
    APPEND
    PROPERTY ADDITIONAL_CLEAN_FILES ${ARTS_BINARY_DIR}/python/src/pyarts3/arts
  )

endfunction()
