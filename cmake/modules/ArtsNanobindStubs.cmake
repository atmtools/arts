function (ARTS_ADD_CPP_STUBS)

  if (ENABLE_PYARTS_STUBS)
    foreach (MODULENAME IN LISTS ARGN)
      nanobind_add_stub(
        pyarts_${MODULENAME}_cpp_stub
        MODULE pyarts3.arts.${MODULENAME}
        OUTPUT ${ARTS_BINARY_DIR}/python/src/pyarts3/arts/${MODULENAME}.pyi
        PYTHON_PATH ${ARTS_BINARY_DIR}/python/src
        DEPENDS pyarts_cpp
      )
      list(APPEND deplist "pyarts_${MODULENAME}_cpp_stub")
    endforeach()

    nanobind_add_stub(
      pyarts_cpp_stub
      MODULE pyarts3.arts
      OUTPUT ${ARTS_BINARY_DIR}/python/src/pyarts3/arts/__init__.pyi
      PYTHON_PATH ${ARTS_BINARY_DIR}/python/src
      DEPENDS pyarts_cpp "${deplist}"
    )
    set_property(
      TARGET pyarts_cpp_stub
      APPEND
      PROPERTY ADDITIONAL_CLEAN_FILES ${ARTS_BINARY_DIR}/python/src/pyarts3/arts
    )
  else()
    add_custom_target(pyarts_cpp_stub)
  endif()

endfunction()
