set(files
  doc.rst
  doc.theory.rst
  doc.user.rst
)

foreach(fil ${files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${fil} ${CMAKE_BINARY_DIR}/python/doc/source/${fil} COPYONLY)
endforeach(fil ${files})
