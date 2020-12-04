project (arts_docs)
find_package (LATEX)

if (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)

  add_custom_target (auto_version_tex DEPENDS auto_version.tex)
  add_dependencies (auto_version_tex UpdateAutoVersion)

  add_custom_command (
    OUTPUT auto_version.tex
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${ARTS_BINARY_DIR}/auto_version.txt auto_version.tex
    DEPENDS ${ARTS_BINARY_DIR}/auto_version.txt
    )

  include (ArtsBuildTexDoc)

  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/arts_user.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/common.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/symbol_defs.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/concept.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/in_and_out.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/atmosphere.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/rte_basics.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/complete_calcs.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/rte.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/faraday.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/transmission.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/cloudradar.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/absorption.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/refr_index.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/ppath.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/winds.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/surface.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/sensor.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/batch.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/clouds.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/scattering.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/wfuns.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/gridded_field_gridnames.txt")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_user" "${TEXFILES}" "${FIGFILES}")
  
  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/arts_developer.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/common.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/symbol_defs.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/development.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/agendas.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/matrix_vector.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/gridded_fields.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/interpolation.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/integration.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/lin_alg.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/cpp_api.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/debug_gui.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_developer" "${TEXFILES}" "${FIGFILES}")
  
  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/arts_theory.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/common.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/symbol_defs.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/formalismn.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/absorption_theory.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/rte_theory.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/polarization.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/clouds_theory.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/ppath_theory.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/doit.tex"
                      "${CMAKE_CURRENT_SOURCE_DIR}/montecarlo.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_theory" "${TEXFILES}" "${FIGFILES}")
  
  add_custom_target (arts_docs
    DEPENDS arts_user.pdf arts_developer.pdf arts_theory.pdf)

else (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)

  if (NOT PDFLATEX_COMPILER)
    message(STATUS "User guide disabled due to missing pdflatex")
  endif (NOT PDFLATEX_COMPILER)

  if (NOT BIBTEX_COMPILER)
    message(STATUS "User guide disabled due to missing bibtex")
  endif (NOT BIBTEX_COMPILER)

  if (NOT MAKEINDEX_COMPILER)
    message(STATUS "User guide disabled due to missing makeindex")
  endif (NOT MAKEINDEX_COMPILER)

  set (ARTS_FEATURE_WARN true)

  mark_as_advanced (ARTS_FEATURE_WARN)
endif (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)

add_test(
  NAME arts.doc.uguide.refs
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_refs.py ${ARTS_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_custom_target(check-uguide
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} DEPENDS arts)

