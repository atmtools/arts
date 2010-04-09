project (arts_docs)
find_package (LATEX)

if (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)

  include ("ArtsBuildTexDoc")

  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/auto_version.tex.in
                  ${CMAKE_CURRENT_BINARY_DIR}/auto_version.tex @ONLY)

  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_user" "${TEXFILES}" "${FIGFILES}")
  
  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_developer" "${TEXFILES}" "${FIGFILES}")
  
  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/Figs/*/*.pdf")
  ARTS_ADD_TEX_DOC ("arts_theory" "${TEXFILES}" "${FIGFILES}")
  
  add_custom_target (arts_docs ALL DEPENDS arts_user.pdf arts_developer.pdf arts_theory.pdf)

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

