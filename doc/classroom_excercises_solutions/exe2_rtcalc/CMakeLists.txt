project (arts_docs)
find_package (LATEX)

if (PDFLATEX_COMPILER)
  include (ArtsBuildTexDoc)

  file (GLOB TEXFILES "${CMAKE_CURRENT_SOURCE_DIR}/exe2_solution.tex")
  file (GLOB FIGFILES "${CMAKE_CURRENT_SOURCE_DIR}/plots/*.pdf")

  ARTS_ADD_SIMPLE_TEX_DOC ("exe2_solution" "${TEXFILES}" "${FIGFILES}")

else (PDFLATEX_COMPILER)

  if (NOT PDFLATEX_COMPILER)
    message(STATUS "User guide disabled due to missing pdflatex")
  endif (NOT PDFLATEX_COMPILER)

endif (PDFLATEX_COMPILER)

