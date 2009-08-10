project (uguide)
find_package (LATEX)

if (PDFLATEX_COMPILER AND BIBTEX_COMPILER)

  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/auto_version.tex.in
                  ${CMAKE_CURRENT_BINARY_DIR}/auto_version.tex
                  @ONLY IMMEDIATE)

  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/uguide.pdf
    COMMAND TEXINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${PDFLATEX_COMPILER} uguide
            && makeindex uguide.idx
            && BSTINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: BIBINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${BIBTEX_COMPILER} uguide
            && TEXINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${PDFLATEX_COMPILER} uguide
            && TEXINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${PDFLATEX_COMPILER} uguide
            && TEXINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${PDFLATEX_COMPILER} uguide
            || rm -f uguide.pdf
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.tex ${CMAKE_CURRENT_BINARY_DIR}/auto_version.tex
  )

  add_custom_target (uguide ALL echo DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/uguide.pdf)

endif (PDFLATEX_COMPILER AND BIBTEX_COMPILER)

set_directory_properties (PROPERTIES ADDITIONAL_CLEAN_FILES *.aux)

