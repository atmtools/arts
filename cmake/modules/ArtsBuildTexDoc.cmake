macro (ARTS_ADD_TEX_DOC TARGET TEXFILES FIGFILES)

  set (PDFLATEX_OPTIONS "-interaction;batchmode;-halt-on-error")
  set (PDFLATEX_LOGGING ">>;${TARGET}.log")
  set (PDFLATEX_COMMAND "TEXINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}:"
    "save_size=15000" ${PDFLATEX_COMPILER} ${PDFLATEX_OPTIONS} ${TARGET}
    ${PDFLATEX_LOGGING})

  get_filename_component(MAKEINDEX_PATH "${MAKEINDEX_COMPILER}" PATH)

  add_custom_command (
    OUTPUT ${TARGET}.pdf
    COMMAND ${CMAKE_COMMAND} -E remove -f ${TARGET}.log
    COMMAND ${PDFLATEX_COMMAND}
            && PATH=$PATH:${MAKEINDEX_PATH} ${MAKEINDEX_COMPILER} -q ${TARGET}.idx ${PDFLATEX_LOGGING}
            && BSTINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: BIBINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}: save_size=15000
              ${BIBTEX_COMPILER} ${TARGET} ${PDFLATEX_LOGGING}
            && ${PDFLATEX_COMMAND}
            && ${PDFLATEX_COMMAND}
            && ${PDFLATEX_COMMAND} || \(cat ${TARGET}.log && exit 1\)
    DEPENDS ${TEXFILES} ${FIGFILES} ${CMAKE_CURRENT_BINARY_DIR}/auto_version.tex
  )

  add_custom_target (${TARGET} ALL DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.pdf)

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.pdf
           DESTINATION share/doc/arts)

  set (TEXCLEANFILES "")
  get_directory_property (TEXCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
  foreach (TEXFILE ${TEXFILES})
    get_filename_component (BASETEXFILE ${TEXFILE} NAME_WE)
    set (TEXCLEANFILES ${TEXCLEANFILES} ${BASETEXFILE}.aux)
  endforeach (TEXFILE)

  set (TEXCLEANFILES ${TEXCLEANFILES}
       ${TARGET}.bbl ${TARGET}.blg ${TARGET}.idx ${TARGET}.ilg
       ${TARGET}.ind ${TARGET}.log ${TARGET}.out ${TARGET}.toc)

  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                            "${TEXCLEANFILES}")

endmacro (ARTS_ADD_TEX_DOC TARGET DEPENDENCIES)

