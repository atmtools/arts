macro (ARTS_ADD_TEX_DOC TARGET TEXFILES FIGFILES)

  set (PDFLATEX_OUTPUTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}")
  set (PDFLATEX_OPTIONS "-interaction;nonstopmode;-halt-on-error;-shell-escape")
  set (PDFLATEX_LOGGING ">>;${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.log")

  set (LOG_OUTPUT "||" "\(cat" "${TARGET}.log" "&&" "exit" "1\)")

  get_filename_component (MAKEINDEX_PATH "${MAKEINDEX_COMPILER}" PATH)

  set (PDFLATEX_COMMAND "TEXINPUTS=.:..:${CMAKE_CURRENT_SOURCE_DIR}:"
    "save_size=15000" ${PDFLATEX_COMPILER} ${PDFLATEX_OPTIONS} ${TARGET}
    ${PDFLATEX_LOGGING} ${LOG_OUTPUT})
  set (MAKEINDEX_COMMAND "PATH=$ENV{PATH}:${MAKEINDEX_PATH}"
    ${MAKEINDEX_COMPILER} "-q" "${TARGET}.idx" ${PDFLATEX_LOGGING})
  set (BIBTEX_COMMAND "BSTINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}:"
    "BIBINPUTS=.:${CMAKE_CURRENT_SOURCE_DIR}:" "save_size=15000"
    ${BIBTEX_COMPILER} ${TARGET} ${PDFLATEX_LOGGING})

  if (${TARGET} STREQUAL "arts_user")
    set (OTHER_DOC1 "arts_developer")
    set (OTHER_DOC2 "arts_theory")
  elseif (${TARGET} STREQUAL "arts_theory")
    set (OTHER_DOC1 "arts_user")
    set (OTHER_DOC2 "arts_developer")
  elseif (${TARGET} STREQUAL "arts_developer")
    set (OTHER_DOC1 "arts_user")
    set (OTHER_DOC2 "arts_theory")
  endif ()

  add_custom_command (
    OUTPUT ${TARGET}.stamp
    COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET}
    COMMAND ${CMAKE_COMMAND} -E touch ${TARGET}.stamp
  )

  add_custom_command (
    OUTPUT ${TARGET}.stage1
    COMMAND ${CMAKE_COMMAND} -E remove -f ${TARGET}.log
    COMMAND ${PDFLATEX_COMMAND}
    COMMAND ${MAKEINDEX_COMMAND}
    COMMAND ${BIBTEX_COMMAND}
    COMMAND ${PDFLATEX_COMMAND}
    COMMAND cp "${PDFLATEX_OUTPUTDIR}/*.aux" "${PDFLATEX_OUTPUTDIR}/../${OTHER_DOC1}/" 
    COMMAND cp "${PDFLATEX_OUTPUTDIR}/*.aux" "${PDFLATEX_OUTPUTDIR}/../${OTHER_DOC2}/" 
    COMMAND ${CMAKE_COMMAND} -E touch ${PDFLATEX_OUTPUTDIR}.stage1
    WORKING_DIRECTORY ${PDFLATEX_OUTPUTDIR}
    DEPENDS ${TEXFILES} ${FIGFILES}
            ${CMAKE_CURRENT_BINARY_DIR}/auto_version.tex
            ${TARGET}.stamp ${OTHER_DOC1}.stamp ${OTHER_DOC2}.stamp
  )

  add_custom_command (
    OUTPUT ${TARGET}.pdf
    COMMAND ${CMAKE_COMMAND} -E remove -f ${TARGET}.log
    COMMAND ${PDFLATEX_COMMAND}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PDFLATEX_OUTPUTDIR}/${TARGET}.pdf" "${CMAKE_CURRENT_BINARY_DIR}"
    WORKING_DIRECTORY ${PDFLATEX_OUTPUTDIR}
    DEPENDS ${TARGET}.stage1 ${OTHER_DOC1}.stage1 ${OTHER_DOC2}.stage1
  )

  add_custom_target (${TARGET}-stamp DEPENDS ${TARGET}.stamp)
  add_custom_target (${TARGET}-stage1 DEPENDS ${TARGET}.stage1)

  add_custom_target (${TARGET} ALL DEPENDS ${TARGET}.pdf)

  add_dependencies (${TARGET}-stage1 auto_version_tex
                    ${TARGET}-stamp ${OTHER_DOC1}-stamp ${OTHER_DOC2}-stamp)

  add_dependencies (${TARGET}
                    ${TARGET}-stage1 ${OTHER_DOC1}-stage1 ${OTHER_DOC2}-stage1)

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.pdf
           DESTINATION share/doc/arts/uguide)

  get_directory_property (TEXCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
  set (TEXCLEANFILES ${TEXCLEANFILES} ${TARGET}.log ${PDFLATEX_OUTPUTDIR})
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                            "${TEXCLEANFILES}")

endmacro (ARTS_ADD_TEX_DOC TARGET DEPENDENCIES)


macro (ARTS_ADD_SIMPLE_TEX_DOC TARGET TEXFILES FIGFILES)

  set (PDFLATEX_OUTPUTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}")
  set (PDFLATEX_OPTIONS "-interaction;nonstopmode;-halt-on-error")
  set (PDFLATEX_LOGGING ">>;${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.log")

  set (LOG_OUTPUT "||" "\(cat" "${TARGET}.log" "&&" "exit" "1\)")

  set (PDFLATEX_COMMAND "TEXINPUTS=.:..:${CMAKE_CURRENT_SOURCE_DIR}:"
    "save_size=15000" ${PDFLATEX_COMPILER} ${PDFLATEX_OPTIONS} ${TARGET}
    ${PDFLATEX_LOGGING} ${LOG_OUTPUT})

  add_custom_command (
    OUTPUT ${TARGET}.stamp
    COMMAND ${CMAKE_COMMAND} -E make_directory ${TARGET}
    COMMAND ${CMAKE_COMMAND} -E touch ${TARGET}.stamp
    DEPENDS ${TEXFILES} ${FIGFILES}
  )

  add_custom_command (
    OUTPUT ${TARGET}.pdf
    COMMAND ${CMAKE_COMMAND} -E remove -f ${TARGET}.log
    COMMAND ${PDFLATEX_COMMAND}
    COMMAND ${PDFLATEX_COMMAND}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${PDFLATEX_OUTPUTDIR}/${TARGET}.pdf" "${CMAKE_CURRENT_BINARY_DIR}"
    WORKING_DIRECTORY ${PDFLATEX_OUTPUTDIR}
    DEPENDS ${TARGET}.stamp
  )

  add_custom_target (${TARGET} ALL DEPENDS ${TARGET}.pdf)

  get_directory_property (TEXCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
  set (TEXCLEANFILES ${TEXCLEANFILES} ${TARGET}.log ${PDFLATEX_OUTPUTDIR})
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                            "${TEXCLEANFILES}")

endmacro (ARTS_ADD_SIMPLE_TEX_DOC TARGET DEPENDENCIES)

