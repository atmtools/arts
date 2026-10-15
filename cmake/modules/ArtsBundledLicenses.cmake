# Licences of bundled code: the code from other projects that a build of ARTS
# compiles in.  Every distribution of ARTS, source or binary, must carry the
# licence notices of the bundled code it contains.  See
# doc/arts/dev.licenses.rst for the policy.
#
#   arts_add_bundled_license(<name> SPDX <expression> FILES <file>...)
#
#     Register a bundled component.  Call it where the component is built, in
#     the branch of the CMake logic that builds or links it, so that exactly
#     the components of the current configuration are registered.
#     <expression> is the SPDX licence expression of the component, or a
#     LicenseRef-<id> for a licence without an SPDX identifier.  The FILES
#     are the component's licence and notice files, relative to the calling
#     directory; configuring fails if one is missing.
#
#   arts_write_bundled_licenses(<directory> <expression-variable>)
#
#     Copy the files of every registered component to <directory>/<name>/,
#     write the index <directory>/THIRD_PARTY.txt, and set
#     <expression-variable> to the licence expression of this build: the
#     licence of ARTS AND the licences of the registered components.  Call it
#     after every component is registered (the Python package does this in
#     python/CMakeLists.txt).  The directory is recreated, so that components
#     of an earlier configuration do not linger.
#
#   arts_write_bundled_license_source(<file.cpp>)
#
#     Write a C++ source that compiles the licence texts of every registered
#     component into ARTS (src/core/licenses/bundled_licenses.h), so that
#     every binary carries them.  The file is only rewritten when its content
#     changes, and editing a licence file reconfigures.
#
#   arts_bundled_license_expression(<variable>)
#
#     Set <variable> to the licence expression of this build.

set(ARTS_LICENSE_EXPRESSION "LGPL-3.0-or-later OR GPL-3.0-or-later")

function(arts_add_bundled_license name)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "SPDX" "FILES")
  if(NOT ARG_SPDX OR NOT ARG_FILES OR ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "arts_add_bundled_license(${name}) takes SPDX <expression> FILES <file>...")
  endif()
  set(files)
  foreach(file IN LISTS ARG_FILES)
    cmake_path(ABSOLUTE_PATH file BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    if(NOT EXISTS "${file}")
      message(FATAL_ERROR "The licence file ${file} of the bundled component ${name} does not exist")
    endif()
    list(APPEND files "${file}")
  endforeach()
  get_property(names GLOBAL PROPERTY ARTS_BUNDLED)
  if("${name}" IN_LIST names)
    message(FATAL_ERROR "The bundled component ${name} is registered twice")
  endif()
  set_property(GLOBAL APPEND PROPERTY ARTS_BUNDLED "${name}")
  set_property(GLOBAL PROPERTY ARTS_BUNDLED_${name}_SPDX "${ARG_SPDX}")
  set_property(GLOBAL PROPERTY ARTS_BUNDLED_${name}_FILES "${files}")
endfunction()

function(arts_write_bundled_licenses directory expression_variable)
  file(REMOVE_RECURSE "${directory}")
  get_property(names GLOBAL PROPERTY ARTS_BUNDLED)
  set(index
      "Code from other projects compiled into this build of ARTS, with the SPDX\n"
      "licence expression and the licence files (in <name>/) of each component.\n"
      "ARTS itself is ${ARTS_LICENSE_EXPRESSION} (LICENSE.txt).\n\n")
  foreach(name IN LISTS names)
    get_property(spdx GLOBAL PROPERTY ARTS_BUNDLED_${name}_SPDX)
    get_property(files GLOBAL PROPERTY ARTS_BUNDLED_${name}_FILES)
    set(filenames)
    foreach(file IN LISTS files)
      cmake_path(GET file FILENAME filename)
      configure_file("${file}" "${directory}/${name}/${filename}" COPYONLY)
      list(APPEND filenames "${filename}")
    endforeach()
    list(JOIN filenames ", " filenames)
    list(APPEND index "${name}: ${spdx} (${filenames})\n")
  endforeach()
  string(JOIN "" index ${index})
  file(WRITE "${directory}/THIRD_PARTY.txt" "${index}")

  arts_bundled_license_expression(expression)
  set(${expression_variable} "${expression}" PARENT_SCOPE)
endfunction()

function(arts_bundled_license_expression variable)
  get_property(names GLOBAL PROPERTY ARTS_BUNDLED)
  set(expressions)
  foreach(name IN LISTS names)
    get_property(spdx GLOBAL PROPERTY ARTS_BUNDLED_${name}_SPDX)
    if(NOT spdx IN_LIST expressions)
      list(APPEND expressions "${spdx}")
    endif()
  endforeach()
  list(SORT expressions)
  set(expression "(${ARTS_LICENSE_EXPRESSION})")
  foreach(spdx IN LISTS expressions)
    # Parenthesize compound expressions so that AND binds as intended
    if(spdx MATCHES " (AND|OR) ")
      set(spdx "(${spdx})")
    endif()
    string(APPEND expression " AND ${spdx}")
  endforeach()
  set(${variable} "${expression}" PARENT_SCOPE)
endfunction()

function(arts_write_bundled_license_source output)
  get_property(names GLOBAL PROPERTY ARTS_BUNDLED)
  arts_bundled_license_expression(expression)
  set(data "")
  set(components "")
  set(i 0)
  set(c 0)
  foreach(name IN LISTS names)
    get_property(spdx GLOBAL PROPERTY ARTS_BUNDLED_${name}_SPDX)
    get_property(files GLOBAL PROPERTY ARTS_BUNDLED_${name}_FILES)
    set(entries "")
    foreach(file IN LISTS files)
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${file}")
      cmake_path(GET file FILENAME filename)
      # Bytes rather than a string literal: no escaping, and no compiler
      # limit on the length of a literal (MSVC's is about 16 kB).
      file(READ "${file}" hex HEX)
      string(LENGTH "${hex}" length)
      math(EXPR length "${length} / 2")
      string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${hex}")
      string(APPEND data "const unsigned char text_${i}[] = {${bytes}0x00};\n")
      string(APPEND entries "    {\"${filename}\", {reinterpret_cast<const char*>(text_${i}), ${length}}},\n")
      math(EXPR i "${i} + 1")
    endforeach()
    string(APPEND data "const license_file files_${c}[] = {\n${entries}};\n")
    string(APPEND components "    {\"${name}\", \"${spdx}\", files_${c}},\n")
    math(EXPR c "${c} + 1")
  endforeach()
  set(content "// Generated by arts_write_bundled_license_source (cmake/modules/ArtsBundledLicenses.cmake)
#include \"bundled_licenses.h\"

namespace arts {
namespace {
${data}
const bundled_component components[] = {
${components}};
}  // namespace

std::span<const bundled_component> bundled_components() { return components; }

std::string_view license_expression() { return \"${expression}\"; }
}  // namespace arts
")
  file(WRITE "${output}.tmp" "${content}")
  configure_file("${output}.tmp" "${output}" COPYONLY)
endfunction()
