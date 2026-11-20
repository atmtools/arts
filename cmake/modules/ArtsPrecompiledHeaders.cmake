# Precompiled headers, see the ENABLE_PCH option
#
# matpack.h brings the standard library, the XML IO and the format helpers that
# nearly every ARTS translation unit needs.  The matpack target precompiles it,
# and every other target that links matpack reuses a precompiled copy instead
# of building its own.
#
# A precompiled header can only be used with the compiler options it was made
# with.  Libraries are compiled as position-independent code (-fPIC) and
# executables as position-independent executables (-fPIE), so a library's copy
# is of no use to an executable: Clang refuses it and GCC cannot use it.
# Executables therefore reuse a second copy, made by the empty executable
# matpack_pch_executable.
#
# A target is left alone if it
#   - has a precompiled header of its own,
#   - has other compile options than the target that makes its copy, or
#   - sets DISABLE_PRECOMPILE_HEADERS, because forcing matpack.h on its sources
#     costs more than it saves.
#
# GCC reports a reused header that it cannot use as an error (-Werror=invalid-pch,
# see CMakeLists.txt).  That also catches a compile definition of a target that
# changes what matpack.h means, which Clang and MSVC do not check.

#! Whether target links matpack, directly or through the targets it links
function(arts_links_matpack target result)
  set(${result} FALSE PARENT_SCOPE)
  set(todo ${target})
  set(seen)
  while(todo)
    list(POP_FRONT todo t)
    if(t IN_LIST seen OR NOT TARGET ${t})
      continue()
    endif()
    list(APPEND seen ${t})
    get_target_property(type ${t} TYPE)
    if(type STREQUAL "INTERFACE_LIBRARY")
      get_target_property(deps ${t} INTERFACE_LINK_LIBRARIES)
    else()
      get_target_property(deps ${t} LINK_LIBRARIES)
    endif()
    if(deps)
      if("matpack" IN_LIST deps)
        set(${result} TRUE PARENT_SCOPE)
        return()
      endif()
      list(APPEND todo ${deps})
    endif()
  endwhile()
endfunction()

#! All targets defined in dir and its subdirectories
function(arts_targets_in dir out)
  get_property(targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
  get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
  foreach(sub IN LISTS subdirs)
    arts_targets_in(${sub} more)
    list(APPEND targets ${more})
  endforeach()
  set(${out} ${targets} PARENT_SCOPE)
endfunction()

#! The compile options of target, empty if it has none
function(arts_compile_options target out)
  get_target_property(options ${target} COMPILE_OPTIONS)
  if(NOT options)
    set(options "")
  endif()
  set(${out} "${options}" PARENT_SCOPE)
endfunction()

#! Let every target below dir that links matpack reuse a precompiled matpack.h
function(arts_reuse_matpack_pch dir)
  set(main ${CMAKE_CURRENT_BINARY_DIR}/matpack_pch_executable.cc)
  file(CONFIGURE OUTPUT ${main} CONTENT "int main() { return 0; }\n")
  add_executable(matpack_pch_executable EXCLUDE_FROM_ALL ${main})
  target_link_libraries(matpack_pch_executable PRIVATE matpack)
  target_precompile_headers(matpack_pch_executable PRIVATE <matpack.h>)

  arts_compile_options(matpack library_options)
  arts_compile_options(matpack_pch_executable executable_options)

  set(library_count 0)
  set(executable_count 0)
  arts_targets_in(${dir} targets)
  foreach(t IN LISTS targets)
    get_target_property(type ${t} TYPE)
    if(type STREQUAL "EXECUTABLE")
      set(kind executable)
      set(source matpack_pch_executable)
    elseif(type MATCHES "^(STATIC|OBJECT|SHARED|MODULE)_LIBRARY$")
      set(kind library)
      set(source matpack)
    else()
      continue()
    endif()

    get_target_property(own ${t} PRECOMPILE_HEADERS)
    get_target_property(reused ${t} PRECOMPILE_HEADERS_REUSE_FROM)
    get_target_property(disabled ${t} DISABLE_PRECOMPILE_HEADERS)
    arts_compile_options(${t} options)
    if(own OR reused OR disabled OR NOT "${options}" STREQUAL "${${kind}_options}")
      continue()
    endif()

    arts_links_matpack(${t} links)
    if(links)
      target_precompile_headers(${t} REUSE_FROM ${source})
      math(EXPR ${kind}_count "${${kind}_count} + 1")
    endif()
  endforeach()
  message(STATUS "Targets reusing the precompiled matpack.h: "
                 "${library_count} libraries, ${executable_count} executables")
endfunction()
