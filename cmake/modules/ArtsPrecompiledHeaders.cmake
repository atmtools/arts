# Precompiled headers, see the ENABLE_PCH option
#
# The matpack target precompiles matpack.h, which brings the standard library,
# the XML IO and the format helpers that nearly every ARTS translation unit
# needs.  Every other target that links matpack, and has no precompiled header
# of its own, reuses that one instead of building its own copy.

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

#! Let every target below dir that links matpack reuse its precompiled matpack.h
function(arts_reuse_matpack_pch dir)
  arts_targets_in(${dir} targets)
  set(users)
  foreach(t IN LISTS targets)
    get_target_property(type ${t} TYPE)
    get_target_property(own ${t} PRECOMPILE_HEADERS)
    if(t STREQUAL "matpack" OR own OR NOT type MATCHES "^(STATIC_LIBRARY|OBJECT_LIBRARY|SHARED_LIBRARY|MODULE_LIBRARY|EXECUTABLE)$")
      continue()
    endif()
    arts_links_matpack(${t} links)
    if(links)
      target_precompile_headers(${t} REUSE_FROM matpack)
      list(APPEND users ${t})
    endif()
  endforeach()
  list(LENGTH users n)
  message(STATUS "Targets reusing the precompiled matpack.h: ${n}")
endfunction()
