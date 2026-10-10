add_subdirectory(wigxjpf)
arts_add_bundled_license(wigxjpf SPDX LGPL-3.0-or-later FILES NOTICE wigxjpf/COPYING.LESSER wigxjpf/COPYING)

if(FASTWIGNER)
  add_subdirectory(fastwigxj)
  arts_add_bundled_license(fastwigxj SPDX LGPL-3.0-or-later FILES NOTICE fastwigxj/COPYING.LESSER fastwigxj/COPYING)
endif()

if(FASTWIGNER)
  add_library(wigner ALIAS fastwigxj)
else()
  add_library(wigner ALIAS wigxjpf)
endif()

if (FASTWIGNER)
  message (STATUS "Fast Wigner calculations enabled")
else()
  message (STATUS "Fast Wigner calculations disabled (use -DFASTWIGNER=1 to enable if your compiler supports it)")
endif()
