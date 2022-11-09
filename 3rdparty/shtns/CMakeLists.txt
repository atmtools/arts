add_subdirectory(SHT)

file(GLOB_RECURSE shtns_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
foreach(f IN LISTS shtns_files)
  configure_file(${f} ${f} COPYONLY)
endforeach()

add_library(
  shtns
  sht_init.c
  sht_kernels_a.c
  sht_kernels_s.c
  sht_odd_nlat.c
  sht_fly.c
  sht_omp.c
  )
target_compile_options(
  shtns
  PRIVATE
  -O2 -march=native -mtune=skylake -O3 -ffast-math
  -fopenmp-simd -fno-tree-loop-distribute-patterns -mno-vzeroupper
  -std=gnu99 -D_GNU_SOURCE -DSHTNS_GIT=\"ARTS\"
)
target_include_directories(shtns PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_dependencies(shtns shtns_generated_files)

