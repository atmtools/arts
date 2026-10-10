# Header-only components, compiled into every build
arts_add_bundled_license(invlib SPDX MIT FILES invlib/LICENSE.txt)
arts_add_bundled_license(mdspan SPDX "Apache-2.0 WITH LLVM-exception" FILES mdspan/LICENSE)

add_subdirectory(wigner)
add_subdirectory(Faddeeva)
add_subdirectory(cdisort)
add_subdirectory(shtns)

add_subdirectory(tmatrix)
add_subdirectory(polradtran)
