# The C++ port in src/core/tmatrix is a derivative of Mishchenko's code, whose
# licence permits use in not-for-profit research and has no SPDX identifier.
# It is registered here, before src/core/licenses collects the registrations.
if(ENABLE_TMATRIX)
  if(ENABLE_ARTS_LGPL)
    message(FATAL_ERROR "T-matrix is unavailable in ENABLE_ARTS_LGPL builds")
  endif()
  arts_add_bundled_license(tmatrix SPDX LicenseRef-Mishchenko-T-matrix FILES License.txt)
endif()
