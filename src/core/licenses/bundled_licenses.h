#pragma once

#include <span>
#include <string_view>

/** The licences of the code from other projects compiled into this build.
 *
 * The contents are generated at configure time from the components
 * registered with arts_add_bundled_license (cmake/modules/ArtsBundledLicenses.cmake),
 * so every binary of ARTS carries the licence notices of the code it
 * contains.  See doc/arts/dev.licenses.rst.
 */
namespace arts {
//! One licence or notice file of a bundled component, as compiled in.
struct license_file {
  std::string_view name;  //!< The file name, e.g. "LICENSE"
  std::string_view text;  //!< The complete file
};

//! Code from another project that this build compiles in.
struct bundled_component {
  std::string_view              name;   //!< E.g. "Faddeeva"
  std::string_view              spdx;   //!< SPDX licence expression, or LicenseRef-...
  std::span<const license_file> files;  //!< Its licence and notice files
};

//! The bundled components of this build.
[[nodiscard]] std::span<const bundled_component> bundled_components();

//! The licence expression of this build: ARTS's licence AND those of the bundled components.
[[nodiscard]] std::string_view license_expression();
}  // namespace arts
