#include "atm.h"
#include "debug.h"
#include "gridded_fields.h"
#include "messages.h"
#include "xml_io.h"
#include "xml_io_arts_types.h"
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <variant>

void atm_fieldInit(AtmField &atm_field, const Numeric& top_of_atmosphere, const Verbosity &) {
  atm_field = AtmField{};
  atm_field.top_of_atmosphere = top_of_atmosphere;
}

void atm_fieldRegularize(AtmField &atm_field, const Vector &z,
                         const Vector &lat, const Vector &lon,
                         const Verbosity &) {
  atm_field.regularize(z, lat, lon);
}

namespace detail {
/** Tries to read a file as if it were some type T
 *
 * Assigns the value of the read to the atm_field at key_val
 *
 * @tparam T The type
 * @param atm_field As WSV
 * @param key_val A key value
 * @param filename A filename
 * @param verbosity As WSV
 * @return true If everything went well
 * @return false If everything went wrong
 */
template <typename T>
bool try_read(AtmField &atm_field, const Atm::KeyVal &key_val,
              const String &filename, const Verbosity &verbosity) {
  try {
    T v;
    xml_read_from_file(filename, v, verbosity);
    std::visit([&](auto &key) { atm_field[key] = v; }, key_val);
    atm_field.throwing_check();

    // Return success state
    return true;
  } catch (...) {
    //! CONTROL FLOW --- ARTS CANNOT READ VARIADICALLY

    // We must clean the atm_field so it can work again in the  future
    std::visit([&](auto &key) { atm_field.erase_key(key); }, key_val);

    // Return failure state
    return false;
  }
}

//! Calls try_read for the listed types
template <typename... T>
bool try_reading(AtmField &atm_field, const Atm::KeyVal &key_val,
                 const String &filename, const Verbosity &verbosity) {
  // Short-circuited read
  return (detail::try_read<T>(atm_field, key_val, filename, verbosity) or ...);
}
} // namespace detail

void atm_fieldSetFromCustomFile(AtmField &atm_field, const Atm::KeyVal &key_val,
                                const String &filename,
                                const Verbosity &verbosity) {
  const bool ok = detail::try_reading<GriddedField3, Tensor3, Numeric>(
      atm_field, key_val, filename, verbosity);

  ARTS_USER_ERROR_IF(not ok, "The file ", std::quoted(filename),
                     " cannot be understood as atmospheric field data.\n"
                     "Please make sure that the file exists, that it is "
                     "possible to read the file, and\n"
                     "that its type is one that can be handled by "
                     "ARTS atmospheric fields")
}

void atm_fieldSetFromCustomFile(AtmField &atm_field,
                                const String &atmospheric_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  atm_fieldSetFromCustomFile(atm_field, Atm::toKeyOrThrow(atmospheric_key),
                             filename, verbosity);
}

void atm_fieldSetFromCustomFile(AtmField &atm_field,
                                const QuantumIdentifier &nlte_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  atm_fieldSetFromCustomFile(atm_field, Atm::KeyVal{nlte_key}, filename,
                             verbosity);
}

void atm_fieldSetFromCustomFile(AtmField &atm_field,
                                const ArrayOfSpeciesTag &spec_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  atm_fieldSetFromCustomFile(atm_field, Atm::KeyVal{spec_key}, filename,
                             verbosity);
}
