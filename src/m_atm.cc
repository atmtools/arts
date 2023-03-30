#include "atm.h"
#include "debug.h"
#include "gridded_fields.h"
#include "matpack_data.h"
#include "messages.h"
#include "quantum_numbers.h"
#include "species_tags.h"
#include "xml_io.h"
#include "xml_io_arts_types.h"
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string>
#include <variant>

void atm_fieldTopOfAtmosphere(AtmField &atm_field,
                              const Numeric &top_of_atmosphere,
                              const Verbosity &) {
  atm_field.top_of_atmosphere = top_of_atmosphere;
}

void atm_fieldInit(AtmField &atm_field, const Numeric &top_of_atmosphere,
                   const Verbosity &verbosity) {
  atm_field = AtmField{};
  atm_fieldTopOfAtmosphere(atm_field, top_of_atmosphere, verbosity);
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

/** Wraps try_read for multiple types
 *
 * The atm_field is updated, the order of short-circuiting is as
 * given by the types, i.e., the execution is (T1 or T2 or T3 or ...)
 *
 * @tparam T... The types
 * @param atm_field As WSV
 * @param key_val A key value
 * @param filename A filename
 * @param verbosity As WSV
 * @return true If everything went well
 * @return false If everything went wrong
 */
template <typename... T>
bool try_reading(AtmField &atm_field, const Atm::KeyVal &key_val,
                 const String &filename, const Verbosity &verbosity) {
  return (try_read<T>(atm_field, key_val, filename, verbosity) or ...);
}

void atm_fieldAddCustomDataFileImpl(AtmField &atm_field,
                                    const Atm::KeyVal &key_val,
                                    const String &filename,
                                    const Verbosity &verbosity) {
  const bool ok = try_reading<GriddedField3, Tensor3, Numeric>(
      atm_field, key_val, filename, verbosity);

  ARTS_USER_ERROR_IF(not ok, "The file ", std::quoted(filename),
                     " cannot be understood as atmospheric field data.\n"
                     "Please make sure that the file exists, that it is "
                     "possible to read the file, and\n"
                     "that its type is one that can be handled by "
                     "ARTS atmospheric fields")
}
} // namespace detail

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const String &atmospheric_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  detail::atm_fieldAddCustomDataFileImpl(
      atm_field, Atm::toKeyOrThrow(atmospheric_key), filename, verbosity);
}

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const QuantumIdentifier &nlte_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  detail::atm_fieldAddCustomDataFileImpl(atm_field, Atm::KeyVal{nlte_key},
                                         filename, verbosity);
}

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const ArrayOfSpeciesTag &spec_key,
                                const String &filename,
                                const Verbosity &verbosity) {
  detail::atm_fieldAddCustomDataFileImpl(atm_field, Atm::KeyVal{spec_key},
                                         filename, verbosity);
}

void atm_fieldAddField(AtmField &atm_field, const String &filename,
                       const Index &set_top_of_atmosphere,
                       const Verbosity &verbosity) {
  AtmField atm_field_other;
  xml_read_from_file(filename, atm_field_other, verbosity);
  for (auto &key : atm_field_other.keys()) {
    atm_field[key] = atm_field_other[key];
  }

  if (set_top_of_atmosphere)
    atm_field.top_of_atmosphere = atm_field_other.top_of_atmosphere;
}

void atm_fieldRead(AtmField &atm_field,
                   const ArrayOfArrayOfSpeciesTag &abs_species,
                   const String &basename, const Numeric &top_of_atmosphere,
                   const Index &read_tp, const Index &read_mag,
                   const Index &read_wind, const Index &read_specs,
                   const Index &read_nlte, const Verbosity &verbosity) {
  using enum Atm::Key;

  CREATE_OUT3;

  // Fix filename
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Reset and initialize
  atm_fieldInit(atm_field, top_of_atmosphere, verbosity);

  if (read_tp) {
    for (auto &key : {t, p}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0, verbosity);
      out3 << key << " field read from file: " << std::quoted(file_name)
           << '\n';
    }
  }

  if (read_mag) {
    for (auto &key : {mag_u, mag_v, mag_w}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0, verbosity);
      out3 << key << " field read from file: " << std::quoted(file_name)
           << '\n';
    }
  }

  if (read_wind) {
    for (auto &key : {wind_u, wind_v, wind_w}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0, verbosity);
      out3 << key << " field read from file: " << std::quoted(file_name)
           << '\n';
    }
  }

  if (read_specs) {
    for (auto &spec : abs_species) {
      const String file_name{
          var_string(tmp_basename, toShortName(spec.Species()), ".xml")};
      atm_fieldAddField(atm_field, file_name, 0, verbosity);
      out3 << std::quoted(var_string(spec))
           << " profile read from file: " << std::quoted(file_name) << "\n";
    }
  }

  if (read_nlte) {
    const String file_name{var_string(tmp_basename, "nlte.xml")};
    atm_fieldAddField(atm_field, file_name, 0, verbosity);
    out3 << "Read all NLTE fields from: " << std::quoted(file_name) << '\n';
  }

  atm_field.throwing_check();
}

void atm_fieldSave(const AtmField &atm_field, const String &basename,
                   const String &filetype, const Index &no_clobber,
                   const Verbosity &verbosity) {
  const auto ftype = string2filetype(filetype);

  // Fix filename
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  const auto keys = atm_field.keys();

  //
  ArrayOfSpecies specs{};
  ArrayOfIndex nspecs{};

  //
  AtmField nlte;
  nlte.top_of_atmosphere = atm_field.top_of_atmosphere;

  for (auto &key : keys) {
    if (std::holds_alternative<QuantumIdentifier>(key)) {
      nlte[key] = atm_field[key];
    } else {
      String keyname;
      if (std::holds_alternative<ArrayOfSpeciesTag>(key)) {
        auto *spec_key = std::get_if<ArrayOfSpeciesTag>(&key);
        auto spec = spec_key->Species();

        keyname = toString(spec);
        if (auto ptr = std::find(specs.begin(), specs.end(), spec);
            ptr == specs.end()) {
          specs.emplace_back(spec);
          nspecs.emplace_back(1);
        } else {
          const auto pos = std::distance(specs.begin(), ptr);
          nspecs[pos]++;
          keyname += var_string(".", nspecs[pos]);
        }
      } else if (std::holds_alternative<Atm::Key>(key)) {
        keyname = toString(*std::get_if<Atm::Key>(&key));
      } else {
        ARTS_ASSERT(false, "Failed to account for key type")
      }

      AtmField out;
      out.top_of_atmosphere = atm_field.top_of_atmosphere;
      out[key] = atm_field[key];
      const String filename{var_string(tmp_basename, keyname, ".xml")};
      xml_write_to_file(filename, out, ftype, no_clobber, verbosity);
    }
  }

  if (nlte.nnlte()) {
    const String filename{var_string(tmp_basename, "nlte.xml")};
    xml_write_to_file(filename, nlte, ftype, no_clobber, verbosity);
  }
}

void atm_fieldAddRegularData(AtmField &atm_field, const String &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
  ARTS_USER_ERROR_IF(data.shape() not_eq atm_field.regularized_shape(),
                     "Size-mismatch: field ",
                     matpack::shape_help{atm_field.regularized_shape()},
                     " vs input ", matpack::shape_help{data.shape()})
  atm_field[Atm::toKeyOrThrow(key)] = data;
}

void atm_fieldAddRegularData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
  ARTS_USER_ERROR_IF(data.shape() not_eq atm_field.regularized_shape(),
                     "Size-mismatch: field ",
                     matpack::shape_help{atm_field.regularized_shape()},
                     " vs input ", matpack::shape_help{data.shape()})
  atm_field[key] = data;
}

void atm_fieldAddRegularData(AtmField &atm_field, const QuantumIdentifier &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
  ARTS_USER_ERROR_IF(data.shape() not_eq atm_field.regularized_shape(),
                     "Size-mismatch: field ",
                     matpack::shape_help{atm_field.regularized_shape()},
                     " vs input ", matpack::shape_help{data.shape()})
  atm_field[key] = data;
}

void atm_fieldAddGriddedData(AtmField &atm_field, const String &key,
                             const GriddedField3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add gridded data to regularized atmospheric field")
  atm_field[Atm::toKeyOrThrow(key)] = data;
}

void atm_fieldAddGriddedData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const GriddedField3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add gridded data to regularized atmospheric field")
  atm_field[key] = data;
}

void atm_fieldAddGriddedData(AtmField &atm_field, const QuantumIdentifier &key,
                             const GriddedField3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add gridded data to regularized atmospheric field")
  atm_field[key] = data;
}

void atm_fieldAddNumericData(AtmField &atm_field, const String &key,
                             const Numeric &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add numeric data to regularized atmospheric field")
  atm_field[Atm::toKeyOrThrow(key)] = data;
}

void atm_fieldAddNumericData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const Numeric &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add numeric data to regularized atmospheric field")
  atm_field[key] = data;
}

void atm_fieldAddNumericData(AtmField &atm_field, const QuantumIdentifier &key,
                             const Numeric &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(atm_field.regularized,
                     "Cannot add numeric data to regularized atmospheric field")
  atm_field[key] = data;
}
