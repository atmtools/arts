#include "atm.h"
#include "debug.h"
#include "gridded_fields.h"
#include "igrf13.h"
#include "matpack_data.h"
#include "quantum_numbers.h"
#include "species_tags.h"
#include "xml_io.h"
#include "xml_io_arts_types.h"
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>

void atm_fieldTopOfAtmosphere(AtmField &atm_field,
                              const Numeric &top_of_atmosphere) {
  atm_field.top_of_atmosphere = top_of_atmosphere;
}

void atm_fieldInit(AtmField &atm_field, const Numeric &top_of_atmosphere) {
  atm_field = AtmField{};
  atm_fieldTopOfAtmosphere(atm_field, top_of_atmosphere);
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
 * @return true If everything went well
 * @return false If everything went wrong
 */
template <typename T>
bool try_read(AtmField &atm_field, const Atm::KeyVal &key_val,
              const String& filename) {
  try {
    T v;
    xml_read_from_file(filename, v);
    std::visit([&](auto &key) { atm_field[key] = v; }, key_val);

    // Return success state
    return true;
  } catch (...) {
    //! CONTROL FLOW --- ARTS CANNOT READ VARIADICALLY

    // We must clean the atm_field so it can work again in the future
    std::visit([&](auto &key) { atm_field.erase(key); }, key_val);

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
 * @return true If everything went well
 * @return false If everything went wrong
 */
template <typename... T>
bool try_reading(AtmField &atm_field, const Atm::KeyVal &key_val,
                 const String& filename) {
  return (try_read<T>(atm_field, key_val, filename) or ...);
}

void atm_fieldAddCustomDataFileImpl(AtmField &atm_field,
                                    const Atm::KeyVal &key_val,
                                    const String &filename,
                                    const Atm::Extrapolation& extrapolation) {
  const bool ok = try_reading<GriddedField3, Numeric>(
      atm_field, key_val, filename);

  ARTS_USER_ERROR_IF(not ok, "The file ", std::quoted(filename),
                     " cannot be understood as atmospheric field data.\n"
                     "Please make sure that the file exists, that it is "
                     "possible to read the file, and\n"
                     "that its type is one that can be handled by "
                     "ARTS atmospheric fields")
  
  atm_field[key_val].alt_low = extrapolation;
  atm_field[key_val].lat_low = extrapolation;
  atm_field[key_val].lon_low = extrapolation;
  atm_field[key_val].alt_upp = extrapolation;
  atm_field[key_val].lat_upp = extrapolation;
  atm_field[key_val].lon_upp = extrapolation;
}
} // namespace detail

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const String &atmospheric_key,
                                const String &filename,
                                const String &extrapolation_type) {
  detail::atm_fieldAddCustomDataFileImpl(
      atm_field, Atm::toKeyOrThrow(atmospheric_key), filename,
      Atm::toExtrapolationOrThrow(extrapolation_type));
}

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const QuantumIdentifier &nlte_key,
                                const String &filename,
                                const String &extrapolation_type) {
  detail::atm_fieldAddCustomDataFileImpl(
      atm_field, Atm::KeyVal{nlte_key}, filename,
      Atm::toExtrapolationOrThrow(extrapolation_type));
}

void atm_fieldAddCustomDataFile(AtmField &atm_field,
                                const ArrayOfSpeciesTag &spec_key,
                                const String &filename,
                                const String &extrapolation_type) {
  detail::atm_fieldAddCustomDataFileImpl(
      atm_field, Atm::KeyVal{spec_key}, filename,
      Atm::toExtrapolationOrThrow(extrapolation_type));
}

void atm_fieldAddField(AtmField &atm_field, const String &filename,
                       const Index &set_top_of_atmosphere) {
  AtmField atm_field_other;
  xml_read_from_file(filename, atm_field_other);
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
                   const Index& read_nlte) {
  using enum Atm::Key;

  // Fix filename
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Reset and initialize
  atm_fieldInit(atm_field, top_of_atmosphere);

  if (read_tp) {
    for (auto &key : {t, p}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0);
    }
  }

  if (read_mag) {
    for (auto &key : {mag_u, mag_v, mag_w}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0);
    }
  }

  if (read_wind) {
    for (auto &key : {wind_u, wind_v, wind_w}) {
      const String file_name{var_string(tmp_basename, key, ".xml")};
      atm_fieldAddField(atm_field, file_name, 0);
    }
  }

  if (read_specs) {
    for (auto &spec : abs_species) {
      const String file_name{
          var_string(tmp_basename, toShortName(spec.Species()), ".xml")};
      atm_fieldAddField(atm_field, file_name, 0);
    }
  }

  if (read_nlte) {
    const String file_name{var_string(tmp_basename, "nlte.xml")};
    atm_fieldAddField(atm_field, file_name, 0);
  }
}

void atm_fieldSave(const AtmField &atm_field, const String &basename,
                   const String &filetype, const Index &no_clobber) {
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
      xml_write_to_file(filename, out, ftype, no_clobber);
    }
  }

  if (nlte.nnlte()) {
    const String filename{var_string(tmp_basename, "nlte.xml")};
    xml_write_to_file(filename, nlte, ftype, no_clobber);
  }
}

void atm_fieldAddGriddedData(AtmField &atm_field, const String &key,
                             const GriddedField3 &data, 
                                const String &extrapolation_type) {
  auto& fld = atm_field[Atm::toKeyOrThrow(key)] = data;

  const auto extrapolation = Atm::toExtrapolationOrThrow(extrapolation_type);
  
  fld.alt_low = extrapolation;
  fld.lat_low = extrapolation;
  fld.lon_low = extrapolation;
  fld.alt_upp = extrapolation;
  fld.lat_upp = extrapolation;
  fld.lon_upp = extrapolation;
}

void atm_fieldAddGriddedData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const GriddedField3 &data,
                                const String& extrapolation_type) {
  auto& fld = atm_field[key] = data;

  const auto extrapolation = Atm::toExtrapolationOrThrow(extrapolation_type);
  
  fld.alt_low = extrapolation;
  fld.lat_low = extrapolation;
  fld.lon_low = extrapolation;
  fld.alt_upp = extrapolation;
  fld.lat_upp = extrapolation;
  fld.lon_upp = extrapolation;
}

void atm_fieldAddGriddedData(AtmField &atm_field, const QuantumIdentifier &key,
                             const GriddedField3 &data,
                                const String& extrapolation_type) {
  auto& fld = atm_field[key] = data;

  const auto extrapolation = Atm::toExtrapolationOrThrow(extrapolation_type);
  
  fld.alt_low = extrapolation;
  fld.lat_low = extrapolation;
  fld.lon_low = extrapolation;
  fld.alt_upp = extrapolation;
  fld.lat_upp = extrapolation;
  fld.lon_upp = extrapolation;
}

void atm_fieldAddNumericData(AtmField &atm_field, const String &key,
                             const Numeric& data) {
  atm_field[Atm::toKeyOrThrow(key)] = data;
}

void atm_fieldAddNumericData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const Numeric& data) {
  atm_field[key] = data;
}

void atm_fieldAddNumericData(AtmField &atm_field, const QuantumIdentifier &key,
                             const Numeric& data) {
  atm_field[key] = data;
}

void atm_fieldIGRF(AtmField &atm_field, const Time &time, const Index &parsafe) {
  using namespace IGRF;

  //! We need explicit planet-size as IGRF requires the radius
  //! This is the WGS84 version of that, with radius and flattening
  static const Vector ell{6378137., 0.081819190842621};

  //! This struct deals with the computations, it's internally saving a mutable
  //! state
  struct res {
    mutable Numeric u{}, v{}, w{};
    mutable Numeric alt{}, lat{}, lon{};
    const Time t;

    res(const Time &x) : t(x) {}

    void comp(Numeric al, Numeric la, Numeric lo) const {
      if (alt != al or lat != la or lo != lon) {
        alt = al;
        lat = la;
        lon = lo;

        auto f = compute(Vector{alt}.reshape(1, 1, 1), {lat}, {lon}, t, ell);
        u = f.u(0, 0, 0);
        v = f.v(0, 0, 0);
        w = f.w(0, 0, 0);
      }
    }

    Numeric get_u(Numeric z, Numeric la, Numeric lo) const {
      comp(z, la, lo);
      return u;
    }

    Numeric get_v(Numeric z, Numeric la, Numeric lo) const {
      comp(z, la, lo);
      return v;
    }

    Numeric get_w(Numeric z, Numeric la, Numeric lo) const {
      comp(z, la, lo);
      return w;
    }
  };

  /** Different methods for parallel and non-parallel safe versions
   * In the paralell safe version each field component holds a copy
   * of the res-struct.  In the unsafe version, the res-struct is
   * shared between the fields.
   *
   * The parallel version is thus safe because it is doing all its
   * computations locally but the non-parallel version is 3X faster
   * because it saves the state of the query for, e.g., u to be
   * reused by v and w (with no regards for order) */
  if (parsafe == 0) {
    const std::shared_ptr igrf_ptr = std::make_shared<res>(time);
    atm_field[Atm::Key::mag_u] = [ptr = igrf_ptr](Numeric h, Numeric lat,
                                                  Numeric lon) {
      return ptr->get_u(h, lat, lon);
    };
    atm_field[Atm::Key::mag_v] = [ptr = igrf_ptr](Numeric h, Numeric lat,
                                                  Numeric lon) {
      return ptr->get_v(h, lat, lon);
    };
    atm_field[Atm::Key::mag_w] = [ptr = igrf_ptr](Numeric h, Numeric lat,
                                                  Numeric lon) {
      return ptr->get_w(h, lat, lon);
    };
  } else {
    atm_field[Atm::Key::mag_u] = [cpy = res(time)](Numeric h, Numeric lat,
                                                   Numeric lon) {
      return res{cpy}.get_u(h, lat, lon);
    };
    atm_field[Atm::Key::mag_v] = [cpy = res(time)](Numeric h, Numeric lat,
                                                   Numeric lon) {
      return res{cpy}.get_v(h, lat, lon);
    };
    atm_field[Atm::Key::mag_w] = [cpy = res(time)](Numeric h, Numeric lat,
                                                   Numeric lon) {
      return res{cpy}.get_w(h, lat, lon);
    };
  }
}
