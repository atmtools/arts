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
#include <variant>

void atm_fieldInit(AtmField &atm_field, const Numeric &top_of_atmosphere,
                   const Verbosity &) {
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


void xml_read_from_stream_helper(istream &is_xml, Atm::KeyVal &key_val,
                                 Atm::Data &data, bifstream *pbifs,
                                 const Verbosity &verbosity);

void atm_fieldAddDataFile(AtmField &atm_field, const String &filename,
                          const Verbosity &verbosity) {
  String fn = filename;
  find_xml_file(fn, verbosity);

  CREATE_OUT2;

  out2 << "  Reading " + filename + '\n';

  // Open input stream:
  std::unique_ptr<istream> ifs;
  if (filename.nelem() > 2 &&
      filename.substr(filename.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::make_unique<igzstream>();
    xml_open_input_file(*static_cast<igzstream *>(ifs.get()), filename,
                        verbosity);
  }
#else
  {
    throw runtime_error("This arts version was compiled without zlib support.\n"
                        "Thus zipped xml files cannot be read.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ifs = std::make_unique<ifstream>();
    xml_open_input_file(*static_cast<ifstream *>(ifs.get()), filename,
                        verbosity);
  }

  // No need to check for error, because xml_open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the
  // filename.
  try {
    FileType ftype;
    NumericType ntype;
    EndianType etype;

    Atm::KeyVal key;
    Atm::Data data;

    xml_read_header_from_stream(*ifs, ftype, ntype, etype, verbosity);
    if (ftype == FILE_TYPE_ASCII) {
      xml_read_from_stream_helper(*ifs, key, data, nullptr, verbosity);
    } else {
      String bfilename = filename + ".bin";
      bifstream bifs(bfilename.c_str());
      xml_read_from_stream_helper(*ifs, key, data, &bifs, verbosity);
    }
    xml_read_footer_from_stream(*ifs, verbosity);

    atm_field[key] = data;
  } catch (const std::runtime_error &e) {
    ostringstream os;
    os << "Error reading file: " << filename << '\n' << e.what();
    throw runtime_error(os.str());
  }
}


void atm_fieldReadAtm(AtmField &atm_field, const String &basename,
                   const ArrayOfArrayOfSpeciesTag &abs_species,
                   const Verbosity &verbosity) {
  CREATE_OUT3;

  // Fix filename
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Read the temperature field:
  String file_name = tmp_basename + "t.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Temperature field read from file: " << std::quoted(file_name) << "\n";

  // Read geometrical altitude field:
  file_name = tmp_basename + "z.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Altitude field read from file: " << std::quoted(file_name) << "\n";

  // We need to read one profile for each tag group.
  for (Index i = 0; i < abs_species.nelem(); i++) {
    // Determine the name.
    file_name = tmp_basename +
                String(Species::toShortName(abs_species[i].Species())) + ".xml";

    // Read the VMR:
    atm_fieldAddDataFile(atm_field, file_name, verbosity);
    // state the source of profile.
    out3 << "  " << Species::toShortName(abs_species[i].Species())
         << " profile read from file: " << std::quoted(file_name) << "\n";
  }
}


void atm_fieldReadWind(AtmField &atm_field, const String &basename,
                   const Verbosity &verbosity) {
  CREATE_OUT3;

  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Read wind field u component:
  String file_name = tmp_basename + "wind_u.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Wind u field read from file: " << std::quoted(file_name) << "\n";

  // Read wind field u component:
  file_name = tmp_basename + "wind_v.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Wind v field read from file: " << std::quoted(file_name) << "\n";

  // Read wind field u component:
  file_name = tmp_basename + "wind_w.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Wind w field read from file: " << std::quoted(file_name) << "\n";
}


void atm_fieldReadMag(AtmField &atm_field, const String &basename,
                   const Verbosity &verbosity) {
  CREATE_OUT3;

  String tmp_basename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Read mag field u component:
  String file_name = tmp_basename + "mag_u.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Mag u field read from file: " << std::quoted(file_name) << "\n";

  // Read mag field u component:
  file_name = tmp_basename + "mag_v.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Mag v field read from file: " << std::quoted(file_name) << "\n";

  // Read mag field u component:
  file_name = tmp_basename + "mag_w.xml";
  atm_fieldAddDataFile(atm_field, file_name, verbosity);
  out3 << "Mag w field read from file: " << std::quoted(file_name) << "\n";
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


void atm_fieldAddRegularData(AtmField &atm_field, const String &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
  atm_field[Atm::toKeyOrThrow(key)] = data;
}

void atm_fieldAddRegularData(AtmField &atm_field, const ArrayOfSpeciesTag &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
  atm_field[key] = data;
}

void atm_fieldAddRegularData(AtmField &atm_field, const QuantumIdentifier &key,
                             const Tensor3 &data, const Verbosity &) {
  ARTS_USER_ERROR_IF(
      not atm_field.regularized,
      "Cannot add regular data to non-regularized atmospheric field")
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
