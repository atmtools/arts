/*!
  \file   m_cia.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \author Stefan Buehler
  \date   2012-12-04

  \brief  Workspace methods for HITRAN CIA data.

*/

#include <workspace.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>

#include "arts_constants.h"
#include "atm.h"
#include "cia.h"
#include "debug.h"
#include "file.h"
#include "jacobian.h"
#include "physics_funcs.h"
#include "species.h"
#include "species_tags.h"
#include "xml_io.h"

inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixAddCIA(  // WS Output:
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    // WS Input:
    const SpeciesEnum& select_species,
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& f_grid,
    const AtmPoint& atm_point,
    const ArrayOfCIARecord& absorption_cia_data,
    // WS Generic Input:
    const Numeric& T_extrapolfac,
    const Index& ignore_errors) {
  ARTS_TIME_REPORT

  // Size of problem
  const Index nf = f_grid.size();
  const Index nq = jacobian_targets.target_count();

  // Possible things that can go wrong in this code (excluding line parameters)
  ARTS_USER_ERROR_IF(static_cast<Index>(propagation_matrix.size()) not_eq nf,
                     "*f_grid* must match *propagation_matrix*")
  ARTS_USER_ERROR_IF(
      propagation_matrix_jacobian.nrows() not_eq nq,
      "*propagation_matrix_jacobian* must match derived form of *jacobian_targets*")
  ARTS_USER_ERROR_IF(
      propagation_matrix_jacobian.ncols() not_eq nf,
      "*propagation_matrix_jacobian* must have frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(any_negative(f_grid.vec()),
                     "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(atm_point.temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(atm_point.pressure <= 0, "Non-positive pressure")

  // Jacobian overhead START
  const auto jac_freqs = jacobian_targets.find_all<Jacobian::AtmTarget>(
      AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w);
  const auto jac_temps = jacobian_targets.find<Jacobian::AtmTarget>(AtmKey::t);

  const bool do_wind_jac =
      std::ranges::any_of(jac_freqs, [](const auto& x) { return x.first; });
  const bool do_temp_jac = jac_temps.first;
  const Numeric dt       = field_perturbation(std::span{&jac_temps, 1});
  const Numeric df       = field_perturbation(std::span{jac_freqs});

  Vector dfreq;
  Vector dabs_t{atm_point.temperature + dt};

  if (do_wind_jac) {
    dfreq.resize(f_grid.size());
    for (Size iv = 0; iv < f_grid.size(); iv++) dfreq[iv] = f_grid[iv] + df;
  }

  Vector dxsec_temp_dF(do_wind_jac ? f_grid.size() : 0),
      dxsec_temp_dT(do_temp_jac ? f_grid.size() : 0);

  ARTS_USER_ERROR_IF(do_wind_jac and !std::isnormal(df),
                     "df must be >0 and not NaN or Inf: {}",
                     df)
  ARTS_USER_ERROR_IF(do_temp_jac and !std::isnormal(dt),
                     "dt must be >0 and not NaN or Inf: {}",
                     dt)
  // Jacobian overhead END

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  // Allocate a vector with dimension frequencies for constructing our
  // cross-sections before adding them (more efficient to allocate this here
  // outside of the loops)
  Vector xsec_temp(f_grid.size());

  // Loop over CIA data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (const auto& this_cia : absorption_cia_data) {
    if (select_species != SpeciesEnum::Bath and
        select_species != this_cia.Species(0))
      continue;

    // We have to multiply with the number density of the second CIA species.
    // We do not have to multiply with the first, since we still
    // want to return a (unary) absorption cross-section, not an
    // absorption coefficient.
    const Numeric nd_sec =
        number_density(atm_point.pressure, atm_point.temperature) *
        atm_point[this_cia.Species(1)];
    // Get the binary absorption cross sections from the CIA data:

    try {
      this_cia.Extract(xsec_temp,
                       f_grid,
                       atm_point.temperature,
                       T_extrapolfac,
                       ignore_errors);
      if (do_wind_jac) {
        this_cia.Extract(dxsec_temp_dF,
                         dfreq,
                         atm_point.temperature,
                         T_extrapolfac,
                         ignore_errors);
      }
      if (do_temp_jac) {
        this_cia.Extract(dxsec_temp_dT,
                         f_grid,
                         atm_point.temperature + dt,
                         T_extrapolfac,
                         ignore_errors);
      }
    } catch (const std::runtime_error& e) {
      ARTS_USER_ERROR("Problem with CIA species {}-{}:\n{}",
                      this_cia.MoleculeName(0),
                      this_cia.MoleculeName(1),
                      e.what())
    }

    const Numeric nd =
        number_density(atm_point.pressure, atm_point.temperature);
    const Numeric dnd_dt =
        dnumber_density_dt(atm_point.pressure, atm_point.temperature);
    const Numeric dnd_dt_sec =
        dnumber_density_dt(atm_point.pressure, atm_point.temperature) *
        atm_point[this_cia.Species(1)];
    for (Size iv = 0; iv < f_grid.size(); iv++) {
      propagation_matrix[iv].A() +=
          nd_sec * xsec_temp[iv] * nd * atm_point[this_cia.Species(0)];

      if (jac_temps.first) {
        const auto iq = jac_temps.second->target_pos;
        propagation_matrix_jacobian[iq, iv].A() +=
            ((nd_sec * (dxsec_temp_dT[iv] - xsec_temp[iv]) / dt +
              xsec_temp[iv] * dnd_dt_sec) *
                 nd +
             xsec_temp[iv] * nd_sec * dnd_dt) *
            atm_point[this_cia.Species(0)];
      }

      for (auto& j : jac_freqs) {
        if (j.first) {
          const auto iq = j.second->target_pos;
          propagation_matrix_jacobian[iq, iv].A() +=
              nd_sec * (dxsec_temp_dF[iv] - xsec_temp[iv]) / df * nd *
              atm_point[this_cia.Species(1)];
        }
      }

      if (const auto j =
              jacobian_targets.find<Jacobian::AtmTarget>(this_cia.Species(0));
          j.first) {
        const auto iq                            = j.second->target_pos;
        propagation_matrix_jacobian[iq, iv].A() += nd_sec * xsec_temp[iv] * nd;
      }

      if (const auto j =
              jacobian_targets.find<Jacobian::AtmTarget>(this_cia.Species(1));
          j.first) {
        const auto iq                            = j.second->target_pos;
        propagation_matrix_jacobian[iq, iv].A() += nd_sec * xsec_temp[iv] * nd;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void CIARecordReadFromFile(  // WS GOutput:
    CIARecord& cia_record,
    // WS Generic Input:
    const String& species_tag,
    const String& filename) {
  ARTS_TIME_REPORT

  SpeciesTag species(species_tag);

  ARTS_USER_ERROR_IF(species.Type() != SpeciesTagType::Cia,
                     "Invalid species tag {}"
                     ".\n"
                     "This is not recognized as a CIA type.\n",
                     species_tag)

  cia_record.SetSpecies(species.Spec(), species.cia_2nd_species);
  cia_record.ReadFromCIA(filename);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_cia_dataAddCIARecord(  // WS Output:
    ArrayOfCIARecord& absorption_cia_data,
    // WS GInput:
    const CIARecord& cia_record,
    const Index& clobber) {
  ARTS_TIME_REPORT

  Index cia_index = cia_get_index(
      absorption_cia_data, cia_record.Species(0), cia_record.Species(1));
  if (cia_index == -1)
    absorption_cia_data.push_back(cia_record);
  else if (clobber)
    absorption_cia_data[cia_index] = cia_record;
  else
    absorption_cia_data[cia_index].AppendDataset(cia_record);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_cia_dataReadFromCIA(  // WS Output:
    ArrayOfCIARecord& absorption_cia_data,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& catalogpath) {
  ARTS_TIME_REPORT

  ArrayOfString subfolders;
  subfolders.push_back("Main-Folder/");
  subfolders.push_back("Alternate-Folder/");

  absorption_cia_data.resize(0);

  // Loop species tag groups to find CIA tags.
  // Index sp loops through the tag groups, index iso through the tags within
  // each group. Despite the name, iso does not denote the isotope!
  for (Size sp = 0; sp < abs_species.size(); sp++) {
    for (Size iso = 0; iso < abs_species[sp].size(); iso++) {
      if (abs_species[sp][iso].Type() != SpeciesTagType::Cia) continue;

      ArrayOfString cia_names;

      Index cia_index = cia_get_index(absorption_cia_data,
                                      abs_species[sp][iso].Spec(),
                                      abs_species[sp][iso].cia_2nd_species);

      // If cia_index is not -1, we have already read this datafile earlier
      if (cia_index != -1) continue;

      cia_names.push_back(
          String(toString<1>(abs_species[sp][iso].Spec())) + "-" +
          String(toString<1>(abs_species[sp][iso].cia_2nd_species)));

      cia_names.push_back(
          String(toString<1>(abs_species[sp][iso].cia_2nd_species)) + "-" +
          String(toString<1>(abs_species[sp][iso].Spec())));

      ArrayOfString checked_dirs;

      bool found = false;
      for (Size fname = 0; !found && fname < cia_names.size(); fname++) {
        String cia_name = cia_names[fname];

        for (Size dir = 0; !found && dir < subfolders.size(); dir++) {
          ArrayOfString files;
          checked_dirs.push_back(catalogpath + "/" + subfolders[dir] +
                                 cia_name + "/");
          try {
            files = list_directory(*(checked_dirs.end() - 1));
          } catch (const std::runtime_error& e) {
            continue;
          }

          for (Index i = files.size() - 1; i >= 0; i--) {
            if (files[i].find(cia_name) != 0 ||
                files[i].rfind(".cia") != files[i].length() - 4) {
              files.erase(files.begin() + i);
            }
          }
          if (files.size()) {
            CIARecord ciar;

            found          = true;
            String catfile = *(checked_dirs.end() - 1) + files[0];

            ciar.SetSpecies(abs_species[sp][iso].Spec(),
                            abs_species[sp][iso].cia_2nd_species);
            ciar.ReadFromCIA(catfile);

            absorption_cia_data.push_back(ciar);
          }
        }
      }

      ARTS_USER_ERROR_IF(!found,
                         "Error: No data file found for CIA species {}"
                         "\n"
                         "Looked in directories: {}",
                         cia_names[0],
                         checked_dirs)
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_cia_dataReadFromXML(  // WS Output:
    ArrayOfCIARecord& absorption_cia_data,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& filename) {
  ARTS_TIME_REPORT

  xml_read_from_file(filename, absorption_cia_data);

  // Check that all CIA tags from abs_species are present in the
  // XML file

  std::vector<String> missing_tags;

  // Loop species tag groups to find CIA tags.
  // Index sp loops through the tag groups, index iso through the tags within
  // each group. Despite the name, iso does not denote the isotope!
  for (Size sp = 0; sp < abs_species.size(); sp++) {
    for (Size iso = 0; iso < abs_species[sp].size(); iso++) {
      if (abs_species[sp][iso].Type() != SpeciesTagType::Cia) continue;

      Index cia_index = cia_get_index(absorption_cia_data,
                                      abs_species[sp][iso].Spec(),
                                      abs_species[sp][iso].cia_2nd_species);

      // If cia_index is -1, this CIA tag was not present in the input file
      if (cia_index == -1) {
        missing_tags.push_back(
            String(toString<1>(abs_species[sp][iso].Spec())) + "-" +
            String(toString<1>(abs_species[sp][iso].cia_2nd_species)));
      }
    }
  }

  if (missing_tags.size()) {
    std::ostringstream os;
    bool first = true;

    os << "Error: The following CIA tag(s) are missing in input file: ";
    for (Size i = 0; i < missing_tags.size(); i++) {
      if (!first)
        os << ", ";
      else
        first = false;
      os << missing_tags[i];
    }
    ARTS_USER_ERROR("{}", os.str());
  }
}

void absorption_cia_dataReadSpeciesSplitCatalog(
    ArrayOfCIARecord& absorption_cia_data,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& basename,
    const Index& ignore_missing_) try {
  ARTS_TIME_REPORT

  const bool ignore_missing = static_cast<bool>(ignore_missing_);

  ArrayOfString names{};
  for (auto& spec : abs_species) {
    for (auto& tag : spec) {
      if (tag.type == SpeciesTagType::Cia) {
        names.emplace_back(std::format("{}-CIA-{}",
                                       toString<1>(tag.Spec()),
                                       toString<1>(tag.cia_2nd_species)));
      }
    }
  }

  names.erase(std::unique(names.begin(), names.end()), names.end());

  const std::filesystem::path inpath{basename.c_str()};
  const bool is_dir =
      basename.back() == '/' or std::filesystem::is_directory(inpath);

  for (auto& name : names) {
    auto fil{inpath};
    if (is_dir) {
      fil /= name + ".xml";
    } else {
      fil += "." + name + ".xml";
    }

    String filename = fil.string();
    if (not find_xml_file_existence(filename)) {
      if (ignore_missing) continue;
      ARTS_USER_ERROR("File {} not found", filename);
    }

    xml_read_from_file(fil.string(), absorption_cia_data.emplace_back());
  }
}
ARTS_METHOD_ERROR_CATCH
