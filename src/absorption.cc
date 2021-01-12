/* Copyright (C) 2000-2012
   Stefan Buehler  <sbuehler@ltu.se>
   Axel von Engeln <engeln@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  \file   absorption.cc

  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Stefan Buehler and Axel von Engeln
*/

#include "absorption.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <map>
#include "arts.h"
#include "auto_md.h"
#include "file.h"
#include "linescaling.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"

#include "global_data.h"
#include "linefunctions.h"

/** Mapping of species auxiliary type names to SpeciesAuxData::AuxType enum */
static const char* SpeciesAuxTypeNames[] = {"NONE",
                                            "ISORATIO",  //Built-in type
                                            "ISOQUANTUM",
                                            "PART_TFIELD",
                                            "PART_COEFF",  //Built-in type
                                            "PART_COEFF_VIBROT"};

void SpeciesAuxData::InitFromSpeciesData() {
  using global_data::species_data;

  mparams.resize(species_data.nelem());
  mparam_type.resize(species_data.nelem());

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    const Index niso = species_data[isp].Isotopologue().nelem();
    mparams[isp].resize(niso);
    mparam_type[isp].resize(niso);
    for (Index iso = 0; iso < niso; iso++) {
      mparams[isp][iso].resize(0);
      mparam_type[isp][iso] = SpeciesAuxData::AT_NONE;
    }
  }
}

void SpeciesAuxData::setParam(const Index species,
                              const Index isotopologue,
                              const AuxType auxtype,
                              const ArrayOfGriddedField1& auxdata) {
  mparam_type[species][isotopologue] = auxtype;
  mparams[species][isotopologue] = auxdata;
}

void SpeciesAuxData::setParam(const String& artstag,
                              const String& auxtype,
                              const ArrayOfGriddedField1& auxdata) {
  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if (!hinit) {
    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];
      for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
        SpecIsoMap indicies(i, j);
        String buf = sr.Name() + "-" + sr.Isotopologue()[j].Name();
        ArtsMap[buf] = indicies;
      }
    }
    hinit = true;
  }

  Index species;
  Index isotopologue;

  // ok, now for the cool index map:
  // is this arts identifier valid?
  const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artstag);
  if (i == ArtsMap.end()) {
    ostringstream os;
    os << "ARTS Tag: " << artstag << " is unknown.";
    throw runtime_error(os.str());
  }

  SpecIsoMap id = i->second;

  // Set mspecies:
  species = id.Speciesindex();

  // Set misotopologue:
  isotopologue = id.Isotopologueindex();

  Index this_auxtype = 0;

  while (this_auxtype < AT_FINAL_ENTRY &&
         auxtype != SpeciesAuxTypeNames[this_auxtype])
    this_auxtype++;

  if (this_auxtype != AT_FINAL_ENTRY) {
    setParam(species, isotopologue, (AuxType)this_auxtype, auxdata);
  } else {
    ostringstream os;
    os << "Unknown SpeciesAuxData type: " << auxtype;
    std::runtime_error(os.str());
  }
}

const ArrayOfGriddedField1& SpeciesAuxData::getParam(
    const Index species, const Index isotopologue) const {
  return mparams[species][isotopologue];
}

Numeric SpeciesAuxData::getIsotopologueRatio(const SpeciesTag& st) const {
  Numeric val = 1.0;
  if (st.Isotopologue() > -1)
    val = getParam(st.Species(), st.Isotopologue())[0].data[0];
  return val;
}

Numeric SpeciesAuxData::getIsotopologueRatio(const QuantumIdentifier& qid) const {
  return getParam(qid.Species(), qid.Isotopologue())[0].data[0];
}

String SpeciesAuxData::getTypeString(const Index species,
                                     const Index isotopologue) const {
  assert(mparam_type[species][isotopologue] < AT_FINAL_ENTRY);
  return SpeciesAuxTypeNames[mparam_type[species][isotopologue]];
}

bool SpeciesAuxData::ReadFromStream(String& artsid,
                                    istream& is,
                                    Index nparams,
                                    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if (!hinit) {
    out3 << "  ARTS index table:\n";

    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];

      for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
        SpecIsoMap indicies(i, j);
        String buf = sr.Name() + "-" + sr.Isotopologue()[j].Name();

        ArtsMap[buf] = indicies;

        // Print the generated data structures (for debugging):
        // The explicit conversion of Name to a c-String is
        // necessary, because setw does not work correctly for
        // stl Strings.
        const Index& i1 = ArtsMap[buf].Speciesindex();
        const Index& i2 = ArtsMap[buf].Isotopologueindex();

        out3 << "  Arts Identifier = " << buf << "   Species = " << setw(10)
             << setiosflags(ios::left) << species_data[i1].Name().c_str()
             << "iso = " << species_data[i1].Isotopologue()[i2].Name().c_str()
             << "\n";
      }
    }
    hinit = true;
  }

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return true;

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return true;

    // @ as first character marks catalogue entry
    char c;
    extract(c, line, 1);

    // check for empty line
    if (c == '@') {
      comment = false;
    }
  }

  // read the arts identifier String
  istringstream icecream(line);

  icecream >> artsid;

  if (artsid.length() != 0) {
    Index mspecies;
    Index misotopologue;

    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
    if (i == ArtsMap.end()) {
      ostringstream os;
      os << "ARTS Tag: " << artsid << " is unknown.";
      throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;

    // Set mspecies:
    mspecies = id.Speciesindex();

    // Set misotopologue:
    misotopologue = id.Isotopologueindex();

    ArrayOfGriddedField1 ratios;
    ratios.resize(1);
    // Extract accuracies:
    try {
      Numeric p = NAN;
      std::vector<Numeric> aux;
      for (Index ip = 0; ip < nparams; ip++) {
        icecream >> double_imanip() >> p;
        aux.push_back(p);
      }

      Vector grid;
      if (aux.size() > 1)
        nlinspace(grid, 1, (Numeric)aux.size(), aux.size());
      else
        grid = Vector(1, .1);

      ratios[0].set_grid(0, grid);
      ratios[0].data = aux;
      mparams[mspecies][misotopologue] = ratios;
    } catch (const runtime_error&) {
      throw runtime_error("Error reading SpeciesAuxData.");
    }
  }

  // That's it!
  return false;
}

void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& isoratios) {
  using global_data::species_data;

  // Check total number of species:
  if (species_data.nelem() != isoratios.nspecies()) {
    ostringstream os;
    os << "Number of species in SpeciesAuxData (" << isoratios.nspecies()
       << ") "
       << "does not fit builtin species data (" << species_data.nelem() << ").";
    throw runtime_error(os.str());
  }

  // For the selected species, we check all isotopes by looping over the
  // species data. (Trying to check only the isotopes actually used gets
  // quite complicated, actually, so we do the simple thing here.)

  // Loop over the absorption species:
  for (Index i = 0; i < abs_species.nelem(); i++) {
    // sp is the index of this species in the internal lookup table
    const Index sp = abs_species[i][0].Species();

    // Get handle on species data for this species:
    const SpeciesRecord& this_sd = species_data[sp];

    // Check number of isotopologues:
    if (this_sd.Isotopologue().nelem() != isoratios.nisotopologues(sp)) {
      ostringstream os;
      os << "Incorrect number of isotopologues in isotopologue data.\n"
         << "Species: " << this_sd.Name() << ".\n"
         << "Number of isotopes in SpeciesAuxData ("
         << isoratios.nisotopologues(sp) << ") "
         << "does not fit builtin species data ("
         << this_sd.Isotopologue().nelem() << ").";
      throw runtime_error(os.str());
    }

    for (Index iso = 0; iso < this_sd.Isotopologue().nelem(); ++iso) {
      // For "real" species (not representing continau) the isotopologue
      // ratio must not be NAN or below zero.
      if (!this_sd.Isotopologue()[iso].isContinuum()) {
        if (std::isnan(isoratios.getParam(sp, iso)[0].data[0]) ||
            isoratios.getParam(sp, iso)[0].data[0] < 0.) {
          ostringstream os;
          os << "Invalid isotopologue ratio.\n"
             << "Species: " << this_sd.Name() << "-"
             << this_sd.Isotopologue()[iso].Name() << "\n"
             << "Ratio:   " << isoratios.getParam(sp, iso)[0].data[0];
          throw runtime_error(os.str());
        }
      }
    }
  }
}

void checkPartitionFunctions(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& partfun) {
  using global_data::species_data;

  // Check total number of species:
  if (species_data.nelem() != partfun.nspecies()) {
    ostringstream os;
    os << "Number of species in SpeciesAuxData (" << partfun.nspecies() << ") "
       << "does not fit builtin species data (" << species_data.nelem() << ").";
    throw runtime_error(os.str());
  }

  // For the selected species, we check all isotopes by looping over the
  // species data. (Trying to check only the isotopes actually used gets
  // quite complicated, actually, so we do the simple thing here.)

  // Loop over the absorption species:
  for (Index i = 0; i < abs_species.nelem(); i++) {
    // sp is the index of this species in the internal lookup table
    const Index sp = abs_species[i][0].Species();

    // Get handle on species data for this species:
    const SpeciesRecord& this_sd = species_data[sp];

    // Check number of isotopologues:
    if (this_sd.Isotopologue().nelem() != partfun.nisotopologues(sp)) {
      ostringstream os;
      os << "Incorrect number of isotopologues in partition function data.\n"
         << "Species: " << this_sd.Name() << ".\n"
         << "Number of isotopes in SpeciesAuxData ("
         << partfun.nisotopologues(sp) << ") "
         << "does not fit builtin species data ("
         << this_sd.Isotopologue().nelem() << ").";
      throw runtime_error(os.str());
    }
  }
}

void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(
    SpeciesAuxData& sad) {
  using global_data::species_data;

  sad.InitFromSpeciesData();

  Vector grid(1, 1.);
  ArrayOfGriddedField1 ratios;
  ratios.resize(1);
  ratios[0].set_name("IsoRatios");
  ratios[0].set_grid_name(0, "Index");
  ratios[0].set_grid(0, grid);
  ratios[0].resize(1);

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem();
         iiso++) {
      ratios[0].data[0] = species_data[isp].Isotopologue()[iiso].Abundance();
      sad.setParam(isp, iiso, SpeciesAuxData::AT_ISOTOPOLOGUE_RATIO, ratios);
    }
  }
}

void fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(
    SpeciesAuxData& sad) {
  using global_data::species_data;

  sad.InitFromSpeciesData();

  ArrayOfGriddedField1 pfuncs;
  pfuncs.resize(2);
  pfuncs[0].set_name("PartitionFunction");
  pfuncs[0].set_grid_name(0, "Coeff");
  pfuncs[1].set_grid_name(0, "Temperature");

  ArrayOfString tgrid;
  tgrid.resize(2);
  tgrid[0] = "Tlower";
  tgrid[1] = "Tupper";

  for (Index isp = 0; isp < species_data.nelem(); isp++) {
    for (Index iiso = 0; iiso < species_data[isp].Isotopologue().nelem();
         iiso++) {
      Vector grid;
      const Vector& coeffs = species_data[isp].Isotopologue()[iiso].GetCoeff();

      assert(coeffs.nelem() >= 2);

      nlinspace(grid, 0, (Numeric)coeffs.nelem() - 1., coeffs.nelem());
      pfuncs[0].set_grid(0, grid);
      pfuncs[0].data = coeffs;

      const Vector& temp_range =
          species_data[isp].Isotopologue()[iiso].GetCoeffGrid();

      // Temperature data should either contain two Ts, lower and upper value of
      // the valid range or be empty
      assert(temp_range.nelem() == 0 || temp_range.nelem() == 2);

      if (temp_range.nelem() == 2) {
        pfuncs[1].set_grid(0, tgrid);
        pfuncs[1].data = temp_range;
      } else {
        pfuncs[1].resize(0);
        pfuncs[1].set_grid(0, ArrayOfString());
      }

      sad.setParam(
          isp, iiso, SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, pfuncs);
    }
  }
}

ostream& operator<<(ostream& os, const SpeciesRecord& sr) {
  for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
    os << sr.Name() << "-" << sr.Isotopologue()[j].Name() << "\n";
  }
  return os;
}

ostream& operator<<(ostream& os, const SpeciesAuxData& sad) {
  using global_data::species_data;
  for (Index sp = 0; sp < sad.nspecies(); sp++) {
    for (Index iso = 0; iso < sad.nisotopologues(sp); iso++) {
      os << species_name_from_species_index(sp) << "-"
         << global_data::species_data[sp].Isotopologue()[iso].Name();
      os << " " << sad.getTypeString(sp, iso) << std::endl;
      for (Index ip = 0; ip < sad.getParam(sp, iso).nelem(); ip++)
        os << "AuxData " << ip << " " << sad.getParam(sp, iso) << std::endl;
    }
  }

  return os;
}

/** A little helper function to convert energy from units of
    wavenumber (cm^-1) to Joule (J). 

    This is used when reading HITRAN or JPL catalogue files, which
    have the lower state energy in cm^-1.

    \return Energy in J.
    \param[in]  e Energy in cm^-1.

    \author Stefan Buehler
    \date   2001-06-26 */
Numeric wavenumber_to_joule(Numeric e) {
  // Planck constant [Js]
  extern const Numeric PLANCK_CONST;

  // Speed of light [m/s]
  extern const Numeric SPEED_OF_LIGHT;

  // Constant to convert lower state energy from cm^-1 to J
  const Numeric lower_energy_const = PLANCK_CONST * SPEED_OF_LIGHT * 1E2;

  return e * lower_energy_const;
}

//======================================================================
//             Functions related to species
//======================================================================

//! Return species index for given species name.
/*! 
  This is useful in connection with other functions that need a species
  index.

  \see find_first_species_tg.

  \param[in] name Species name.

  \return Species index, -1 means not found.

  \author Stefan Buehler
  \date   2003-01-13
*/
Index species_index_from_species_name(String name) {
  using global_data::SpeciesMap;

  // For the return value:
  Index mspecies;

  // Trim whitespace
  name.trim();

  //  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<String, Index>::const_iterator mi = SpeciesMap.find(name);
  if (mi != SpeciesMap.end()) {
    // Ok, we've found the species. Set mspecies.
    mspecies = mi->second;
  } else {
    // The species does not exist!
    mspecies = -1;
  }

  return mspecies;
}

//! Return species name for given species index.
/*!
 This is useful in connection with other functions that use a species
 index.
 
 Does an assertion that the index really corresponds to a species.
 
 \param[in] spec_ind Species index.
 
 \return Species name
 
 \author Stefan Buehler
 \date   2013-01-04
 */
String species_name_from_species_index(const Index spec_ind) {
  // Species lookup data:
  using global_data::species_data;

  // Assert that spec_ind is inside species data. (This is an assertion,
  // because species indices should never be user input, but set by the
  // program automatically, based on species names.)
  assert(spec_ind >= 0);
  assert(spec_ind < species_data.nelem());

  // A reference to the relevant record of the species data:
  const SpeciesRecord& spr = species_data[spec_ind];

  return spr.Name();
}

//!  set_abs_from_first_species.
/*!
 Returns vmr for the profile of the first tag group containing
 the given species.

 \author Oliver Lemke

 \param[out] vmr          Volume mixing ratio
 \param[in]  species_name Species Name
 \param[in]  abs_species  WS Input
 \param[in]  abs_vmrs     WS Input
 */
void set_vmr_from_first_species(Vector& vmr,
                                const String& species_name,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Matrix& abs_vmrs) {
  const Index index = find_first_species_tg(
      abs_species, species_index_from_species_name (species_name));

  vmr.resize(abs_vmrs.ncols());
  if (index < 0)
    vmr = -99;
  else
    vmr = abs_vmrs(index, Range(joker));
}

const SpeciesRecord& SpeciesDataOfBand(const AbsorptionLines& band)
{
  return global_data::species_data[band.Species()];
}

void xsec_species(Matrix& xsec,
                  Matrix& source,
                  Matrix& phase,
                  ArrayOfMatrix& dxsec_dx,
                  ArrayOfMatrix& dsource_dx,
                  ArrayOfMatrix& dphase_dx,
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const ArrayOfIndex& jacobian_propmat_positions,
                  const Vector& f_grid,
                  const Vector& abs_p,
                  const Vector& abs_t,
                  const EnergyLevelMap& abs_nlte,
                  const Matrix& abs_vmrs,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const AbsorptionLines& band,
                  const Numeric& isot_ratio,
                  const SpeciesAuxData::AuxType& partfun_type,
                  const ArrayOfGriddedField1& partfun_data) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = band.NumLines();  // number of lines in the catalog
  const Index nj =
      jacobian_propmat_positions.nelem();  // number of partial derivatives
  const Index nt = source.nrows();         // number of energy levels in NLTE

  // Type of problem
  const bool do_nonlte = nt;

  Linefunctions::InternalData scratch(nf, nj);
  Linefunctions::InternalData sum(nf, nj);
  
  // Test if the size of the problem is 0
  if (not np or not nf or not nl) return;
  
  // Constant for all lines
  const Numeric QT0 = single_partition_function(band.T0(), partfun_type, partfun_data);

  ArrayOfString fail_msg;
  bool do_abort = false;

#pragma omp parallel for if (!arts_omp_in_parallel() && np > 1) \
    firstprivate(scratch, sum)
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      // Constants for this level
      const Numeric& temperature = abs_t[ip];
      const Numeric& pressure = abs_p[ip];

      // Constants for this level
      const Numeric QT =
          single_partition_function(temperature, partfun_type, partfun_data);
      const Numeric dQTdT = dsingle_partition_function_dT(
          temperature,
          partfun_type,
          partfun_data);
      const Numeric DC =
          Linefunctions::DopplerConstant(temperature, band.SpeciesMass());
      const Numeric dDCdT = Linefunctions::dDopplerConstant_dT(temperature, DC);
      const Vector line_shape_vmr =
          band.BroadeningSpeciesVMR(abs_vmrs(joker, ip), abs_species);

      Linefunctions::set_cross_section_of_band(scratch,
                                               sum,
                                               f_grid,
                                               band,
                                               jacobian_quantities,
                                               jacobian_propmat_positions,
                                               line_shape_vmr,
                                               abs_nlte[ip],
                                               pressure,
                                               temperature,
                                               isot_ratio,
                                               0,
                                               DC,
                                               dDCdT,
                                               QT,
                                               dQTdT,
                                               QT0,
                                               false);

      // absorption cross-section
      MapToEigen(xsec).col(ip).noalias() += sum.F.real();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dxsec_dx[j]).col(ip).noalias() += sum.dF.col(j).real();

      // phase cross-section
      if (not phase.empty()) {
        MapToEigen(phase).col(ip).noalias() += sum.F.imag();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dphase_dx[j]).col(ip).noalias() += sum.dF.col(j).imag();
      }

      // source ratio cross-section
      if (do_nonlte) {
        MapToEigen(source).col(ip).noalias() += sum.N.real();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dsource_dx[j]).col(ip).noalias() += sum.dN.col(j).real();
      }
    } catch (const std::runtime_error& e) {
      ostringstream os;
      os << "Runtime-error in cross-section calculation at p_abs index " << ip
         << ": \n";
      os << e.what();
#pragma omp critical(xsec_species_cross_sections)
      {
        do_abort = true;
        fail_msg.push_back(os.str());
      }
    }
  }

  if (do_abort) {
    std::ostringstream os;
    os << "Error messages from failed cases:\n";
    for (const auto& msg : fail_msg) {
      os << msg << '\n';
    }
    throw std::runtime_error(os.str());
  }
}
