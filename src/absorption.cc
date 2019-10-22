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
#include "interpolation_poly.h"
#include "linescaling.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"

#include "global_data.h"
#include "linefunctions.h"
#include "partial_derivatives.h"

/** Mapping of species auxiliary type names to SpeciesAuxData::AuxType enum */
static const char* SpeciesAuxTypeNames[] = {"NONE",
                                            "ISORATIO",  //Built-in type
                                            "ISOQUANTUM",
                                            "PART_TFIELD",
                                            "PART_COEFF",  //Built-in type
                                            "PART_COEFF_VIBROT"};

// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromCoeff_OldUnused(
    Numeric temperature) const {
  Numeric result = 0.;
  Numeric exponent = 1.;

  Vector::const_iterator it;

  //      cout << "T: " << temperature << "\n";
  for (it = mqcoeff.begin(); it != mqcoeff.end(); ++it) {
    result += *it * exponent;
    exponent *= temperature;
    //      cout << "it: " << it << "\n";
    //      cout << "res: " << result << ", exp: " << exponent << "\n";
  }
  return result;
}

// member fct of isotopologuerecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopologueRecord::CalculatePartitionFctAtTempFromData_OldUnused(
    Numeric temperature) const {
  GridPosPoly gp;
  gridpos_poly(gp, mqcoeffgrid, temperature, mqcoeffinterporder);
  Vector itw;
  interpweights(itw, gp);
  return interp(itw, mqcoeff, gp);
}

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

//======================================================================
//        Functions to convert the accuracy index to ARTS units
//======================================================================

// ********* for HITRAN database *************
//convert HITRAN index for line position accuracy to ARTS
//units (Hz).

void convHitranIERF(Numeric& mdf, const Index& df) {
  switch (df) {
    case 0: {
      mdf = -1;
      break;
    }
    case 1: {
      mdf = 30 * 1E9;
      break;
    }
    case 2: {
      mdf = 3 * 1E9;
      break;
    }
    case 3: {
      mdf = 300 * 1E6;
      break;
    }
    case 4: {
      mdf = 30 * 1E6;
      break;
    }
    case 5: {
      mdf = 3 * 1E6;
      break;
    }
    case 6: {
      mdf = 0.3 * 1E6;
      break;
    }
  }
}

//convert HITRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convHitranIERSH(Numeric& mdh, const Index& dh) {
  switch (dh) {
    case 0: {
      mdh = -1;
      break;
    }
    case 1: {
      mdh = -1;
      break;
    }
    case 2: {
      mdh = -1;
      break;
    }
    case 3: {
      mdh = 30;
      break;
    }
    case 4: {
      mdh = 20;
      break;
    }
    case 5: {
      mdh = 10;
      break;
    }
    case 6: {
      mdh = 5;
      break;
    }
    case 7: {
      mdh = 2;
      break;
    }
    case 8: {
      mdh = 1;
      break;
    }
  }
  mdh = mdh / 100;
}

// ********* for MYTRAN database *************
//convert MYTRAN index for intensity and halfwidth accuracy to ARTS
//units (relative difference).
void convMytranIER(Numeric& mdh, const Index& dh) {
  switch (dh) {
    case 0: {
      mdh = 200;
      break;
    }
    case 1: {
      mdh = 100;
      break;
    }
    case 2: {
      mdh = 50;
      break;
    }
    case 3: {
      mdh = 30;
      break;
    }
    case 4: {
      mdh = 20;
      break;
    }
    case 5: {
      mdh = 10;
      break;
    }
    case 6: {
      mdh = 5;
      break;
    }
    case 7: {
      mdh = 2;
      break;
    }
    case 8: {
      mdh = 1;
      break;
    }
    case 9: {
      mdh = 0.5;
      break;
    }
  }
  mdh = mdh / 100;
}

ostream& operator<<(ostream& os, const LineshapeSpec& lsspec) {
  os << "LineshapeSpec Index: " << lsspec.Ind_ls()
     << ", NormIndex: " << lsspec.Ind_lsn() << ", Cutoff: " << lsspec.Cutoff()
     << endl;

  return os;
}

/*! cross-section replacement computer 
 *  
 * This will work as the interface for all line-by-line computations 
 * lacking special demands
 *  
 *  \param[in,out] xsec         Cross section of one tag group. This is now the
 *                              true attenuation cross section in units of m^2.
 *  \param[in,out] source       Cross section of one tag group. This is now the
 *                              true source cross section in units of m^2.
 *  \param[in,out] phase        Cross section of one tag group. This is now the
 *                              true phase cross section in units of m^2.
 *  \param[in,out] dxsec        Partial derivatives of xsec.
 *  \param[in,out] dsource      Partial derivatives of source.
 *  \param[in,out] dphase       Partial derivatives of phase.
 * 
 *  \param[in] flag_partials        Partial derivatives flags.
 *  \param[in] f_grid               Frequency grid.
 *  \param[in] abs_p                Pressure grid.
 *  \param[in] abs_t                Temperatures associated with abs_p.
 *  \param[in] abs_t_nlte           Non-lte temperatures for various energy levels.
 *  \param[in] all_vmrs             Gas volume mixing ratios [nspecies, np].
 *  \param[in] abs_species          Species tags for all species.
 *  \param[in] this_species         Index of the current species in abs_species.
 *  \param[in] abs_lines            The spectroscopic line list.
 *  \param[in] Z_DF                 The Zeeman line center shift over the magnitude of the magnetic field.
 *  \param[in] H_magntitude_Zeeman  The magnitude of the magnetic field required by Zeeman effect.
 *  \param[in] lm_p_lim             Line mixing pressure limit
 *  \param[in] isotopologue_ratios  Isotopologue ratios.
 *  \param[in] partition_functions  Partition functions.
 *  \param[in] verbosity            Verbosity level.
 * 
 *  \author Richard Larsson
 *  \date   2013-04-24
 * 
 */
void xsec_species2(Matrix& xsec,
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
                   const Matrix& abs_nlte,
                   const Matrix& all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const ArrayOfLineRecord& abs_lines,
                   const SpeciesAuxData& isotopologue_ratios,
                   const SpeciesAuxData& partition_functions) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = abs_lines.nelem();  // number of lines in the catalog
  const Index nj =
      jacobian_propmat_positions.nelem();  // number of partial derivatives
  const Index nt = source.nrows();         // number of energy levels in NLTE

  // Type of problem
  const bool do_nonlte = nt;
  const bool do_jacobi = nj;
  const bool do_temperature = do_temperature_jacobian(jacobian_quantities);

  // Test if the size of the problem is 0
  if (not np or not nf or not nl) return;

  // Move the problem to Eigen-library types
  const auto f_grid_eigen = MapToEigen(f_grid);

  // Parallelization data holder (skip if in parallel or there are too few threads)
  const Index nthread =
      arts_omp_in_parallel()
          ? 1
          : arts_omp_get_max_threads() < nl ? arts_omp_get_max_threads() : 1;
  std::vector<Eigen::VectorXcd> Fsum(nthread, Eigen::VectorXcd(nf));
  std::vector<Eigen::MatrixXcd> dFsum(nthread, Eigen::MatrixXcd(nf, nj));
  std::vector<Eigen::VectorXcd> Nsum(nthread, Eigen::VectorXcd(nf));
  std::vector<Eigen::MatrixXcd> dNsum(nthread, Eigen::MatrixXcd(nf, nj));

  for (Index ip = 0; ip < np; ip++) {
    // Constants for this level
    const Numeric& temperature = abs_t[ip];
    const Numeric& pressure = abs_p[ip];

    // Quasi-constants for this level, defined here to speed up later computations
    Index this_iso = -1;  // line isotopologue number
    Numeric t0 = -1;      // line temperature
    Numeric dc = 0, ddc_dT = 0, qt = 0, qt0 = 1,
            dqt_dT = 0;  // Doppler and partition functions

    // Line shape constants
    LineShape::Model line_shape_default_model;
    LineShape::Model& line_shape_model{line_shape_default_model};
    Vector line_shape_vmr(0);

    // Reset sum-operators
    for (Index ithread = 0; ithread < nthread; ithread++) {
      Fsum[ithread].setZero();
      dFsum[ithread].setZero();
      Nsum[ithread].setZero();
      dNsum[ithread].setZero();
    }

#pragma omp parallel for if (nthread > 1) schedule(guided, 4) \
    firstprivate(this_iso,                                    \
                 t0,                                          \
                 dc,                                          \
                 ddc_dT,                                      \
                 qt,                                          \
                 qt0,                                         \
                 dqt_dT,                                      \
                 line_shape_model,                            \
                 line_shape_vmr)
    for (Index il = 0; il < abs_lines.nelem(); il++) {
      const auto& line = abs_lines[il];

      // Local compute variables
      thread_local Eigen::VectorXcd F(nf);
      thread_local Eigen::MatrixXcd dF(nf, nj);
      thread_local Eigen::VectorXcd N(nf);
      thread_local Eigen::MatrixXcd dN(nf, nj);
      thread_local Eigen::
          Matrix<Complex, Eigen::Dynamic, Linefunctions::ExpectedDataSize()>
              data(nf, Linefunctions::ExpectedDataSize());
      thread_local Index start, nelem;

      // Partition function depends on isotopologue and line temperatures.
      // Both are commonly constant in a single catalog.  They are, however,
      // allowed to change so we must check that they do not
      if (line.Isotopologue() not_eq this_iso or line.Ti0() not_eq t0) {
        t0 = line.Ti0();

        partition_function(
            qt0,
            qt,
            t0,
            temperature,
            partition_functions.getParamType(line.Species(),
                                             line.Isotopologue()),
            partition_functions.getParam(line.Species(), line.Isotopologue()));

        if (do_temperature)
          dpartition_function_dT(dqt_dT,
                                 qt,
                                 temperature,
                                 temperature_perturbation(jacobian_quantities),
                                 partition_functions.getParamType(
                                     line.Species(), line.Isotopologue()),
                                 partition_functions.getParam(
                                     line.Species(), line.Isotopologue()));

        if (line.Isotopologue() not_eq this_iso) {
          this_iso = line.Isotopologue();
          dc = Linefunctions::DopplerConstant(temperature,
                                              line.IsotopologueData().Mass());
          if (do_temperature)
            ddc_dT = Linefunctions::dDopplerConstant_dT(temperature, dc);
        }
      }

      if (not line_shape_model.same_broadening_species(
              line.GetLineShapeModel())) {
        line_shape_model = line.GetLineShapeModel();
        line_shape_vmr = line_shape_model.vmrs(
            all_vmrs(joker, ip), abs_species, line.QuantumIdentity());
      }

      Linefunctions::set_cross_section_for_single_line(
          F,
          dF,
          N,
          dN,
          data,
          start,
          nelem,
          f_grid_eigen,
          line,
          jacobian_quantities,
          jacobian_propmat_positions,
          line_shape_vmr,
          nt ? abs_nlte(joker, ip) : Vector(0),
          pressure,
          temperature,
          dc,
          isotopologue_ratios.getParam(line.Species(), this_iso)[0].data[0],
          0.0,
          0.0,
          ddc_dT,
          qt,
          dqt_dT,
          qt0);

      auto ithread = arts_omp_get_thread_num();
      Fsum[ithread].segment(start, nelem).noalias() += F.segment(start, nelem);
      if (do_jacobi)
        dFsum[ithread].middleRows(start, nelem).noalias() +=
            dF.middleRows(start, nelem);
      if (do_nonlte)
        Nsum[ithread].segment(start, nelem).noalias() +=
            N.segment(start, nelem);
      if (do_nonlte and do_jacobi)
        dNsum[ithread].middleRows(start, nelem).noalias() +=
            dN.middleRows(start, nelem);
    }

    // Sum all the threaded results
    for (Index ithread = 0; ithread < nthread; ithread++) {
      // absorption cross-section
      MapToEigen(xsec).col(ip).noalias() += Fsum[ithread].real();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dxsec_dx[j]).col(ip).noalias() +=
            dFsum[ithread].col(j).real();

      // phase cross-section
      if (not phase.empty()) {
        MapToEigen(phase).col(ip).noalias() += Fsum[ithread].imag();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dphase_dx[j]).col(ip).noalias() +=
              dFsum[ithread].col(j).imag();
      }

      // source ratio cross-section
      if (do_nonlte) {
        MapToEigen(source).col(ip).noalias() += Nsum[ithread].real();
        for (Index j = 0; j < nj; j++)
          MapToEigen(dsource_dx[j]).col(ip).noalias() +=
              dNsum[ithread].col(j).real();
      }
    }
  }
}


const SpeciesRecord& SpeciesDataOfLines(const AbsorptionLines& lines)
{
  return global_data::species_data[lines.Species()];
}


/** Cross-section algorithm
 * 
 *  \author Richard Larsson
 *  \date   2019-10-10
 * 
 */
void xsec_species3(Matrix& xsec,
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
                   const Matrix& all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const AbsorptionLines& lines,
                   const Numeric& isot_ratio,
                   const SpeciesAuxData::AuxType& partfun_type,
                   const ArrayOfGriddedField1& partfun_data) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = lines.NumLines();  // number of lines in the catalog
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
  const Numeric QT0 = single_partition_function(lines.T0(), partfun_type, partfun_data);
  
  for (Index ip = 0; ip < np; ip++) {
    // Constants for this level
    const Numeric& temperature = abs_t[ip];
    const Numeric& pressure = abs_p[ip];
    
    // Constants for this level
    const Numeric QT = single_partition_function(temperature, partfun_type, partfun_data);
    const Numeric dQTdT = dsingle_partition_function_dT(QT, temperature, temperature_perturbation(jacobian_quantities), partfun_type, partfun_data);
    const Numeric DC = Linefunctions::DopplerConstant(temperature, lines.SpeciesMass());
    const Numeric dDCdT = Linefunctions::dDopplerConstant_dT(temperature, DC);
    const Vector line_shape_vmr = lines.BroadeningSpeciesVMR(all_vmrs(joker, ip), abs_species);
    
    Linefunctions::set_cross_section_of_band(
      scratch,
      sum,
      f_grid,
      lines,
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
      MapToEigen(dxsec_dx[j]).col(ip).noalias() +=
      sum.dF.col(j).real();
    
    // phase cross-section
    if (not phase.empty()) {
      MapToEigen(phase).col(ip).noalias() += sum.F.imag();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dphase_dx[j]).col(ip).noalias() +=
        sum.dF.col(j).imag();
    }
    
    // source ratio cross-section
    if (do_nonlte) {
      MapToEigen(source).col(ip).noalias() += sum.N.real();
      for (Index j = 0; j < nj; j++)
        MapToEigen(dsource_dx[j]).col(ip).noalias() +=
        sum.dN.col(j).real();
    }
  }
}
