/* Copyright (C) 2019
 Richard Larsson <ric.larsson@gmail.com>

 
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

/** Contains the absorption namespace
 * @file   absorptionlines.cc
 * @author Richard Larsson
 * @date   2019-09-07
 * 
 * @brief  Contains the absorption lines implementation
 * 
 * This namespace contains classes to deal with absorption lines
 **/

#include "absorptionlines.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <ostream>
#include <string>

#include "absorption.h"
#include "arts_conversions.h"
#include "debug.h"
#include "enums.h"
#include "file.h"
#include "hitran_species.h"
#include "jpl_species.h"
#include "linescaling.h"
#include "lineshapemodel.h"
#include "matpack_math.h"
#include "quantum_numbers.h"
#include "rational.h"
#include "wigner_functions.h"

LineShape::Output
Absorption::Lines::ShapeParameters(size_t k, Numeric T, Numeric P,
                                   const Vector &vmrs) const ARTS_NOEXCEPT {
  auto &lineshape = lines[k].lineshape;

  using namespace LineShape;
  return Output{lineshape.G0(T, T0, P, vmrs),  lineshape.D0(T, T0, P, vmrs),
                lineshape.G2(T, T0, P, vmrs),  lineshape.D2(T, T0, P, vmrs),
                lineshape.FVC(T, T0, P, vmrs), lineshape.ETA(T, T0, P, vmrs),
                lineshape.Y(T, T0, P, vmrs),   lineshape.G(T, T0, P, vmrs),
                lineshape.DV(T, T0, P, vmrs)}
      .no_linemixing(not DoLineMixing(P));
}

LineShape::Output Absorption::Lines::ShapeParameters(size_t k,
                                                     Numeric T,
                                                     Numeric P,
                                                     size_t pos) const ARTS_NOEXCEPT {
  auto &lineshape = lines[k].lineshape[pos];

  return lineshape.at(T, T0, P).no_linemixing(not DoLineMixing(P));
}

LineShape::Output Absorption::Lines::ShapeParameters_dT(
    size_t k, Numeric T, Numeric P, const Vector& vmrs) const ARTS_NOEXCEPT {
  auto &lineshape = lines[k].lineshape;

  using namespace LineShape;
  return Output{lineshape.dG0dT(T, T0, P, vmrs),  lineshape.dD0dT(T, T0, P, vmrs),
                lineshape.dG2dT(T, T0, P, vmrs),  lineshape.dD2dT(T, T0, P, vmrs),
                lineshape.dFVCdT(T, T0, P, vmrs), lineshape.dETAdT(T, T0, P, vmrs),
                lineshape.dYdT(T, T0, P, vmrs),   lineshape.dGdT(T, T0, P, vmrs),
                lineshape.dDVdT(T, T0, P, vmrs)}
      .no_linemixing(not DoLineMixing(P));
}

LineShape::Output Absorption::Lines::ShapeParameters_dT(
    size_t k, Numeric T, Numeric P, size_t pos) const ARTS_NOEXCEPT {
  auto &lineshape = lines[k].lineshape[pos];

  return lineshape.dT(T, T0, P).no_linemixing(not DoLineMixing(P));
}

Index Absorption::Lines::LineShapePos(
    const Species::Species spec) const ARTS_NOEXCEPT {
  // Is always first if this is self and self broadening exists
  if (selfbroadening and spec == quantumidentity.Species()) {
    return 0;
  }

  // First and last might be artificial so they should not be checked
  const Index s = selfbroadening;
  const Index e = broadeningspecies.nelem() - bathbroadening;
  for (Index i = s; i < e; i++) {
    if (spec == broadeningspecies[i]) {
      return i;
    }
  }

  // At this point, the ID is not explicitly among the broadeners, but bath broadening means its VMR still might matter
  if (bathbroadening) return broadeningspecies.nelem() - 1;
  return -1;
}

LineShape::Output Absorption::Lines::ShapeParameters_dVMR(
    size_t k,
    Numeric T,
    Numeric P,
    const QuantumIdentifier& vmr_qid) const ARTS_NOEXCEPT {
  auto &lineshape = lines[k].lineshape;

  const Index pos = LineShapePos(vmr_qid.Species());

  LineShape::Output out{};
  if (pos >= 0) {
    out = lineshape[pos].at(T, T0, P);
    const Index bath = lineshape.nelem() - 1;
    if (bathbroadening and pos not_eq bath) {
      out -= lineshape[bath].at(T, T0, P);
    }
  }
  return out;
}

Absorption::SingleLineExternal Absorption::ReadFromArtscat3Stream(istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

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

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0) {
    // Set the species
    const auto isotopologue = Species::Tag(Species::update_isot_name(artsid));
    ARTS_USER_ERROR_IF(
        isotopologue.is_joker() or
            isotopologue.type not_eq Species::TagType::Plain,
        "A line catalog species can only be of the form \"Plain\", meaning it\nhas the form SPECIES-ISONUM.\n"
        "Your input contains: ",
        artsid,
        ". which we cannot interpret as a plain species")
    data.quantumidentity = QuantumIdentifier{
        Species::find_species_index(isotopologue.Isotopologue())};

    // Extract center frequency:
    icecream >> double_imanip() >> data.line.F0;

    Numeric psf;
    // Extract pressure shift:
    icecream >> double_imanip() >> psf;

    // Extract intensity:
    icecream >> double_imanip() >> data.line.I0;

    // Extract reference temperature for Intensity in K:
    icecream >> double_imanip() >> data.T0;

    // Extract lower state energy:
    icecream >> double_imanip() >> data.line.E0;

    // Extract air broadening parameters:
    Numeric agam, sgam;
    icecream >> double_imanip() >> agam;
    icecream >> double_imanip() >> sgam;

    // Extract temperature coefficient of broadening parameters:
    Numeric nair, nself;
    icecream >> double_imanip() >> nair;
    icecream >> double_imanip() >> nself;

    // Extract reference temperature for broadening parameter in K:
    Numeric tgam;
    icecream >> double_imanip() >> tgam;

    // Extract the aux parameters:
    Index naux;
    icecream >> naux;

    // resize the aux array and read it
    ArrayOfNumeric maux;
    maux.resize(naux);

    for (Index j = 0; j < naux; j++) {
      icecream >> double_imanip() >> maux[j];
      //cout << "maux" << j << " = " << maux[j] << "\n";
    }

    // Extract accuracies:
    try {
      Numeric unused_numeric;
      icecream >> double_imanip() >> unused_numeric /*mdf*/;
      icecream >> double_imanip() >> unused_numeric /*mdi0*/;
      icecream >> double_imanip() >> unused_numeric /*dagam*/;
      icecream >> double_imanip() >> unused_numeric /*dsgam*/;
      icecream >> double_imanip() >> unused_numeric /*dnair*/;
      icecream >> double_imanip() >> unused_numeric /*dnself*/;
      icecream >> double_imanip() >> unused_numeric /*dpsf*/;
    } catch (const std::runtime_error&) {
    }

    // Fix if tgam is different from ti0
    if (tgam != data.T0) {
      agam = agam * pow(tgam / data.T0, nair);
      sgam = sgam * pow(tgam / data.T0, nself);
      psf = psf * pow(tgam / data.T0, (Numeric).25 + (Numeric)1.5 * nair);
    }

    // Set line shape computer
    data.line.lineshape = LineShape::Model(sgam, nself, agam, nair, psf);
  }

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromArtscat4Stream(istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = false;
  data.lineshapetype = LineShape::Type::VP;

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

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

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0) {
    // Set line ID
    const auto isotopologue = Species::Tag(Species::update_isot_name(artsid));
    ARTS_USER_ERROR_IF(
        isotopologue.is_joker() or
            isotopologue.type not_eq Species::TagType::Plain,
        "A line catalog species can only be of the form \"Plain\", meaning it\nhas the form SPECIES-ISONUM.\n"
        "Your input contains: ",
        artsid,
        ". which we cannot interpret as a plain species")
    data.quantumidentity = QuantumIdentifier{
        Species::find_species_index(isotopologue.Isotopologue())};

    // Extract center frequency:
    icecream >> double_imanip() >> data.line.F0;

    // Extract intensity:
    icecream >> double_imanip() >> data.line.I0;

    // Extract reference temperature for Intensity in K:
    icecream >> double_imanip() >> data.T0;

    // Extract lower state energy:
    icecream >> double_imanip() >> data.line.E0;

    // Extract Einstein A-coefficient:
    icecream >> double_imanip() >> data.line.A;

    // Extract upper state stat. weight:
    icecream >> double_imanip() >> data.line.gupp;

    // Extract lower state stat. weight:
    icecream >> double_imanip() >> data.line.glow;

    LineShape::from_artscat4(icecream,
                             data.lineshapetype,
                             data.selfbroadening,
                             data.bathbroadening,
                             data.line.lineshape,
                             data.species,
                             data.quantumidentity);
  }

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromArtscat5Stream(istream& is) {
  // Default data and values for this type
  SingleLineExternal data;

  LineShape::Model line_mixing_model;
  bool lmd_found = false;

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

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

  try {
    String artsid;
    icecream >> artsid;

    if (artsid.length() != 0) {
      // Set line ID:
      const auto isotopologue = Species::Tag(Species::update_isot_name(artsid));
      ARTS_USER_ERROR_IF(
          isotopologue.is_joker() or
              isotopologue.type not_eq Species::TagType::Plain,
          "A line catalog species can only be of the form \"Plain\", meaning it\nhas the form SPECIES-ISONUM.\n"
          "Your input contains: ",
          artsid,
          ". which we cannot interpret as a plain species")
      data.quantumidentity = QuantumIdentifier{
          Species::find_species_index(isotopologue.Isotopologue())};

      // Extract center frequency:
      icecream >> double_imanip() >> data.line.F0;

      // Extract intensity:
      icecream >> double_imanip() >> data.line.I0;

      // Extract reference temperature for Intensity in K:
      icecream >> double_imanip() >> data.T0;

      // Extract lower state energy:
      icecream >> double_imanip() >> data.line.E0;

      // Extract Einstein A-coefficient:
      icecream >> double_imanip() >> data.line.A;

      // Extract upper state stat. weight:
      icecream >> double_imanip() >> data.line.gupp;

      // Extract lower state stat. weight:
      icecream >> double_imanip() >> data.line.glow;

      String token;
      Index nelem;
      icecream >> token;

      while (icecream) {
        // Read pressure broadening (LEGACY)
        if (token == "PB") {
          LineShape::from_pressurebroadeningdata(icecream,
                                                 data.lineshapetype,
                                                 data.selfbroadening,
                                                 data.bathbroadening,
                                                 data.line.lineshape,
                                                 data.species,
                                                 data.quantumidentity);
          icecream >> token;
        } else if (token == "QN") {
          // Quantum numbers

          icecream >> token;
          ARTS_USER_ERROR_IF(
              token != "UP", "Unknown quantum number tag: ", token)

          icecream >> token;
          String r;
          while (icecream and token not_eq "LO") {
            auto qn = Quantum::Number::toType(token);
            icecream >> r;
            icecream >> token;

            if (good_enum(qn)) {
              if (data.quantumidentity.val.has(qn)) {
                Quantum::Number::Value val = data.quantumidentity.val[qn];
                val.set(r, true);
                data.quantumidentity.val.set(val);
              } else {
                data.quantumidentity.val.add(qn).set(r, true);
              }
            }
          }

          ARTS_USER_ERROR_IF(
              !is || token != "LO",
              "Error in catalog. Lower quantum number tag 'LO' not found.")

          icecream >> token;
          auto qn = Quantum::Number::toType(token);
          while (icecream and
                 not(token == "LM" or token == "LF" or token == "ZM" or
                     token == "LSM" or token == "PB")) {
            icecream >> r;
            icecream >> token;

            if (good_enum(qn)) {
              if (data.quantumidentity.val.has(qn)) {
                Quantum::Number::Value val = data.quantumidentity.val[qn];
                val.set(r, false);
                data.quantumidentity.val.set(val);
              } else {
                data.quantumidentity.val.add(qn).set(r, false);
              }
            }

            qn = Quantum::Number::toType(token);
          }
        } else if (token == "LM") {
          LineShape::from_linemixingdata(icecream, line_mixing_model);
          icecream >> token;
          lmd_found = true;
        } else if (token == "LF") {
          LineShape::from_linefunctiondata(icecream,
                                           data.lineshapetype,
                                           data.selfbroadening,
                                           data.bathbroadening,
                                           data.line.lineshape,
                                           data.species);
          // if (data.selfbroadening) data.species[0] = data.quantumidentity.Species();
          icecream >> token;
        } else if (token == "ZM") {
          // Zeeman effect
          icecream >> data.line.zeeman;
          icecream >> token;
        } else if (token == "LSM") {
          // Line shape modifications

          // Starts with the number of modifications
          icecream >> nelem;
          for (Index lsm = 0; lsm < nelem; lsm++) {
            icecream >> token;

            // cutoff frequency
            if (token == "CUT") {
              icecream >> double_imanip() >> data.cutofffreq;
              data.cutoff = CutoffType::ByLine;
            }

            // linemixing pressure limit
            if (token == "LML") {
              icecream >> double_imanip() >> data.linemixinglimit;
            }

            // mirroring
            else if (token == "MTM") {
              icecream >> data.mirroring;
            }

            // line normalization
            else if (token == "LNT") {
              icecream >> data.normalization;
            }

            else {
              ARTS_USER_ERROR("Unknown line modifications given: ", token)
            }
          }
          icecream >> token;
        } else {
          ARTS_USER_ERROR("Unknown line data tag in legacy reading routine: ",
                          token)
        }
      }
    }
  } catch (const std::runtime_error& e) {
    ARTS_USER_ERROR("Parse error in catalog line: \n", line, '\n', e.what())
  }

  if (lmd_found)
    data.line.lineshape.SetLineMixingModel(line_mixing_model.Data()[0]);

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromHitran2004Stream(
    istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);
  data.species[1] = Species::Species::Bath;

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    comment = false;
  }

  // Extract isotopologue:
  char iso;
  extract(iso, line, 1);

  // Set line data
  data.quantumidentity = Hitran::id_from_lookup(mo, iso);
  data.species[0] = data.quantumidentity.Species();

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    constexpr Numeric w2Hz = Constant::c * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    data.line.F0 = v * w2Hz;
  }

  // Intensity.
  {
    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    constexpr Numeric hi2arts = 1e-2 * Constant::c;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s, line, 10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    data.line.I0 = s * hi2arts;
    // Take out isotopologue ratio:
    data.line.I0 /= Hitran::ratio_from_lookup(mo, iso);
  }

  // Einstein coefficient
  { extract(data.line.A, line, 10); }

  // Air broadening parameters.
  Numeric agam, sgam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(data.line.E0, line, 10);

    // Convert to Joule:
    data.line.E0 = wavenumber_to_joule(data.line.E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  // Upper state global quanta
  {
    Index eu;
    extract(eu, line, 15);
  }

  // Lower state global quanta
  {
    Index el;
    extract(el, line, 15);
  }

  // Upper state local quanta
  {
    Index eul;
    extract(eul, line, 15);
  }

  // Lower state local quanta
  {
    Index ell;
    extract(ell, line, 15);
  }

  // Accuracy index for frequency
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for air-broadened halfwidth
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for self-broadened half-width
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for temperature-dependence exponent for agam
  {
    Index dn;
    // Extract HITRAN value:
    extract(dn, line, 1);
  }

  // Accuracy index for pressure shift
  {
    Index dpsfi;
    // Extract HITRAN value (given in cm-1):
    extract(dpsfi, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN 2004. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  data.T0 = 296.0;

  // Set line shape computer
  data.line.lineshape = LineShape::hitran_model(sgam, nself, agam, nair, psf);
  {
    Index garbage;
    extract(garbage, line, 13);

    // The statistical weights
    extract(data.line.gupp, line, 7);
    extract(data.line.glow, line, 7);
  }

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromHitranOnlineStream(
    istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);
  data.species[1] = Species::Species::Bath;

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    comment = false;
  }

  // Extract isotopologue:
  char iso;
  extract(iso, line, 1);

  // Set line data
  data.quantumidentity = Hitran::id_from_lookup(mo, iso);
  data.species[0] = data.quantumidentity.Species();

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    constexpr Numeric w2Hz = Constant::c * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    data.line.F0 = v * w2Hz;
  }

  // Intensity.
  {
    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    constexpr Numeric hi2arts = 1e-2 * Constant::c;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s, line, 10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    data.line.I0 = s * hi2arts;
    // Take out isotopologue ratio:
    data.line.I0 /= Hitran::ratio_from_lookup(mo, iso);
  }

  // Einstein coefficient
  { extract(data.line.A, line, 10); }

  // Air broadening parameters.
  Numeric agam, sgam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(data.line.E0, line, 10);

    // Convert to Joule:
    data.line.E0 = wavenumber_to_joule(data.line.E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  // Upper state global quanta
  {
    Index eu;
    extract(eu, line, 15);
  }

  // Lower state global quanta
  {
    Index el;
    extract(el, line, 15);
  }

  // Upper state local quanta
  {
    Index eul;
    extract(eul, line, 15);
  }

  // Lower state local quanta
  {
    Index ell;
    extract(ell, line, 15);
  }

  // Accuracy index for frequency
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for air-broadened halfwidth
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for self-broadened half-width
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // Accuracy index for temperature-dependence exponent for agam
  {
    Index dn;
    // Extract HITRAN value:
    extract(dn, line, 1);
  }

  // Accuracy index for pressure shift
  {
    Index dpsfi;
    // Extract HITRAN value (given in cm-1):
    extract(dpsfi, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN 2004. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  data.T0 = 296.0;

  // Set line shape computer
  data.line.lineshape = LineShape::hitran_model(sgam, nself, agam, nair, psf);
  {
    Index garbage;
    extract(garbage, line, 13);

    // The statistical weights
    extract(data.line.gupp, line, 7);
    extract(data.line.glow, line, 7);
  }

  // ADD QUANTUM NUMBER PARSING HERE!
  String upper, lower;
  std::stringstream ss;
  ss.str(line);
  ss >> upper >> lower;
  data.quantumidentity.val = Quantum::Number::from_hitran(upper, lower);

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromHitran2001Stream(
    istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);
  data.species[1] = Species::Species::Bath;

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    comment = false;
  }

  // Extract isotopologue:
  char iso;
  extract(iso, line, 1);

  // Set line data
  data.quantumidentity = Hitran::id_from_lookup(mo, iso);
  data.species[0] = data.quantumidentity.Species();

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    constexpr Numeric w2Hz = Constant::c * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    data.line.F0 = v * w2Hz;
    //    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    constexpr Numeric hi2arts = 1e-2 * Constant::c;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s, line, 10);
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    data.line.I0 = s * hi2arts;
    // Take out isotopologue ratio:
    data.line.I0 /= Hitran::ratio_from_lookup(mo, iso);
  }

  // Skip transition probability:
  {
    Numeric r;
    extract(r, line, 10);
  }

  // Air broadening parameters.
  Numeric agam, sgam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(data.line.E0, line, 10);

    // Convert to Joule:
    data.line.E0 = wavenumber_to_joule(data.line.E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  //Skip upper state global quanta index
  {
    Index eu;
    extract(eu, line, 3);
  }

  //Skip lower state global quanta index
  {
    Index el;
    extract(el, line, 3);
  }

  //Skip upper state local quanta
  {
    Index eul;
    extract(eul, line, 9);
  }

  //Skip lower state local quanta
  {
    Index ell;
    extract(ell, line, 9);
  }

  // Accuracy index for frequency reference
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity reference
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for halfwidth reference
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  data.T0 = 296.0;

  // Set line shape computer
  data.line.lineshape = LineShape::hitran_model(sgam, nself, agam, nair, psf);

  // That's it!
  data.bad = false;
  return data;
}

Absorption::SingleLineExternal Absorption::ReadFromLBLRTMStream(istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // The first item is the molecule number:
  Index mo;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);
    if (line[0] == '>' or line[0] == '%') continue;

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // If the catalogue is in dos encoding, throw away the
    // additional carriage return
    if (line[line.nelem() - 1] == 13) {
      line.erase(line.nelem() - 1, 1);
    }

    mo = 0;
    // Initialization of mo is important, because mo stays the same
    // if line is empty.
    extract(mo, line, 2);
    comment = false;
  }

  // Extract isotopologue:
  char iso;
  extract(iso, line, 1);

  // Set line data
  data.quantumidentity = Hitran::id_from_lookup(mo, iso);

  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    constexpr Numeric w2Hz = Constant::c * 100.;

    // Extract HITRAN postion:
    extract(v, line, 12);

    // ARTS position in Hz:
    data.line.F0 = v * w2Hz;
    //    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function).
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c).
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4).
    // 3. Take out the isotopologue ratio.

    constexpr Numeric hi2arts = 1e-2 * Constant::c;

    Numeric s;
    if (line[6] == 'D') line[6] = 'E';
    // Extract HITRAN intensity:
    extract(s,
            line,
            10);  // NOTE:  If error shooting, FORTRAN "D" is not read properly.
    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    data.line.I0 = s * hi2arts;
    // Take out isotopologue ratio:
    data.line.I0 /= Hitran::ratio_from_lookup(mo, iso);
  }

  // Skip transition probability:
  {
    Numeric r;
    extract(r, line, 10);
  }

  // Air broadening parameters.
  Numeric sgam, agam;
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN AGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    agam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam, line, 5);

    // ARTS parameter in Hz/Pa:
    sgam = gam * hi2arts;

    // If zero, set to agam:
    if (0 == sgam) sgam = agam;

    //    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }

  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(data.line.E0, line, 10);

    // Convert to Joule:
    data.line.E0 = wavenumber_to_joule(data.line.E0);
  }

  // Temperature coefficient of broadening parameters.
  Numeric nair, nself;
  {
    // This is dimensionless, we can also extract directly.
    extract(nair, line, 4);

    // Set self broadening temperature coefficient to the same value:
    nself = nair;
    //    cout << "mnair = " << mnair << endl;
  }

  // Pressure shift.
  Numeric psf;
  {
    // HITRAN value in cm^-1 / atm. So the conversion goes exactly as
    // for the broadening parameters.
    Numeric d;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    constexpr Numeric w2Hz = Constant::c * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    constexpr Numeric hi2arts = w2Hz / Conversion::atm2pa(1);

    // Extract HITRAN value:
    extract(d, line, 8);

    // ARTS value in Hz/Pa
    psf = d * hi2arts;
  }
  // Set the accuracies using the definition of HITRAN
  // indices. If some are missing, they are set to -1.

  //Skip upper state global quanta index
  {
    Index eu;
    extract(eu, line, 3);
  }

  //Skip lower state global quanta index
  {
    Index el;
    extract(el, line, 3);
  }

  //Skip upper state local quanta
  {
    Index eul;
    extract(eul, line, 9);
  }

  //Skip lower state local quanta
  {
    Index ell;
    extract(ell, line, 9);
  }

  // Accuracy index for frequency reference
  {
    Index df;
    // Extract HITRAN value:
    extract(df, line, 1);
  }

  // Accuracy index for intensity reference
  {
    Index di0;
    // Extract HITRAN value:
    extract(di0, line, 1);
  }

  // Accuracy index for halfwidth reference
  {
    Index dgam;
    // Extract HITRAN value:
    extract(dgam, line, 1);
  }

  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for HITRAN)
  data.T0 = 296.0;

  // Skip four
  {
    Index four;
    extract(four, line, 4);
  }

  // This is the test for the last two characters of the line
  {
    /* 
     *   0 is nothing, 
     *  -1 is linemixing on the next line, 
     *  -3 is the non-resonant line 
     */
    Index test;
    extract(test, line, 2);
    //If the tag is as it should be, then a minus one means that more should be read
    if (test == -1 || test == -3)
      getline(is, line);
    else  // the line is done and we are happy to leave
    {
      data.line.lineshape =
          LineShape::hitran_model(sgam, nself, agam, nair, psf);

      data.bad = false;
      return data;
    }
  }

  // In case we are unable to leave, the next line is a line mixing parameter line

  // First is the molecular number.  This should be the same as above.
  {
    Index mo2;
    extract(mo2, line, 2);
    // Skip one

    ARTS_USER_ERROR_IF(mo != mo2, "There is an error in the line mixing\n");
  }

  Vector Y(4), G(4), T(4);

  // These are constants for AER but should be included because we need their grid.
  T[0] = 200;
  T[1] = 250;
  T[2] = 296;
  T[3] = 340;

  // Next is the Y  and G at various temperatures
  {
    Numeric Y_200K;
    extract(Y_200K, line, 13);
    Y[0] = Y_200K;
  }
  {
    Numeric G_200K;
    extract(G_200K, line, 11);
    G[0] = G_200K;
  }
  {
    Numeric Y_250K;
    extract(Y_250K, line, 13);
    Y[1] = Y_250K;
  }
  {
    Numeric G_250K;
    extract(G_250K, line, 11);
    G[1] = G_250K;
  }
  {
    Numeric Y_296K;
    extract(Y_296K, line, 13);
    Y[2] = Y_296K;
  }
  {
    Numeric G_296K;
    extract(G_296K, line, 11);
    G[2] = G_296K;
  }
  {
    Numeric Y_340K;
    extract(Y_340K, line, 13);
    Y[3] = Y_340K;
  }
  {
    Numeric G_340K;
    extract(G_340K, line, 11);
    G[3] = G_340K;
  }

  // Cpnvert from per Atm and per Atm^2
  Y /= Conversion::atm2pa(1);
  G /= Math::pow2(Conversion::atm2pa(1));

  // ARTS uses (1-iY) as line-mixing factor, LBLRTM uses (1+iY), so we must change sign
  Y *= -1;

  // Test that this is the end
  {
    Index test;
    extract(test, line, 2);
    if (test == -1) {
      data.line.lineshape = LineShape::lblrtm_model(sgam,
                                                    nself,
                                                    agam,
                                                    nair,
                                                    psf,
                                                    {T[0],
                                                     T[1],
                                                     T[2],
                                                     T[3],
                                                     Y[0],
                                                     Y[1],
                                                     Y[2],
                                                     Y[3],
                                                     G[0],
                                                     G[1],
                                                     G[2],
                                                     G[3]});

      data.bad = false;
      return data;
    }
    if (test == -3) {
      data.line.lineshape = LineShape::lblrtm_model(sgam,
                                                    nself,
                                                    agam,
                                                    nair,
                                                    psf,
                                                    {T[0],
                                                     T[1],
                                                     T[2],
                                                     T[3],
                                                     Y[0],
                                                     Y[1],
                                                     Y[2],
                                                     Y[3],
                                                     G[0],
                                                     G[1],
                                                     G[2],
                                                     G[3]});

      data.bad = false;
      return data;
    }
    return data;
  }
}

std::vector<Absorption::Lines> Absorption::split_list_of_external_lines(
    std::vector<SingleLineExternal>& external_lines,
    const std::vector<QuantumNumberType>& localquantas,
    const std::vector<QuantumNumberType>& globalquantas) {
  std::vector<Lines> lines(0);

  // Loop but make copies because we will need to modify some of the data
  while (external_lines.size()) {
    auto& sle = external_lines.back();

    // Adapt broadening to fit with line catalog
    if (sle.selfbroadening) sle.species.front() = sle.quantumidentity.Species();
    if (sle.bathbroadening) sle.species.back() = Species::Species::Bath;

    Quantum::Number::GlobalState global_id{
        sle.quantumidentity.isotopologue_index};
    for (auto qn : globalquantas)
      if (sle.quantumidentity.val.has(qn))
        global_id.val.set(sle.quantumidentity.val[qn]);

    Quantum::Number::LocalState local_id{};
    for (auto qn : localquantas)
      if (sle.quantumidentity.val.has(qn))
        local_id.val.set(sle.quantumidentity.val[qn]);

    // Set the line
    auto line = sle.line;
    line.localquanta = local_id;

    // Either find a line like this in the list of lines or start a new Lines
    auto band = std::find_if(lines.begin(), lines.end(), [&](const Lines& li) {
      return li.MatchWithExternal(sle, global_id);
    });
    if (band not_eq lines.end()) {
      band->AppendSingleLine(line);
    } else {
      lines.push_back(Lines(sle.selfbroadening,
                            sle.bathbroadening,
                            sle.cutoff,
                            sle.mirroring,
                            sle.population,
                            sle.normalization,
                            sle.lineshapetype,
                            sle.T0,
                            sle.cutofffreq,
                            sle.linemixinglimit,
                            global_id,
                            sle.species,
                            {line}));
    }
    external_lines.pop_back();
  }

  return lines;
}

namespace Absorption {
std::ostream& operator<<(std::ostream& os, const Absorption::Lines& lines) {
  for (auto& line : lines.lines) os << line << '\n';
  return os;
}

std::istream& operator>>(std::istream& is, Lines& lines) {
  for (auto& line : lines.lines) is >> line;
  return is;
}

std::ostream& operator<<(std::ostream& os, const Absorption::SingleLine& line) {
  os << line.F0 << ' ' << line.I0 << ' ' << line.E0 << ' ' << line.glow << ' '
     << line.gupp << ' ' << line.A << ' ' << line.zeeman << ' '
     << line.lineshape;
  if (line.localquanta.val.nelem()) os << ' ' << line.localquanta.values();
  return os;
}

std::istream& operator>>(std::istream& is, Absorption::SingleLine& line) {
  is >> double_imanip() >> line.F0 >> line.I0 >> line.E0 >> line.glow >>
      line.gupp >> line.A;
  return is >> line.zeeman >> line.lineshape >> line.localquanta;
}
}  // namespace Absorption

String Absorption::Lines::SpeciesName() const noexcept {
  return quantumidentity.Isotopologue().FullName();
}

String Absorption::Lines::MetaData() const {
  std::ostringstream os;

  os << "\nLines meta-data:\n";
  os << '\t' << "Species identity:\n";
  os << "\t\tSpecies: " << SpeciesName() << '\n';
  os << "\t\tIdentity: " << quantumidentity << '\n';
  os << '\t' << cutofftype2metadatastring(cutoff, cutofffreq);
  os << '\t' << populationtype2metadatastring(population);
  os << '\t' << normalizationtype2metadatastring(normalization);
  os << '\t' << LineShape::shapetype2metadatastring(lineshapetype);
  os << '\t' << mirroringtype2metadatastring(mirroring);
  os << '\t' << "The reference temperature for all line parameters is " << T0
     << " K.\n";
  if (linemixinglimit < 0)
    os << '\t' << "If applicable, there is no line mixing limit.\n";
  else
    os << '\t' << "If applicable, there is a line mixing limit at "
       << linemixinglimit << " Pa.\n";

  if (not NumLines()) {
    os << "\tNo line data is available.\n";
  } else {
    os << "\tThere are " << NumLines() << " lines available.\n";

    auto& line = lines.front();
    os << "\tThe front line has:\n";
    os << "\t\t"
       << "f0: " << line.F0 << " Hz\n";
    os << "\t\t"
       << "i0: " << line.I0 << " m^2/Hz\n";
    os << "\t\t"
       << "e0: " << line.E0 << " J\n";
    os << "\t\t"
       << "Lower stat. weight: " << line.glow << " [-]\n";
    os << "\t\t"
       << "Upper stat. weight: " << line.gupp << " [-]\n";
    os << "\t\t"
       << "A: " << line.A << " 1/s\n";
    os << "\t\t"
       << "Zeeman splitting of lower state: " << line.zeeman.gl() << " [-]\n";
    os << "\t\t"
       << "Zeeman splitting of upper state: " << line.zeeman.gu() << " [-]\n";
    os << "\t\t"
       << "Local quantum numbers: " << line.localquanta.val << "\n";

    ArrayOfString ls_meta = LineShape::ModelMetaDataArray(
        line.lineshape, selfbroadening, broadeningspecies, T0);
    os << "\t\t"
       << "Line shape parameters (are normalized by sum(VMR)):\n";
    for (auto& ls_form : ls_meta) os << "\t\t\t" << ls_form << "\n";
  }

  return os.str();
}

void Absorption::Lines::RemoveLine(Index i) noexcept {
  lines.erase(lines.begin() + i);
}

Absorption::SingleLine Absorption::Lines::PopLine(Index i) noexcept {
  auto line = lines[i];
  RemoveLine(i);
  return line;
}

void Absorption::Lines::ReverseLines() noexcept {
  std::reverse(lines.begin(), lines.end());
}

Numeric Absorption::Lines::SpeciesMass() const noexcept {
  return quantumidentity.Isotopologue().mass;
}

Vector Absorption::Lines::BroadeningSpeciesVMR(
    const AtmPoint& atm_point) const {
  return LineShape::vmrs(atm_point, broadeningspecies);
}

Vector Absorption::Lines::BroadeningSpeciesVMR(
    const ConstVectorView& atm_vmrs,
    const ArrayOfArrayOfSpeciesTag& atm_spec) const {
  return LineShape::vmrs(atm_vmrs, atm_spec, broadeningspecies);
}

Vector Absorption::Lines::BroadeningSpeciesMass(
    const ConstVectorView& atm_vmrs,
    const ArrayOfArrayOfSpeciesTag& atm_spec,
    const SpeciesIsotopologueRatios& ir,
    const Numeric& bath_mass) const {
  Vector mass = LineShape::mass(atm_vmrs, atm_spec, broadeningspecies, ir);
  if (bathbroadening and bath_mass > 0) mass[mass.nelem() - 1] = bath_mass;
  return mass;
}

Numeric Absorption::Lines::SelfVMR(
    const ConstVectorView& atm_vmrs,
    const ArrayOfArrayOfSpeciesTag& atm_spec) const {
  ARTS_USER_ERROR_IF(atm_vmrs.nelem() not_eq atm_spec.nelem(),
                     "Bad species and vmr lists");

  for (Index i = 0; i < atm_spec.nelem(); i++)
    if (atm_spec[i].nelem() and atm_spec[i].Species() == Species())
      return atm_vmrs[i];
  return 0;
}

Numeric Absorption::reduced_rovibrational_dipole(
    Rational Jf, Rational Ji, Rational lf, Rational li, Rational k) {
  if (not iseven(Jf + lf + 1))
    return -sqrt(2 * Jf + 1) * wigner3j(Jf, k, Ji, li, lf - li, -lf);
  return +sqrt(2 * Jf + 1) * wigner3j(Jf, k, Ji, li, lf - li, -lf);
}

Numeric Absorption::reduced_magnetic_quadrapole(Rational Jf,
                                                Rational Ji,
                                                Rational N) {
  if (not iseven(Jf + N))
    return -sqrt(6 * (2 * Jf + 1) * (2 * Ji + 1)) *
           wigner6j(1, 1, 1, Ji, Jf, N);
  return +sqrt(6 * (2 * Jf + 1) * (2 * Ji + 1)) * wigner6j(1, 1, 1, Ji, Jf, N);
}

Absorption::SingleLineExternal Absorption::ReadFromJplStream(istream& is) {
  // Default data and values for this type
  SingleLineExternal data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  data.species.resize(2);

  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // Because of the fixed FORTRAN format, we need to break up the line
    // explicitly in apropriate pieces. Not elegant, but works!

    // Extract center frequency:
    // Initialization of v is important, because v stays the same
    // if line is empty.
    // JPL position in MHz:
    Numeric v = 0.0;

    // Extract JPL position:
    extract(v, line, 13);

    // check for empty line
    if (v != 0.0) {
      // ARTS position in Hz:
      data.line.F0 = v * 1E6;

      comment = false;
    }
  }

  // Accuracy for line position
  {
    Numeric df;
    extract(df, line, 8);
  }

  // Intensity.
  {
    // JPL has log (10) of intensity in nm2 MHz at 300 Kelvin.
    //
    // We need to do the following:
    // 1. take 10^intensity
    // 2. convert to cm-1/(molecule * cm-2): devide by c * 1e10
    // 3. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 4. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)

    Numeric s;

    // Extract JPL intensity:
    extract(s, line, 8);

    // remove log
    s = pow((Numeric)10., s);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    data.line.I0 = s / 1E12;
  }

  // Degrees of freedom
  {
    Index dr;

    // Extract degrees of freedom
    extract(dr, line, 2);
  }

  // Lower state energy.
  {
    // JPL parameter is in wavenumbers (cm^-1).
    // We have to convert this to the ARTS unit Joule.

    // Extract from Catalogue line
    extract(data.line.E0, line, 10);

    // Convert to Joule:
    data.line.E0 = wavenumber_to_joule(data.line.E0);
  }

  // Upper state degeneracy
  {
    Index gup;

    // Extract upper state degeneracy
    extract(gup, line, 3);
  }

  // Tag number
  Index tag;
  {
    // Extract Tag number
    extract(tag, line, 7);

    // make sure tag is not negative (damned jpl cat):
    tag = tag > 0 ? tag : -tag;
  }

  // Set line ID
  data.quantumidentity = Jpl::id_from_lookup(tag);

  // Air broadening parameters: unknown to jpl, use old iup forward
  // model default values, which is mostly set to 0.0025 GHz/hPa, even
  // though for some lines the pressure broadening is given explicitly
  // in the program code. The explicitly given values are ignored and
  // only the default value is set. Self broadening was in general not
  // considered in the old forward model.
  Numeric agam, sgam;
  {
    // ARTS parameter in Hz/Pa:
    agam = 2.5E4;

    // ARTS parameter in Hz/Pa:
    sgam = agam;
  }

  // Temperature coefficient of broadening parameters. Was set to 0.75
  // in old forward model, even though for some lines the parameter is
  // given explicitly in the program code. The explicitly given values
  // are ignored and only the default value is set. Self broadening
  // not considered.
  Numeric nair, nself;
  {
    nair = 0.75;
    nself = 0.0;
  }

  // Reference temperature for broadening parameter in K, was
  // generally set to 300 K in old forward model, with the exceptions
  // as already mentioned above: //DEPRECEATED but is same as for mti0 so moving on
  //   {
  //     mtgam = 300.0;
  //   }

  // Pressure shift: not given in JPL, set to 0
  Numeric psf;
  { psf = 0.0; }

  // These were all the parameters that we can extract from
  // JPL. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  data.T0 = 300.0;

  // Set line shape computer
  data.line.lineshape = LineShape::Model(sgam, nself, agam, nair, psf);

  // That's it!
  data.bad = false;
  return data;
}

bool Absorption::Lines::OK() const ARTS_NOEXCEPT {
  const Index nb = broadeningspecies.nelem();

  // Check that the isotopologue is ok
  if (not Isotopologue().OK()) return false;

  // Test that self and bath is covered by the range if set positive
  if (nb < (Index(selfbroadening) + Index(bathbroadening))) return false;

  // Test that the temperature is physical
  if (T0 <= 0) return false;

  // Test that all lines have the correct sized line shape model
  if (std::any_of(lines.cbegin(), lines.cend(), [nb](auto& line) {
        return line.LineShapeElems() != nb;
      }))
    return false;

  // Test that all lines have the correct sized local quantum numbers
  if (std::any_of(lines.cbegin(), lines.cend(), [&](auto& line) {
        return line.LocalQuantumElems() != lines.front().LocalQuantumElems();
      }))
    return false;

  // Otherwise everything is fine!
  return true;
}

Numeric Absorption::Lines::DopplerConstant(Numeric T) const noexcept {
  return std::sqrt(Constant::doppler_broadening_const_squared * T /
                   SpeciesMass());
}

namespace Absorption {
String cutofftype2metadatastring(CutoffType in, Numeric cutoff) {
  std::ostringstream os;
  switch (in) {
    case CutoffType::None:
      os << "No cut-off will be applied.\n";
      break;
    case CutoffType::ByLine:
      os << "The lines will be cut-off " << cutoff
         << " Hz from the line center + D0.\n";
      break;
    case CutoffType::FINAL:
      break;
  }
  return os.str();
}

void SingleLine::SetAutomaticZeeman(QuantumIdentifier qid) {
  for (auto& val : localquanta.val) qid.val.set(val);  // Copy to same struct

  zeeman = Zeeman::Model(qid);
}

void SingleLine::SetLineMixing2SecondOrderData(const Vector& d) {
  lineshape.SetLineMixingModel(
      LineShape::LegacyLineMixingData::vector2modellm(
          d, LineShape::LegacyLineMixingData::TypeLM::LM_2NDORDER)
          .Data()[0]);
}

void SingleLine::SetLineMixing2AER(const Vector& d) {
  const LineShape::ModelParameters Y = {
      LineShape::TemperatureModel::LM_AER, d[4], d[5], d[6], d[7]};
  const LineShape::ModelParameters G = {
      LineShape::TemperatureModel::LM_AER, d[8], d[9], d[10], d[11]};
  for (auto& sm : lineshape.Data()) {
    sm.Y() = Y;
    sm.G() = G;
  }
}

bifstream& SingleLine::read(bifstream& bif) {
  /** Standard parameters */
  bif >> F0 >> I0 >> E0 >> glow >> gupp >> A >> zeeman;

  /** Line shape model */
  lineshape.read(bif);

  /** Lower level quantum numbers */
  for (auto& val : localquanta.val) val.read(bif);

  return bif;
}

bofstream& SingleLine::write(bofstream& bof) const {
  /** Standard parameters */
  bof << F0 << I0 << E0 << glow << gupp << A << zeeman;

  /** Line shape model */
  lineshape.write(bof);

  /** Lower level quantum numbers */
  for (auto& val : localquanta.val) val.write(bof);

  return bof;
}

void Lines::AppendSingleLine(SingleLine&& sl) {
  ARTS_USER_ERROR_IF(
      NumLines() not_eq 0 and NumLocalQuanta() not_eq sl.LocalQuantumElems(),
      "Error calling appending function, bad size of quantum numbers\n"
      "Type of quantum numbers in band: ",
      lines.front().localquanta.val,
      "\nType of quantum numbers in new line:",
      sl.localquanta.val);

  ARTS_USER_ERROR_IF(
      NumLines() not_eq 0 and
          sl.LineShapeElems() not_eq lines[0].LineShapeElems(),
      "Error calling appending function, bad size of broadening species");

  lines.push_back(std::move(sl));
}

void Lines::AppendSingleLine(const SingleLine& sl) {
  AppendSingleLine(SingleLine{sl});
}

bool Lines::MatchWithExternal(const SingleLineExternal& sle,
                              const QuantumIdentifier& qid) const
    ARTS_NOEXCEPT {
  if (sle.bad) return false;
  if (sle.selfbroadening not_eq selfbroadening) return false;
  if (sle.bathbroadening not_eq bathbroadening) return false;
  if (sle.cutoff not_eq cutoff) return false;
  if (sle.mirroring not_eq mirroring) return false;
  if (sle.population not_eq population) return false;
  if (sle.normalization not_eq normalization) return false;
  if (sle.lineshapetype not_eq lineshapetype) return false;
  if (sle.T0 not_eq T0) return false;
  if (sle.cutofffreq not_eq cutofffreq) return false;
  if (sle.linemixinglimit not_eq linemixinglimit) return false;
  if (quantumidentity not_eq qid) return false;
  if (not std::equal(sle.species.cbegin(),
                     sle.species.cend(),
                     broadeningspecies.cbegin(),
                     broadeningspecies.cend()))
    return false;
  if (NumLines() not_eq 0 and
      not lines.front().localquanta.same_types_as(sle.line.localquanta))
    return false;
  if (NumLines() not_eq 0 and
      not sle.line.lineshape.Match(lines.front().lineshape).first)
    return false;
  return true;
}

std::pair<bool, bool> Lines::Match(const Lines& l) const noexcept {
  // Note: The pair here is first: matching and second: nullable
  if (l.selfbroadening not_eq selfbroadening) return {false, false};
  if (l.bathbroadening not_eq bathbroadening) return {false, false};
  if (l.cutoff not_eq cutoff) return {false, false};
  if (l.mirroring not_eq mirroring) return {false, false};
  if (l.population not_eq population) return {false, false};
  if (l.normalization not_eq normalization) return {false, false};
  if (l.lineshapetype not_eq lineshapetype) return {false, false};
  if (l.T0 not_eq T0) return {false, false};
  if (l.cutofffreq not_eq cutofffreq) return {false, false};
  if (l.linemixinglimit not_eq linemixinglimit) return {false, false};
  if (l.quantumidentity not_eq quantumidentity) return {false, false};
  if (not std::equal(l.broadeningspecies.cbegin(),
                     l.broadeningspecies.cend(),
                     broadeningspecies.cbegin(),
                     broadeningspecies.cend()))
    return {false, false};
  if (NumLines() not_eq 0 and l.NumLines() not_eq 0 and
      not lines.front().localquanta.same_types_as(l.lines.front().localquanta))
    return {false, false};
  if (NumLines() not_eq 0 and l.NumLines() not_eq 0) {
    if (auto matchpair =
            l.lines.front().lineshape.Match(lines.front().lineshape);
        not matchpair.first)
      return matchpair;
  }

  return {true, true};
}

void Lines::sort_by_frequency() {
  std::sort(
      lines.begin(), lines.end(), [](const SingleLine& a, const SingleLine& b) {
        return a.F0 < b.F0;
      });
}

void Lines::sort_by_einstein() {
  std::sort(lines.begin(),
            lines.end(),
            [](const SingleLine& a, const SingleLine& b) { return a.A < b.A; });
}

String Lines::LineShapeMetaData() const noexcept {
  return NumLines() ? LineShape::ModelShape2MetaData(lines[0].lineshape) : "";
}

Species::Species Lines::Species() const noexcept {return quantumidentity.Species();}

Species::IsotopeRecord Lines::Isotopologue() const noexcept {return quantumidentity.Isotopologue();}

Index Lines::NumLines() const noexcept {return Index(lines.size());}

Index Lines::NumBroadeners() const ARTS_NOEXCEPT {return Index(broadeningspecies.nelem());}

Index Lines::NumLocalQuanta() const noexcept {
  return lines.size() ? lines.front().localquanta.val.nelem() : 0;
}

bool Lines::OnTheFlyLineMixing() const noexcept {
  return population == PopulationType::ByMakarovFullRelmat or
         population == PopulationType::ByRovibLinearDipoleLineMixing;
}

bool Lines::DoLineMixing(Numeric P) const noexcept {
  return linemixinglimit < 0 ? true : linemixinglimit > P;
}

bool Lines::DoVmrDerivative(const QuantumIdentifier &qid) const noexcept {
  return qid.Isotopologue() == quantumidentity.Isotopologue() or
         (qid.Isotopologue().joker() and
          qid.Species() == quantumidentity.Species()) or
         std::any_of(broadeningspecies.begin(), broadeningspecies.end(),
                     [s = qid.Species()](auto &a) { return a == s; });
}

bool Lines::AnyLinemixing() const noexcept {
  for (auto &line : lines) {
    for (auto &shape : line.lineshape.Data()) {
      if (shape.Y().type not_eq LineShape::TemperatureModel::None or
          shape.G().type not_eq LineShape::TemperatureModel::None or
          shape.DV().type not_eq LineShape::TemperatureModel::None) {
        return true;
      }
    }
  }
  return false;
}

Index Lines::BroadeningSpeciesPosition(Species::Species spec) const noexcept {
  if (auto ptr =
          std::find(broadeningspecies.cbegin(), broadeningspecies.cend(), spec);
      ptr not_eq broadeningspecies.cend())
    return std::distance(broadeningspecies.cbegin(), ptr);
  return -1;
}

const Quantum::Number::Value& get(const Quantum::Number::LocalState& qns)
    ARTS_NOEXCEPT {
  ARTS_ASSERT(qns.val.has(QuantumNumberType::F) or
              qns.val.has(QuantumNumberType::J))
  return qns.val.has(QuantumNumberType::F) ? qns.val[QuantumNumberType::F]
                                           : qns.val[QuantumNumberType::J];
}

Index Lines::ZeemanCount(size_t k, Zeeman::Polarization type) const ARTS_NOEXCEPT {
  if (type == Zeeman::Polarization::None) return 1;

  // Select F before J but assume one of them exist
  auto& val = get(lines[k].localquanta);
  return Zeeman::nelem(val.upp(), val.low(), type);
}

Numeric Lines::ZeemanStrength(size_t k,
                              Zeeman::Polarization type,
                              Index i) const ARTS_NOEXCEPT {
  if (type == Zeeman::Polarization::None) return 1.0;

  // Select F before J but assume one of them exist
  auto& val = get(lines[k].localquanta);
  return lines[k].zeeman.Strength(val.upp(), val.low(), type, i);
}

Numeric Lines::ZeemanSplitting(size_t k,
                               Zeeman::Polarization type,
                               Index i) const ARTS_NOEXCEPT {
  if (type == Zeeman::Polarization::None) return 0.0;

  // Select F before J but assume one of them exist
  auto& val = get(lines[k].localquanta);
  return lines[k].zeeman.Splitting(val.upp(), val.low(), type, i);
}

void Lines::SetAutomaticZeeman() noexcept {
  for (auto& line : lines) line.SetAutomaticZeeman(quantumidentity);
}

Numeric Lines::F_mean(const ConstVectorView& wgts) const noexcept {
  const Numeric val =
      std::inner_product(lines.cbegin(),
                         lines.cend(),
                         wgts.begin(),
                         0.0,
                         std::plus<>(),
                         [](const auto& a, const auto& b) { return a.F0 * b; });
  const Numeric div = sum(wgts);
  return val / div;
}

Numeric Lines::F_mean(Numeric T) const noexcept {
  if (T <= 0) T = T0;

  const Index n = NumLines();
  const Numeric QT = single_partition_function(T, Isotopologue());
  const Numeric QT0 = single_partition_function(T0, Isotopologue());
  const Numeric ratiopart = QT0 / QT;

  Vector wgts(n);
  for (Index i = 0; i < n; i++) {
    const Numeric pop0 =
        (lines[i].gupp / QT0) * boltzman_factor(T0, lines[i].E0);
    const Numeric pop = pop0 * ratiopart * boltzman_ratio(T, T0, lines[i].E0);
    const Numeric dip_squared =
        -lines[i].I0 /
        (pop0 * lines[i].F0 *
         std::expm1(-(Constant::h * lines[i].F0) / (Constant::k * T0)));
    wgts[i] = pop * dip_squared;
  }

  return F_mean(wgts);
}

Numeric Lines::CutoffFreq(size_t k, Numeric shift) const noexcept {
  switch (cutoff) {
    case CutoffType::ByLine:
      return lines[k].F0 + cutofffreq +
             (mirroring == MirroringType::Manual ? -shift : shift);
    case CutoffType::None:
      return std::numeric_limits<Numeric>::max();
    case CutoffType::FINAL:
      break;
  }
  return std::numeric_limits<Numeric>::max();
}

Numeric Lines::CutoffFreqMinus(size_t k, Numeric shift) const noexcept {
  switch (cutoff) {
    case CutoffType::ByLine:
      return lines[k].F0 - cutofffreq +
             (mirroring == MirroringType::Manual ? -shift : shift);
    case CutoffType::None:
      return std::numeric_limits<Numeric>::lowest();
    case CutoffType::FINAL:
      break;
  }
  
  return std::numeric_limits<Numeric>::lowest();
}

bifstream& Lines::read(bifstream& is) {
  for (auto& line : lines) line.read(is);
  return is;
}

bofstream& Lines::write(bofstream& os) const {
  for (auto& line : lines) line.write(os);
  return os;
}

QuantumIdentifier Lines::QuantumIdentityOfLine(Index k) const noexcept {
  QuantumIdentifier qid = quantumidentity;
  for (auto& a : lines[k].localquanta.val) qid.val.set(a);
  return qid;
}

void Lines::MakeLineShapeModelCommon() {
  const Index n = NumLines();
  const Index m = NumBroadeners();
  if (not n) return;

  for (Index j = 0; j < m; j++) {
    for (Index i = 0; i < LineShape::nVars; i++) {
      // Find a common type (same or none) or throw if there is none
      LineShape::TemperatureModel t = LineShape::TemperatureModel::None;
      for (auto& line : lines) {
        if (auto& data = line.lineshape[j].Data()[i];
            not LineShape::modelparameterEmpty(data)) {
          if (t == LineShape::TemperatureModel::None) t = data.type;
          ARTS_USER_ERROR_IF(
              t not_eq data.type,
              "Cannot make a common line shape model for the band as there are multiple non-empty types: ",
              data.type,
              " and ",
              t)
        }
      }

      // Set the common type
      for (auto& line : lines) {
        if (auto& data = line.lineshape[j].Data()[i]; data.type not_eq t) {
          data = LineShape::modelparameterGetEmpty(t);
        }
      }
    }
  }
}
}  // namespace Absorption

AbsorptionPopulationTagTypeStatus::AbsorptionPopulationTagTypeStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      switch (band.population) {
        case Absorption::PopulationType::LTE:
          LTE = true;
          break;
        case Absorption::PopulationType::NLTE:
          NLTE = true;
          break;
        case Absorption::PopulationType::VibTemps:
          VibTemps = true;
          break;
        case Absorption::PopulationType::ByHITRANFullRelmat:
          ByHITRANFullRelmat = true;
          break;
        case Absorption::PopulationType::ByHITRANRosenkranzRelmat:
          ByHITRANRosenkranzRelmat = true;
          break;
        case Absorption::PopulationType::ByMakarovFullRelmat:
          ByMakarovFullRelmat = true;
          break;
        case Absorption::PopulationType::ByRovibLinearDipoleLineMixing:
          ByRovibLinearDipoleLineMixing = true;
          break;
        case Absorption::PopulationType::FINAL: { /* leave list */
        }
      }
    }
  }
}

AbsorptionCutoffTagTypeStatus::AbsorptionCutoffTagTypeStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      switch (band.cutoff) {
        case Absorption::CutoffType::None:
          None = true;
          break;
        case Absorption::CutoffType::ByLine:
          ByLine = true;
          break;
        case Absorption::CutoffType::FINAL: { /* leave list */
        }
      }
    }
  }
}

AbsorptionLineShapeTagTypeStatus::AbsorptionLineShapeTagTypeStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      switch (band.lineshapetype) {
        case LineShape::Type::DP:
          DP = true;
          break;
        case LineShape::Type::LP:
          LP = true;
          break;
        case LineShape::Type::VP:
          VP = true;
          break;
        case LineShape::Type::SDVP:
          SDVP = true;
          break;
        case LineShape::Type::HTP:
          HTP = true;
          break;
        case LineShape::Type::SplitLP:
          SplitLP = true;
          break;
        case LineShape::Type::SplitVP:
          SplitVP = true;
          break;
        case LineShape::Type::SplitSDVP:
          SplitSDVP = true;
          break;
        case LineShape::Type::SplitHTP:
          SplitHTP = true;
          break;
        case LineShape::Type::FINAL: { /* leave list */
        }
      }
    }
  }
}

AbsorptionMirroringTagTypeStatus::AbsorptionMirroringTagTypeStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      switch (band.mirroring) {
        case Absorption::MirroringType::None:
          None = true;
          break;
        case Absorption::MirroringType::Lorentz:
          Lorentz = true;
          break;
        case Absorption::MirroringType::SameAsLineShape:
          SameAsLineShape = true;
          break;
        case Absorption::MirroringType::Manual:
          Manual = true;
          break;
        case Absorption::MirroringType::FINAL: { /* leave list */
        }
      }
    }
  }
}

AbsorptionNormalizationTagTypeStatus::AbsorptionNormalizationTagTypeStatus(
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      switch (band.normalization) {
        case Absorption::NormalizationType::None:
          None = true;
          break;
        case Absorption::NormalizationType::VVH:
          VVH = true;
          break;
        case Absorption::NormalizationType::VVW:
          VVW = true;
          break;
        case Absorption::NormalizationType::RQ:
          RQ = true;
          break;
        case Absorption::NormalizationType::SFS:
          SFS = true;
          break;
        case Absorption::NormalizationType::FINAL: { /* leave list */
        }
      }
    }
  }
}

namespace {
template <typename T>
constexpr std::size_t req_spaces(T my_enum) {
  constexpr std::size_t n = []() {
    std::size_t longest = 0;
    for (auto& x : Absorption::enumstrs::CutoffTypeNames) {
      longest = std::max(x.length(), longest);
    }
    for (auto& x : Absorption::enumstrs::MirroringTypeNames) {
      longest = std::max(x.length(), longest);
    }
    for (auto& x : Absorption::enumstrs::NormalizationTypeNames) {
      longest = std::max(x.length(), longest);
    }
    for (auto& x : Absorption::enumstrs::PopulationTypeNames) {
      longest = std::max(x.length(), longest);
    }
    for (auto& x : LineShape::enumstrs::TypeNames) {
      longest = std::max(x.length(), longest);
    }
    return longest + 1;
  }();
  return n - toString(my_enum).length();
}

auto spaces(std::size_t n) { return std::basic_string(n, ' '); }
}  // namespace

std::ostream& operator<<(std::ostream& os, AbsorptionCutoffTagTypeStatus val) {
  // Trick to never forget to update
  Absorption::CutoffType x{Absorption::CutoffType::FINAL};
  switch (x) {
    case Absorption::CutoffType::FINAL:
      os << "Cutoff tag types:\n";
      [[fallthrough]];
    case AbsorptionCutoffType::ByLine:
      os << "    ByLine:" << spaces(req_spaces(AbsorptionCutoffType::ByLine))
         << val.ByLine << '\n';
      [[fallthrough]];
    case AbsorptionCutoffType::None:
      os << "    None:" << spaces(req_spaces(AbsorptionCutoffType::None))
         << val.None;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         AbsorptionMirroringTagTypeStatus val) {
  // Trick to never forget to update
  Absorption::MirroringType x{Absorption::MirroringType::FINAL};
  switch (x) {
    case Absorption::MirroringType::FINAL:
      os << "Mirroring tag types:\n";
      [[fallthrough]];
    case Absorption::MirroringType::None:
      os << "    None:" << spaces(req_spaces(Absorption::MirroringType::None))
         << val.None << '\n';
      [[fallthrough]];
    case Absorption::MirroringType::Lorentz:
      os << "    Lorentz:"
         << spaces(req_spaces(Absorption::MirroringType::Lorentz))
         << val.Lorentz << '\n';
      [[fallthrough]];
    case Absorption::MirroringType::SameAsLineShape:
      os << "    SameAsLineShape:"
         << spaces(req_spaces(Absorption::MirroringType::SameAsLineShape))
         << val.SameAsLineShape << '\n';
      [[fallthrough]];
    case Absorption::MirroringType::Manual:
      os << "    Manual:"
         << spaces(req_spaces(Absorption::MirroringType::Manual)) << val.Manual;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         AbsorptionNormalizationTagTypeStatus val) {
  // Trick to never forget to update
  Absorption::NormalizationType x{Absorption::NormalizationType::FINAL};
  switch (x) {
    case Absorption::NormalizationType::FINAL:
      os << "Normalization tag types:\n";
      [[fallthrough]];
    case Absorption::NormalizationType::None:
      os << "    None:"
         << spaces(req_spaces(Absorption::NormalizationType::None)) << val.None
         << '\n';
      [[fallthrough]];
    case Absorption::NormalizationType::VVH:
      os << "    VVH:" << spaces(req_spaces(Absorption::NormalizationType::VVH))
         << val.VVH << '\n';
      [[fallthrough]];
    case Absorption::NormalizationType::VVW:
      os << "    VVW:" << spaces(req_spaces(Absorption::NormalizationType::VVW))
         << val.VVW << '\n';
      [[fallthrough]];
    case Absorption::NormalizationType::RQ:
      os << "    RQ:" << spaces(req_spaces(Absorption::NormalizationType::RQ))
         << val.RQ << '\n';
      [[fallthrough]];
    case Absorption::NormalizationType::SFS:
      os << "    SFS:" << spaces(req_spaces(Absorption::NormalizationType::SFS))
         << val.SFS;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         AbsorptionPopulationTagTypeStatus val) {
  // Trick to never forget to update
  Absorption::PopulationType x{Absorption::PopulationType::FINAL};
  switch (x) {
    case Absorption::PopulationType::FINAL:
      os << "Population tag type:\n";
      [[fallthrough]];
    case Absorption::PopulationType::LTE:
      os << "    LTE:" << spaces(req_spaces(Absorption::PopulationType::LTE))
         << val.LTE << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::NLTE:
      os << "    NLTE:" << spaces(req_spaces(Absorption::PopulationType::NLTE))
         << val.NLTE << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::VibTemps:
      os << "    VibTemps:"
         << spaces(req_spaces(Absorption::PopulationType::VibTemps))
         << val.VibTemps << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::ByHITRANRosenkranzRelmat:
      os << "    ByHITRANRosenkranzRelmat:"
         << spaces(req_spaces(
                Absorption::PopulationType::ByHITRANRosenkranzRelmat))
         << val.ByHITRANRosenkranzRelmat << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::ByHITRANFullRelmat:
      os << "    ByHITRANFullRelmat:"
         << spaces(req_spaces(Absorption::PopulationType::ByHITRANFullRelmat))
         << val.ByHITRANFullRelmat << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::ByMakarovFullRelmat:
      os << "    ByMakarovFullRelmat:"
         << spaces(req_spaces(Absorption::PopulationType::ByMakarovFullRelmat))
         << val.ByMakarovFullRelmat << '\n';
      [[fallthrough]];
    case Absorption::PopulationType::ByRovibLinearDipoleLineMixing:
      os << "    ByRovibLinearDipoleLineMixing:"
         << spaces(req_spaces(
                Absorption::PopulationType::ByRovibLinearDipoleLineMixing))
         << val.ByRovibLinearDipoleLineMixing;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         AbsorptionLineShapeTagTypeStatus val) {
  // Trick to never forget to update
  LineShapeType x{LineShapeType::FINAL};
  switch (x) {
    case LineShapeType::FINAL:
      os << "Line shape tag type:\n";
      [[fallthrough]];
    case LineShapeType::DP:
      os << "    DP:" << spaces(req_spaces(LineShapeType::DP)) << val.DP
         << '\n';
      [[fallthrough]];
    case LineShapeType::LP:
      os << "    LP:" << spaces(req_spaces(LineShapeType::LP)) << val.LP
         << '\n';
      [[fallthrough]];
    case LineShapeType::VP:
      os << "    VP:" << spaces(req_spaces(LineShapeType::VP)) << val.VP
         << '\n';
      [[fallthrough]];
    case LineShapeType::SDVP:
      os << "    SDVP:" << spaces(req_spaces(LineShapeType::SDVP)) << val.SDVP
         << '\n';
      [[fallthrough]];
    case LineShapeType::HTP:
      os << "    HTP:" << spaces(req_spaces(LineShapeType::HTP)) << val.HTP
         << '\n';
      [[fallthrough]];
    case LineShapeType::SplitLP:
      os << "    SplitLP:" << spaces(req_spaces(LineShapeType::SplitLP)) << val.SplitLP
         << '\n';
      [[fallthrough]];
    case LineShapeType::SplitVP:
      os << "    SplitVP:" << spaces(req_spaces(LineShapeType::SplitVP)) << val.SplitVP
         << '\n';
      [[fallthrough]];
    case LineShapeType::SplitSDVP:
      os << "    SplitSDVP:" << spaces(req_spaces(LineShapeType::SplitSDVP)) << val.SplitSDVP
         << '\n';
      [[fallthrough]];
    case LineShapeType::SplitHTP:
      os << "    SplitHTP:" << spaces(req_spaces(LineShapeType::SplitHTP)) << val.SplitHTP;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, AbsorptionTagTypesStatus val) {
  return os << "Catalog tag summary:\n  " << val.cutoff << "\n  "
            << val.lineshapetype << "\n  " << val.mirroring << "\n  "
            << val.normalization << "\n  " << val.population;
}

AbsorptionSpeciesBandIndex flat_index(
    Index i,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  const Index n = abs_species.nelem();

  Index ispec{0};
  while (ispec < n and i >= abs_lines_per_species[ispec].nelem()) {
    i -= abs_lines_per_species[ispec].nelem();
    ++ispec;
  }

  return {ispec, i};
}

namespace Absorption {
bool any_cutoff(const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  return std::any_of(
      abs_lines_per_species.cbegin(),
      abs_lines_per_species.cend(),
      [](auto& abs_lines) {
        return std::any_of(
            abs_lines.cbegin(), abs_lines.cend(), [](auto& band) {
              return band.cutoff not_eq CutoffType::None;
            });
      });
}

Rational Lines::max(QuantumNumberType x) const {
  ARTS_USER_ERROR_IF(Quantum::Number::common_value_type(x) == Quantum::Number::ValueType::S, "Cannot get a rational out from quantum number type ", x)
  if (quantumidentity.val.has(x)) {
    auto& val = quantumidentity.val[x];
    return std::max(val.low(), val.upp());
  }

  Rational out{std::numeric_limits<Index>::lowest()};
  for (auto& line : lines) {
    ARTS_USER_ERROR_IF(not line.localquanta.val.has(x), "No ", x, " in some line(s)")
    auto& val = line.localquanta.val[x];
    out = std::max(std::max(val.low(), val.upp()), out);
  }
  return out;
}

/** Number of lines */
Index nelem(const Lines &l) { return l.NumLines(); }

/** Number of lines in list */
Index nelem(const Array<Lines> &l) {
  Index n = 0;
  for (auto &x : l)
    n += nelem(x);
  return n;
}

/** Number of lines in lists */
Index nelem(const Array<Array<Lines>> &l) {
  Index n = 0;
  for (auto &x : l)
    n += nelem(x);
  return n;
}
} // namespace Absorption
