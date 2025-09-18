#include "xml_io_old.h"

#include <double_imanip.h>
#include <string_extract.h>
#include <xml_io.h>

#include <algorithm>
#include <istream>
#include <ostream>

namespace {
ArtscatMeta ReadFromArtscat3Stream(std::istream& is) {
  // Default data and values for this type
  ArtscatMeta output{};
  output.data.ls.one_by_one = false;
  output.data.ls.single_models.reserve(2);

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return output;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.size() == 0 && is.eof()) return output;

    // @ as first character marks catalogue entry
    char c;
    extract(c, line, 1);

    // check for empty line
    if (c == '@') {
      comment = false;
    }
  }

  // read the arts identifier String
  std::istringstream icecream(line);

  String artsid;
  icecream >> artsid;

  try {
    if (artsid.length() != 0) {
      // Set the species
      const auto isotopologue =
          SpeciesIsotope(Species::update_isot_name(artsid));
      ARTS_USER_ERROR_IF(
          isotopologue.is_joker(),
          "A line catalog species can only be of the form \"Plain\", meaning it\nhas the form SPECIES-ISONUM.\n"
          "Your input contains: ",
          artsid,
          ". which we cannot interpret as a plain species")
      output.quantumidentity = QuantumIdentifier{isotopologue};

      // Extract center frequency:
      icecream >> double_imanip() >> output.data.f0;

      Numeric psf;
      // Extract pressure shift:
      icecream >> double_imanip() >> psf;

      Numeric I0;
      // Extract intensity:
      icecream >> double_imanip() >> I0;

      // Extract reference temperature for Intensity in K:
      icecream >> double_imanip() >> output.data.ls.T0;

      // Extract lower state energy:
      icecream >> double_imanip() >> output.data.e0;

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
      if (tgam != output.data.ls.T0) {
        agam = agam * pow(tgam / output.data.ls.T0, nair);
        sgam = sgam * pow(tgam / output.data.ls.T0, nself);
        psf  = psf *
              pow(tgam / output.data.ls.T0, (Numeric).25 + (Numeric)1.5 * nair);
      }

      // Set line shape computer
      output.data.ls.single_models[isotopologue.spec]
          .data[LineShapeModelVariable::G0] =
          lbl::temperature::data{LineShapeModelType::T1, Vector{sgam, nair}};
      output.data.ls.single_models[isotopologue.spec]
          .data[LineShapeModelVariable::D0] =
          lbl::temperature::data{LineShapeModelType::T5, Vector{psf, nair}};
      output.data.ls.single_models[SpeciesEnum::Bath]
          .data[LineShapeModelVariable::G0] =
          lbl::temperature::data{LineShapeModelType::T1, Vector{agam, nair}};
      output.data.ls.single_models[SpeciesEnum::Bath]
          .data[LineShapeModelVariable::D0] =
          lbl::temperature::data{LineShapeModelType::T5, Vector{psf, nair}};

      if (not(output.data.gu > 0.0)) output.data.gu = 1.;
      if (not(output.data.gl > 0.0)) output.data.gl = 1.;
      output.data.a =
          output.data.compute_a(I0, isotopologue, output.data.ls.T0);
    }
  } catch (std::exception& e) {
    throw std::runtime_error(std::format(R"(Cannot parse line: 

{}

The error is:

{}
)",
                                         line,
                                         e.what()));
  }

  // That's it!
  output.bad = false;
  return output;
}

lbl::line_shape::model from_artscat4(std::istream& is,
                                     const Numeric T0_,
                                     const QuantumIdentifier& qid) {
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;
  using enum LineShapeModelCoefficient;

  lbl::line_shape::model out;
  out.one_by_one = false;
  out.T0         = T0_;
  out.single_models.reserve(7);

  const std::array<SpeciesEnum, 7> species = {qid.isot.spec,
                                              to<SpeciesEnum>("N2"),
                                              to<SpeciesEnum>("O2"),
                                              to<SpeciesEnum>("H2O"),
                                              to<SpeciesEnum>("CO2"),
                                              to<SpeciesEnum>("H2"),
                                              to<SpeciesEnum>("He")};

  for (auto& spec : species) {
    out.single_models[spec].data[G0] = lbl::temperature::data{T1, Vector{0, 0}};
    out.single_models[spec].data[D0] = lbl::temperature::data{T5, Vector{0, 0}};
  }

  // G0 main coefficient
  for (auto& spec : species) {
    is >> double_imanip() >> out.single_models[spec].data[G0].X(X0);
  };

  // G0 exponent is same as D0 exponent
  for (auto& spec : species) {
    is >> double_imanip() >> out.single_models[spec].data[G0].X(X1);
    out.single_models[spec].data[D0].X(X1) =
        out.single_models[spec].data[G0].X(X1);
  };

  // D0 coefficient
  out.single_models[species[0]].data[D0].X(X0) = 0;
  for (auto& spec : species | stdv::drop(1)) {
    is >> double_imanip() >> out.single_models[spec].data[D0].X(X0);
  }

  return out;
}

ArtscatMeta ReadFromArtscat4Stream(std::istream& is) {
  // Default data and values for this type
  ArtscatMeta output{};

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return output;

    // Throw runtime_error if stream is bad:
    ARTS_USER_ERROR_IF(!is, "Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.size() == 0 && is.eof()) return output;

    // @ as first character marks catalogue entry
    char c;
    extract(c, line, 1);

    // check for empty line
    if (c == '@') {
      comment = false;
    }
  }

  // read the arts identifier String
  std::istringstream icecream(line);

  String artsid;
  icecream >> artsid;

  try {
    if (artsid.length() != 0) {
      // Set line ID
      const auto isotopologue =
          SpeciesIsotope(Species::update_isot_name(artsid));
      ARTS_USER_ERROR_IF(
          isotopologue.is_joker(),
          "A line catalog species can only be of the form \"Plain\", meaning it\nhas the form SPECIES-ISONUM.\n"
          "Your input contains: ",
          artsid,
          ". which we cannot interpret as a plain species")
      output.quantumidentity = QuantumIdentifier{isotopologue};

      // Extract center frequency:
      icecream >> double_imanip() >> output.data.f0;

      // Extract intensity:
      Numeric I0;
      icecream >> double_imanip() >> I0;

      // Extract reference temperature for Intensity in K:
      Numeric T0;
      icecream >> double_imanip() >> T0;

      // Extract lower state energy:
      icecream >> double_imanip() >> output.data.e0;

      // Extract Einstein A-coefficient:
      icecream >> double_imanip() >> output.data.a;

      // Extract upper state stat. weight:_
      icecream >> double_imanip() >> output.data.gu;

      // Extract lower state stat. weight:
      icecream >> double_imanip() >> output.data.gl;

      output.data.ls = from_artscat4(icecream, T0, output.quantumidentity);

      if (not(output.data.gu > 0.0)) output.data.gu = 1.;
      if (not(output.data.gl > 0.0)) output.data.gl = 1.;

      output.data.a = output.data.compute_a(I0, isotopologue, T0);
    }
  } catch (std::exception& e) {
    throw std::runtime_error(std::format(R"(Cannot parse line:

{}

The error is:

{}
)",
                                         line,
                                         e.what()));
  }

  // That's it!
  output.bad = false;
  return output;
}
}  // namespace

void xml_io_stream<ArrayOfArtscatMeta>::write(std::ostream&,
                                              const ArrayOfArtscatMeta&,
                                              bofstream*,
                                              std::string_view) {
  assert(false);
}

void xml_io_stream<ArrayOfArtscatMeta>::read(std::istream& is_xml,
                                             ArrayOfArtscatMeta& meta,
                                             bifstream*) {
  meta.clear();

  XMLTag tag;

  String version;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("ArrayOfLineRecord");

  tag.get_attribute_value("version", version);
  tag.get_attribute_value("nelem", nelem);

  if (version == "ARTSCAT-3") {
    for (Index i = 0; i < nelem; ++i) {
      meta.emplace_back() = ReadFromArtscat3Stream(is_xml);
      ARTS_USER_ERROR_IF(
          meta.back().bad,
          "Error reading line data from ARTSCAT-3.  Check the input file.");
    }
  } else if (version == "ARTSCAT-4") {
    for (Index i = 0; i < nelem; ++i) {
      meta.emplace_back() = ReadFromArtscat4Stream(is_xml);
      ARTS_USER_ERROR_IF(
          meta.back().bad,
          "Error reading line data from ARTSCAT-4.  Check the input file.");
    }
  } else {
    ARTS_USER_ERROR(
        "Unknown version of ARTSCAT: {}, supported versions are ARTSCAT-3 and ARTSCAT-4",
        version);
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/ArrayOfLineRecord");
}
