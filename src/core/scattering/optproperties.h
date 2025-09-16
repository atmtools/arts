/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   optproperties.h
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   2003-03-06
  
  \brief  Scattering database structure and functions.
  
   This file contains the definition of the SingleScatteringData structure
   and the functions in optproperties.cc that are of interest elsewhere.
*/

#ifndef optproperties_h
#define optproperties_h

#include <matpack.h>
#include <mystring.h>
#include <rtepack.h>

//! An attribute to classify the particle type (ptype) of a SingleScatteringData
//structure (a scattering element).
/*! 
  PTYPE_GENERAL      General case
  PTYPE_TOTAL_RND    Totally randomly oriented particles
  PTYPE_AZIMUTH_RND  Azimuthally randomly oriented particles

  A detailed description of the different cases can be found in AUG.

*/
enum PType : Index {
  PTYPE_GENERAL     = 300,
  PTYPE_AZIMUTH_RND = 200,
  PTYPE_TOTAL_RND   = 100,
};

//! An attribute to classify the method to be used for SingleScatteringData
/*!
  NONE         Dummy value
  TMATRIX      T-Matrix method
*/
enum ParticleSSDMethod {
  PARTICLE_SSDMETHOD_NONE    = 0,
  PARTICLE_SSDMETHOD_TMATRIX = 1,
};

template <>
struct std::formatter<ParticleSSDMethod> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ParticleSSDMethod& v,
                              FmtContext& ctx) const {
    const std::string_view quote = tags.quote();

    switch (v) {
      case PARTICLE_SSDMETHOD_NONE:
        tags.format(ctx, quote, "PARTICLE_SSDMETHOD_NONE"sv, quote);
        break;
      case PARTICLE_SSDMETHOD_TMATRIX:
        tags.format(ctx, quote, "PARTICLE_SSDMETHOD_TMATRIX"sv, quote);
        break;
      default:
        tags.format(ctx, quote, "PARTICLE_SSDMETHOD_UNKNOWN"sv, quote);
        break;
    }

    return ctx.out();
  }
};

/*===========================================================================
  === The SingleScatteringData structure
  ===========================================================================*/

/*!
   Structure which describes the single scattering properties of a scattering
   element (a single particle or a particle bulk).

   The fields of the structure are described in the ARTS user guide (AUG).
   It is listed as a sub-entry to "data structures".  
*/
struct SingleScatteringData {
  PType ptype;
  String description;
  Vector f_grid;
  Vector T_grid;
  Vector za_grid;
  Vector aa_grid;
  Tensor7 pha_mat_data;
  Tensor5 ext_mat_data;
  Tensor5 abs_vec_data;
};

typedef Array<SingleScatteringData> ArrayOfSingleScatteringData;
typedef Array<Array<SingleScatteringData> > ArrayOfArrayOfSingleScatteringData;

/*===========================================================================
  === The ScatteringMetaData structure
  ===========================================================================*/
/*!
   Structure which holds the meta data of a scattering element (a single
   particle or an ensemble of particles), mainly microphysical parameter
   necessary for calculating size and shape distributions.

  For a description of the structure members see built-in documentation of
  scat_meta_single in workspace.cc 
*/
struct ScatteringMetaData {
  String description;
  String source;
  String refr_index;
  Numeric mass;
  Numeric diameter_max;
  Numeric diameter_volume_equ;
  Numeric diameter_area_equ_aerodynamical;
};

typedef Array<ScatteringMetaData> ArrayOfScatteringMetaData;
typedef Array<Array<ScatteringMetaData> > ArrayOfArrayOfScatteringMetaData;

template <>
struct std::formatter<ScatteringMetaData> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ScatteringMetaData& v,
                              FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');

    const std::string_view sep   = tags.sep();
    const std::string_view quote = tags.quote();

    tags.format(ctx,
                quote,
                v.description,
                quote,
                sep,
                quote,
                v.source,
                quote,
                sep,
                quote,
                v.refr_index,
                quote,
                sep,
                v.mass,
                sep,
                v.diameter_max,
                sep,
                v.diameter_volume_equ,
                sep,
                v.diameter_area_equ_aerodynamical);

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<PType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const PType& v, FmtContext& ctx) const {
    const std::string_view quote = tags.quote();

    switch (v) {
      case PTYPE_GENERAL:
        tags.format(ctx, quote, "PTYPE_GENERAL"sv, quote);
        break;
      case PTYPE_AZIMUTH_RND:
        tags.format(ctx, quote, "PTYPE_AZIMUTH_RND"sv, quote);
        break;
      case PTYPE_TOTAL_RND:
        tags.format(ctx, quote, "PTYPE_TOTAL_RND"sv, quote);
        break;
      default: tags.format(ctx, quote, "PTYPE_UNKNOWN"sv, quote); break;
    }

    return ctx.out();
  }
};

template <>
struct std::formatter<SingleScatteringData> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  [[nodiscard]] std::string to_string(const SingleScatteringData& v) const;

  template <class FmtContext>
  FmtContext::iterator format(const SingleScatteringData& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, to_string(v));
  }
};

std::string_view PTypeToString(const PType ptype);

PType PType2FromString(const std::string_view ptype_string);

PType PTypeFromString(const std::string_view ptype_string);

template <>
struct xml_io_stream<SingleScatteringData> {
  static constexpr std::string_view type_name = "SingleScatteringData"sv;

  static void write(std::ostream& os,
                    const SingleScatteringData& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SingleScatteringData& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<ScatteringMetaData> {
  static constexpr std::string_view type_name = "ScatteringMetaData"sv;

  static void write(std::ostream& os,
                    const ScatteringMetaData& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   ScatteringMetaData& x,
                   bifstream* pbifs = nullptr);
};

#endif  //optproperties_h
