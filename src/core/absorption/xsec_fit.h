/*!
  \file   xsec_fit.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/

#ifndef HITRAN_XSEC_H
#define HITRAN_XSEC_H

#include <array.h>
#include <matpack.h>
#include <mystring.h>
#include <species.h>
#include <xml.h>

#include <memory>

/** Hitran crosssection class.
 *
 * Stores the coefficients from our model for hitran crosssection data and
 * applies them to calculate the crossections.
 */
class XsecRecord {
 public:
  /** Return species index */
  [[nodiscard]] const SpeciesEnum& Species() const;

  /** Set species name */
  void SetSpecies(const SpeciesEnum species);

  /** Return species name */
  [[nodiscard]] String SpeciesName() const;

  /** Return species index */
  [[nodiscard]] static constexpr Index Version() { return mversion; };

  /** Set species name */
  void SetVersion(Index version);

  /** Calculate hitran cross section data.

     Calculate crosssections at each frequency for given pressure and
     temperature.

     \param[out] result     Crosssections for given frequency grid.
     \param[in] f_grid      Frequency grid.
     \param[in] pressure    Scalar pressure.
     \param[in] temperature Scalar temperature.
     */
  void Extract(VectorView result,
               const Vector& f_grid,
               Numeric pressure,
               Numeric temperature) const;

  /************ VERSION 2 *************/
  /** Get mininum pressures from fit */
  [[nodiscard]] const Vector& FitMinPressures() const;

  /** Get maximum pressures from fit */
  [[nodiscard]] const Vector& FitMaxPressures() const;

  /** Get mininum temperatures from fit */
  [[nodiscard]] const Vector& FitMinTemperatures() const;

  /** Get maximum temperatures */
  [[nodiscard]] const Vector& FitMaxTemperatures() const;

  /** Get coefficients */
  [[nodiscard]] const ArrayOfGriddedField1Named& FitCoeffs() const;

  /** Get mininum pressures from fit */
  [[nodiscard]] Vector& FitMinPressures();

  /** Get maximum pressures from fit */
  [[nodiscard]] Vector& FitMaxPressures();

  /** Get mininum temperatures from fit */
  [[nodiscard]] Vector& FitMinTemperatures();

  /** Get maximum temperatures */
  [[nodiscard]] Vector& FitMaxTemperatures();

  /** Get coefficients */
  [[nodiscard]] ArrayOfGriddedField1Named& FitCoeffs();

  friend std::ostream& operator<<(std::ostream& os, const XsecRecord& xd);

 private:
  /** Calculate crosssections */
  void CalcXsec(VectorView xsec,
                const Index dataset,
                const Numeric pressure,
                const Numeric temperature) const;

  // /** Calculate temperature derivative of crosssections */
  // void CalcDT(VectorView xsec_dt,
  //             Index dataset,
  //             Numeric temperature) const;

  // /** Calculate pressure derivative of crosssections */
  // void CalcDP(VectorView xsec_dp,
  //             Index dataset,
  //             Numeric pressure) const;

 public:
  static constexpr Index P00 = 0;
  static constexpr Index P10 = 1;
  static constexpr Index P01 = 2;
  static constexpr Index P20 = 3;

  static constexpr Index mversion{2};
  SpeciesEnum mspecies;
  /* VERSION 2 */
  Vector mfitminpressures;
  Vector mfitmaxpressures;
  Vector mfitmintemperatures;
  Vector mfitmaxtemperatures;
  ArrayOfGriddedField1Named mfitcoeffs;
};

using ArrayOfXsecRecord = Array<XsecRecord>;

std::ostream& operator<<(std::ostream& os, const ArrayOfXsecRecord& x);

Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            SpeciesEnum species);

XsecRecord* hitran_xsec_get_data(
    const std::shared_ptr<std::vector<XsecRecord>>& xsec_data,
    const SpeciesEnum species);

template <>
struct std::formatter<XsecRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const XsecRecord& v, FmtContext& ctx) const {
    if (tags.short_str) {
      return tags.format(ctx, "XsecRecord("sv, v.Species(), ")"sv);
    }

    const std::string_view sep = tags.sep(true);
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.Species(),
                sep,
                v.FitMinPressures(),
                sep,
                v.FitMaxPressures(),
                sep,
                v.FitMinTemperatures(),
                sep,
                v.FitMaxTemperatures(),
                sep,
                v.FitCoeffs());
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct xml_io_stream<XsecRecord> {
  static constexpr std::string_view type_name = "XsecRecord"sv;

  static void write(std::ostream& os,
                    const XsecRecord& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, XsecRecord& x, bifstream* pbifs = nullptr);
};

#endif  // HITRAN_XSEC_H
