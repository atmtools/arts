#pragma once

#include <atm.h>
#include <lbl.h>
#include <lbl_data.h>
#include <rtepack.h>
#include <xml.h>

#include <memory>

#include "fwd_cia.h"
#include "fwd_hxsec.h"
#include "fwd_predef.h"

namespace fwd {
struct propmat {
  std::shared_ptr<AtmPoint> atm{};

  lbl::fwd::line_storage lines{};

  cia::full cia{};
  predef::full predef{};
  hxsec::full xsec{};

  propmat();
  propmat(const propmat&);
  propmat(propmat&&) noexcept;
  propmat& operator=(const propmat&);
  propmat& operator=(propmat&&) noexcept;

  propmat(std::shared_ptr<AtmPoint> atm,
          std::shared_ptr<AbsorptionBands> lines,
          std::shared_ptr<CIARecords> cia,
          std::shared_ptr<XsecRecords> xsec,
          std::shared_ptr<PredefinedModelData> predef,
          Numeric ciaextrap = {},
          Index ciarobust   = {});

  std::pair<Propmat, Stokvec> operator()(const Numeric frequency,
                                         const Vector2 los) const;

  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_ciaextrap(Numeric extrap);
  void set_ciarobust(Index robust);
  void set_bands(std::shared_ptr<AbsorptionBands> lines);
  void set_cia(std::shared_ptr<CIARecords> cia);
  void set_predef(std::shared_ptr<PredefinedModelData> predef);
  void set_model(std::shared_ptr<XsecRecords> xsec);
};  // struct propmat
}  // namespace fwd

template <>
struct xml_io_stream<fwd::propmat> {
  static constexpr std::string_view type_name = "ForwardPropmat"sv;

  static void write(std::ostream& os,
                    const fwd::propmat& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   fwd::propmat& x,
                   bifstream* pbifs = nullptr);
};
