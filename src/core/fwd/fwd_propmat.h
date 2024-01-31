#pragma once

#include <lbl.h>

#include <memory>

#include "atm.h"
#include "fwd_cia.h"
#include "fwd_hxsec.h"
#include "fwd_predef.h"
#include "lbl_data.h"
#include "rtepack.h"

namespace fwd {
class propmat_operator {
  lbl::fwd::line_storage lines{};
  cia::full cia{};
  predef::full predef{};
  hxsec::full hxsec{};

  std::shared_ptr<AtmPoint> atm{};
  lbl::zeeman::pol pol{lbl::zeeman::pol::no};

  void adapt();

 public:
  propmat_operator() = default;
  propmat_operator(const propmat_operator&) = default;
  propmat_operator(propmat_operator&&) = default;
  propmat_operator& operator=(const propmat_operator&) = default;
  propmat_operator& operator=(propmat_operator&&) = default;

  propmat_operator(std::shared_ptr<AbsorptionBands> bands,
                   std::shared_ptr<AtmPoint> atm,
                   std::shared_ptr<ArrayOfCIARecord> cia,
                   std::shared_ptr<ArrayOfXsecRecord> xsec,
                   lbl::zeeman::pol pol = lbl::zeeman::pol::no);

  std::pair<Propmat, Stokvec> operator()(const Numeric frequency,
                                         const Vector2 los) const;
};  // struct propmat_operator
}  // namespace fwd
