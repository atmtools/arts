#include "fwd_propmat.h"

#include "lbl_zeeman.h"

namespace fwd {
propmat_operator::propmat_operator(std::shared_ptr<AtmPoint> atm_,
                                   std::shared_ptr<AbsorptionBands> lines_,
                                   std::shared_ptr<ArrayOfCIARecord> cia,
                                   std::shared_ptr<ArrayOfXsecRecord> xsec,
                                   std::shared_ptr<PredefinedModelData> predef,
                                   Numeric ciaextrap,
                                   Index ciarobust,
                                   lbl::zeeman::pol pol_)
    : atm(std::move(atm_)),
      pol(pol_),
      lines(atm, std::move(lines_), pol),
      cia(atm, std::move(cia), ciaextrap, ciarobust),
      predef(atm, std::move(predef)),
      xsec(atm, std::move(xsec)) {}

std::pair<Propmat, Stokvec> propmat_operator::operator()(
    const Numeric f, const Vector2 los) const {
  Propmat propmat{};
  Stokvec stokvec{};

  const auto [a, s] = lines(f);

  const auto propmat_pol = lbl::zeeman::norm_view(pol, atm->mag, los);

  propmat += lbl::zeeman::scale(propmat_pol, s);
  propmat.A() += cia(f).real() + predef(f).real() + xsec(f).real();

  const auto srcvec = lbl::zeeman::scale(propmat_pol, s);
  stokvec += {srcvec.A(), srcvec.B(), srcvec.C(), srcvec.D()};

  return {propmat, stokvec};
}

void propmat_operator::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  lines.set_atm(atm);
  cia.set_atm(atm);
  predef.set_atm(atm);
  xsec.set_atm(atm);
}

void propmat_operator::set_pol(lbl::zeeman::pol pol_) {
  pol = pol_;
  lines.set_pol(pol);
}

void propmat_operator::set_ciaextrap(Numeric extrap) { cia.set_extrap(extrap); }

void propmat_operator::set_ciarobust(Index robust) { cia.set_robust(robust); }

void propmat_operator::set_bands(std::shared_ptr<AbsorptionBands> lines_) {
  lines.set_model(std::move(lines_));
}

void propmat_operator::set_cia(std::shared_ptr<ArrayOfCIARecord> cia_) {
  cia.set_model(std::move(cia_));
}

void propmat_operator::set_predef(
    std::shared_ptr<PredefinedModelData> predef_) {
  predef.set_model(std::move(predef_));
}

void propmat_operator::set_model(std::shared_ptr<ArrayOfXsecRecord> xsec_) {
  xsec.set_model(std::move(xsec_));
}
}  // namespace fwd
