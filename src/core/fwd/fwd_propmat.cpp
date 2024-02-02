#include "fwd_propmat.h"

#include "lbl_zeeman.h"

namespace fwd {
propmat::propmat(std::shared_ptr<AtmPoint> atm_,
                 std::shared_ptr<AbsorptionBands> lines_,
                 std::shared_ptr<ArrayOfCIARecord> cia,
                 std::shared_ptr<ArrayOfXsecRecord> xsec,
                 std::shared_ptr<PredefinedModelData> predef,
                 Numeric ciaextrap,
                 Index ciarobust)
    : atm(std::move(atm_)),
      lines(atm, lines_, lbl::zeeman::pol::no),
      sm_lines(atm, lines_, lbl::zeeman::pol::sm),
      pi_lines(atm, lines_, lbl::zeeman::pol::pi),
      sp_lines(atm, std::move(lines_), lbl::zeeman::pol::sp),
      cia(atm, std::move(cia), ciaextrap, ciarobust),
      predef(atm, std::move(predef)),
      xsec(atm, std::move(xsec)) {}

std::pair<Propmat, Stokvec> propmat::operator()(const Numeric f,
                                                const Vector2 los) const {
  using namespace lbl::zeeman;

  Propmat propmat{};
  Stokvec stokvec{};

  const auto [noa, nos] = lines(f);
  const auto [spa, sps] = sp_lines(f);
  const auto [pia, pis] = pi_lines(f);
  const auto [sma, sms] = sm_lines(f);

  propmat.A() += noa.real() + cia(f).real() + predef(f).real() + xsec(f).real();

  const auto sp_propmat = norm_view(pol::sp, atm->mag, los);
  const auto pi_propmat = norm_view(pol::pi, atm->mag, los);
  const auto sm_propmat = norm_view(pol::sm, atm->mag, los);
  propmat +=
      scale(sp_propmat, spa) + scale(pi_propmat, pia) + scale(sm_propmat, sma);

  const auto sp_srcvec = absvec(scale(sp_propmat, sps));
  const auto pi_srcvec = absvec(scale(pi_propmat, pis));
  const auto sm_srcvec = absvec(scale(sm_propmat, sms));
  stokvec += sp_srcvec + pi_srcvec + sm_srcvec + nos.imag();

  return {propmat, stokvec};
}

void propmat::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  lines.set_atm(atm);
  sp_lines.set_atm(atm);
  sm_lines.set_atm(atm);
  pi_lines.set_atm(atm);
  cia.set_atm(atm);
  predef.set_atm(atm);
  xsec.set_atm(atm);
}

void propmat::set_ciaextrap(Numeric extrap) { cia.set_extrap(extrap); }

void propmat::set_ciarobust(Index robust) { cia.set_robust(robust); }

void propmat::set_bands(std::shared_ptr<AbsorptionBands> lines_) {
  lines.set_model(lines_);
  sp_lines.set_model(lines_);
  pi_lines.set_model(lines_);
  sm_lines.set_model(std::move(lines_));
}

void propmat::set_cia(std::shared_ptr<ArrayOfCIARecord> cia_) {
  cia.set_model(std::move(cia_));
}

void propmat::set_predef(std::shared_ptr<PredefinedModelData> predef_) {
  predef.set_model(std::move(predef_));
}

void propmat::set_model(std::shared_ptr<ArrayOfXsecRecord> xsec_) {
  xsec.set_model(std::move(xsec_));
}
}  // namespace fwd
