#include "fwd_propmat.h"

#include <functional>
#include <numeric>

#include "lbl_zeeman.h"
#include "rtepack.h"

namespace fwd {
propmat::propmat(std::shared_ptr<AtmPoint> atm_,
                 std::shared_ptr<AbsorptionBands> lines_,
                 std::shared_ptr<ArrayOfCIARecord> cia,
                 std::shared_ptr<ArrayOfXsecRecord> xsec,
                 std::shared_ptr<PredefinedModelData> predef,
                 Numeric ciaextrap,
                 Index ciarobust)
    : atm(std::move(atm_)),
      lines(atm, std::move(lines_)),
      cia(atm, std::move(cia), ciaextrap, ciarobust),
      predef(atm, std::move(predef)),
      xsec(atm, std::move(xsec)) {}

std::pair<Propmat, Stokvec> propmat::operator()(const Numeric f,
                                                const Vector2 los) const {
  using namespace lbl::zeeman;

  const auto [ano, sno] = lines(f, pol::no);

  const std::array zres{
      lines(f, pol::sm), lines(f, pol::pi), lines(f, pol::sp)};

  const std::array zpol{
      norm_view(pol::sm, atm->mag, los),
      norm_view(pol::pi, atm->mag, los),
      norm_view(pol::sp, atm->mag, los),
  };

  return {std::transform_reduce(
              zpol.begin(),
              zpol.end(),
              zres.begin(),
              Propmat{cia(f).real() + predef(f).real() + xsec(f).real() +
                      ano.real()},
              std::plus<>(),
              [](const Propmat& a, const std::pair<Complex, Complex>& b) {
                return scale(a, b.first);
              }),
          std::transform_reduce(
              zpol.begin(),
              zpol.end(),
              zres.begin(),
              Stokvec{sno.real()},
              std::plus<>(),
              [](const Propmat& a, const std::pair<Complex, Complex>& b) {
                return absvec(scale(a, b.second));
              })};
}

void propmat::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  lines.set_atm(atm);
  cia.set_atm(atm);
  predef.set_atm(atm);
  xsec.set_atm(atm);
}

void propmat::set_ciaextrap(Numeric extrap) { cia.set_extrap(extrap); }

void propmat::set_ciarobust(Index robust) { cia.set_robust(robust); }

void propmat::set_bands(std::shared_ptr<AbsorptionBands> lines_) {
  lines.set_model(std::move(lines_));
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
