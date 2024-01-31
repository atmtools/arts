#include "fwd_propmat.h"

#include "lbl_zeeman.h"

namespace fwd {
std::pair<Propmat, Stokvec> propmat_operator::operator()(
    const Numeric f, const Vector2 los) const {
  Propmat propmat{};
  Stokvec stokvec{};

  const auto [a, s] = lines(f);

  const auto propmat_pol = lbl::zeeman::norm_view(pol, atm->mag, los);

  propmat += lbl::zeeman::scale(propmat_pol, s);
  propmat.A() += cia(f).real() + predef(f).real() + hxsec(f).real();

  const auto srcvec = lbl::zeeman::scale(propmat_pol, s);
  stokvec += {srcvec.A(), srcvec.B(), srcvec.C(), srcvec.D()};

  return {propmat, stokvec};
}
}  // namespace fwd
