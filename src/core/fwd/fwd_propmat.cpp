#include "fwd_propmat.h"

#include "lbl_zeeman.h"

namespace fwd {
std::pair<Propmat, Stokvec> propmat_operator::operator()(
    const Numeric f, const Vector2 los) const {
  Propmat propmat{};
  Stokvec stokvec{};

  const auto [a, s] = lines(f);

  const auto pol = lines.polarization(los);

  propmat += lbl::zeeman::scale(pol, s);
  propmat.A() += cia(f).real() + predef(f).real() + hxsec(f).real();

  const auto srcvec = lbl::zeeman::scale(pol, s);
  stokvec += {srcvec.A(), srcvec.B(), srcvec.C(), srcvec.D()};

  return {propmat, stokvec};
}
}  // namespace fwd
