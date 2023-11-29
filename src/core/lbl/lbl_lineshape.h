#include <lbl_data.h>

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl {
void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               matpack::matpack_view<Stokvec, 2, false, true> dsv,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const std::span<const lbl::band>& bnds,
               const AtmPoint& atm,
               const Vector2 los = {});
}  // namespace lbl
