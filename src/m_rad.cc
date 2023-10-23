#include <atm.h>
#include <new_jacobian.h>
#include <physics_funcs.h>
#include <ppath_struct.h>
#include <rtepack.h>
#include <surf.h>

void dradEmpty(StokvecMatrix &drad,
               const Vector &f_grid,
               const JacobianTargets &jacobian_targets) {
  drad.resize(jacobian_targets.size(), f_grid.nelem());
  drad = Stokvec{0.0, 0.0, 0.0, 0.0};
}

void dradFromPropagation(StokvecMatrix &drad,
                         const ArrayOfStokvecMatrix &ppvar_drad,
                         const MuelmatVector &cumulative_transmission,
                         const StokvecMatrix &background_drad,
                         const JacobianTargets &jacobian_targets,
                         const AtmField &atm_field,
                         const Ppath &ppath) {
  const auto np = ppvar_drad.size();
  const auto nj = jacobian_targets.size();
  const auto nf = cumulative_transmission.size();

  ARTS_USER_ERROR_IF(static_cast<Size>(ppath.pos.nrows()) != np,
                     "ppath.pos must have same number of rows as the size of ppvar_drad")
  ARTS_USER_ERROR_IF(ppath.pos.ncols() != 3, "ppath.pos must have 3 columns")

  ARTS_USER_ERROR_IF(
      static_cast<Size>(background_drad.nrows()) != nj,
      "background_drad must have same number of rows as the "
      "size of jacobian_targets")
  ARTS_USER_ERROR_IF(
      background_drad.ncols() != nf,
      "background_drad must have same number of colums as the "
      "size of cumulative_transmission")

  for (auto &dr : ppvar_drad) {
    ARTS_USER_ERROR_IF(
        static_cast<Size>(dr.nrows()) != nj,
        "ppvar_drad elements must have same number of rows as the size of "
        "jacobian_targets")
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf,
        "ppvar_drad elements must have same number of columns as the size of "
        "cumulative_transmission")
  }

  //! The radiance derivative shape is the background shape
  drad.resize(background_drad.shape());

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Size i = 0; i < nj; i++) {
    std::transform(cumulative_transmission.begin(),
                   cumulative_transmission.end(),
                   background_drad[i].begin(),
                   drad[i].begin(),
                   std::multiplies<>());
  }

  //! The altitude, latitude and longitude vectors must be copied because of how atm_field works
  const auto [alt, lat, lon] = [&]() {
    const auto n = ppath.pos.nrows();
    std::array<Vector, 3> g;
    for (auto &v : g) v.resize(n);
    for (Index i = 0; i < n; i++) {
      g[0][i] = ppath.pos[i][0];
      g[1][i] = ppath.pos[i][1];
      g[2][i] = ppath.pos[i][2];
    }
    return g;
  }();

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm) {
    ARTS_USER_ERROR_IF(not atm_field.contains(atm_block.type),
                       "No ",
                       atm_block.type,
                       " in atm_field but in jacobian_targets")
    const auto &data = atm_field[atm_block.type];
    const auto weights = data.flat_weights(alt, lat, lon);
    ARTS_ASSERT(weights.size() == np)

    for (Size j = 0; j < np; j++) {
      for (auto &w : weights[j]) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.start;
          ARTS_ASSERT(i < nj)
          std::transform(
              ppvar_drad[j][i].begin(),
              ppvar_drad[j][i].end(),
              drad[i].begin(),
              drad[i].begin(),
              [x = w.second](auto &a, auto &b) { return x * a + b; });
        }
      }
    }
  }
}

void radFromPropagation(StokvecVector& rad,
                        const ArrayOfStokvecVector& ppvar_rad) {
  ARTS_USER_ERROR_IF(ppvar_rad.empty(), "Empty ppvar_rad")
  rad = ppvar_rad.front();
}
