#include <propagationmatrix.h>

namespace Absorption::PredefinedModel {
namespace MPM2020 {
void compute(PropagationMatrix& propmat_clearsky,
             const Vector& f_grid,
             const Numeric& p_pa,
             const Numeric& t,
             const Numeric& oxygen_vmr) noexcept;
}  // namespace MPM2020

namespace CKDMT350 {
void compute_self_h2o(PropagationMatrix& propmat_clearsky,
                      const Vector& f_grid,
                      const Numeric& p,
                      const Numeric& Tave,
                      const Numeric& vmrh2o) noexcept;
void compute_foreign_h2o(PropagationMatrix& propmat_clearsky,
                         const Vector& f_grid,
                         const Numeric& p,
                         const Numeric& Tave,
                         const Numeric& vmrh2o) noexcept;
}  // namespace CKDMT350
}  // namespace Absorption::PredefinedModel