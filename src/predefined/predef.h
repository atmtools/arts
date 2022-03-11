#include <propagationmatrix.h>

namespace Absorption::PredefinedModel {
namespace MPM2020 {
void compute(PropagationMatrix& propmat_clearsky,
             const Vector& f_grid,
             const Numeric& p_pa,
             const Numeric& t,
             const Numeric& oxygen_vmr) noexcept;
}  // namespace MPM2020

namespace PWR2021 {
void compute_h2o(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric& p_pa,
                 const Numeric& t,
                 const Numeric& h2o_vmr) noexcept;
void compute_o2(PropagationMatrix& propmat_clearsky,
                const Vector& f_grid,
                const Numeric& p_pa,
                const Numeric& t,
                const Numeric& o2_vmr,
                const Numeric& h2o_vmr) noexcept;
}  // namespace PWR2021

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