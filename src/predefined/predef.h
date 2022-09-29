#include <propagationmatrix.h>

#include "predef_data.h"

namespace Absorption::PredefinedModel {

namespace ELL07 {
void compute(PropagationMatrix& propmat_clearsky,
             const Vector& f_grid,
             const Numeric t,
             const Numeric lwc);
}  // namespace ELL07

namespace MPM89 {
void water(PropagationMatrix& propmat_clearsky,
           const Vector& f_grid,
           const Numeric p_pa,
           const Numeric t,
           const Numeric h2o) noexcept;

void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric o2,
            const Numeric h2o);
}  // namespace MPM89

namespace MPM2020 {
void compute(PropagationMatrix& propmat_clearsky,
             const Vector& f_grid,
             const Numeric& p_pa,
             const Numeric& t,
             const Numeric& oxygen_vmr) noexcept;
}  // namespace MPM2020
namespace PWR98 {
void water(PropagationMatrix& propmat_clearsky,
           const Vector& f_grid,
           const Numeric p_pa,
           const Numeric t,
           const Numeric h2o) noexcept;

void oxygen(PropagationMatrix& propmat_clearsky,
            const Vector& f_grid,
            const Numeric p_pa,
            const Numeric t,
            const Numeric o2,
            const Numeric h2o);
}  // namespace PWR98

namespace MT_CKD100 {
void oxygen_cia(PropagationMatrix& propmat_clearsky,
                const Vector& f_grid,
                const Numeric p,
                const Numeric Tave,
                const Numeric vmr);
void oxygen_v0v0(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric p,
                 const Numeric Tave,
                 const Numeric vmr,
                 const Numeric n2);
void oxygen_v0v1(PropagationMatrix& propmat_clearsky,
                 const Vector& f_grid,
                 const Numeric p,
                 const Numeric Tave,
                 const Numeric vmr);
}  // namespace MT_CKD100

namespace MT_CKD252 {
void carbon_dioxide(PropagationMatrix& propmat_clearsky,
                    const Vector& f_grid,
                    const Numeric p,
                    const Numeric Tave,
                    const Numeric vmr);
void oxygen_vis(PropagationMatrix& propmat_clearsky,
                const Vector& f_grid,
                const Numeric p,
                const Numeric Tave,
                const Numeric vmr);
void nitrogen_fun(PropagationMatrix& propmat_clearsky,
                  const Vector& f_grid,
                  const Numeric p_pa,
                  const Numeric t,
                  const Numeric vmr,
                  const Numeric h2o,
                  const Numeric o2);
void nitrogen_rot(PropagationMatrix& propmat_clearsky,
                  const Vector& f_grid,
                  const Numeric p_pa,
                  const Numeric t,
                  const Numeric vmr,
                  const Numeric h2o,
                  const Numeric o2);
}  // namespace MT_CKD252

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

namespace Hitran::MTCKD {
void compute_foreign_h2o(PropagationMatrix& propmat_clearsky,
                         const Vector& f_grid,
                         const Numeric& P,
                         const Numeric& T,
                         const Numeric& vmrh2o,
                         const WaterData& data);
void compute_self_h2o(PropagationMatrix& propmat_clearsky,
                      const Vector& f_grid,
                      const Numeric& P,
                      const Numeric& T,
                      const Numeric& vmrh2o,
                      const WaterData& data);
}  // namespace Hitran::MTCKD

}  // namespace Absorption::PredefinedModel