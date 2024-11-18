#include <atm.h>
#include <rtepack.h>

#include "predef_data.h"

namespace Absorption::PredefinedModel {

namespace ELL07 {
void compute(PropmatVector& propmat_clearsky,
             const Vector& f_grid,
             const AtmPoint& atm_point);
}  // namespace ELL07

namespace MPM89 {
void water(PropmatVector& propmat_clearsky,
           const Vector& f_grid,
           const AtmPoint& atm_point);

void oxygen(PropmatVector& propmat_clearsky,
            const Vector& f_grid,
            const AtmPoint& atm_point);
}  // namespace MPM89

namespace MPM93 {
void nitrogen(PropmatVector& propmat_clearsky,
              const Vector& f_grid,
              const AtmPoint& atm_point);
}  // namespace MPM93

namespace MPM2020 {
void compute(PropmatVector& propmat_clearsky,
             const Vector& f_grid,
             const AtmPoint& atm_point);
}  // namespace MPM2020
namespace PWR98 {
void water(PropmatVector& propmat_clearsky,
           const Vector& f_grid,
           const AtmPoint& atm_point);

void oxygen(PropmatVector& propmat_clearsky,
            const Vector& f_grid,
            const AtmPoint& atm_point);
}  // namespace PWR98
namespace TRE05 {
void oxygen(PropmatVector& propmat_clearsky,
            const Vector& f_grid,
            const AtmPoint& atm_point);
}  // namespace TRE05
namespace Standard {
void water_self(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const AtmPoint& atm_point);

void water_foreign(PropmatVector& propmat_clearsky,
                   const Vector& f_grid,
                   const AtmPoint& atm_point);

void nitrogen(PropmatVector& propmat_clearsky,
              const Vector& f_grid,
              const AtmPoint& atm_point);

void oxygen(PropmatVector& propmat_clearsky,
            const Vector& f_grid,
            const AtmPoint& atm_point);
}  // namespace Standard

namespace MT_CKD100 {
void oxygen_cia(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const AtmPoint& atm_point);

void oxygen_v0v0(PropmatVector& propmat_clearsky,
                 const Vector& f_grid,
                 const AtmPoint& atm_point);

void oxygen_v0v1(PropmatVector& propmat_clearsky,
                 const Vector& f_grid,
                 const AtmPoint& atm_point);
}  // namespace MT_CKD100

namespace MT_CKD252 {
void carbon_dioxide(PropmatVector& propmat_clearsky,
                    const Vector& f_grid,
                    const AtmPoint& atm_point);

void oxygen_vis(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const AtmPoint& atm_point);

void nitrogen_fun(PropmatVector& propmat_clearsky,
                  const Vector& f_grid,
                  const AtmPoint& atm_point);

void nitrogen_rot(PropmatVector& propmat_clearsky,
                  const Vector& f_grid,
                  const AtmPoint& atm_point);
}  // namespace MT_CKD252

namespace PWR20xx {
void compute_h2o_2021(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const AtmPoint& atm_point);

void compute_o2_2021(PropmatVector& propmat_clearsky,
                     const Vector& f_grid,
                     const AtmPoint& atm_point);

void compute_h2o_2022(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const AtmPoint& atm_point);

void compute_o2_2022(PropmatVector& propmat_clearsky,
                     const Vector& f_grid,
                     const AtmPoint& atm_point);

void compute_n2(PropmatVector& propmat_clearsky,
                const Vector& f_grid,
                const AtmPoint& atm_point);
}  // namespace PWR20xx

namespace CKDMT350 {
void compute_self_h2o(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const AtmPoint& atm_point);

void compute_foreign_h2o(PropmatVector& propmat_clearsky,
                         const Vector& f_grid,
                         const AtmPoint& atm_point);
}  // namespace CKDMT350

namespace CKDMT320 {
void compute_self_h2o(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const AtmPoint& atm_point);

void compute_foreign_h2o(PropmatVector& propmat_clearsky,
                         const Vector& f_grid,
                         const AtmPoint& atm_point);
}  // namespace CKDMT320

namespace MT_CKD400 {
void compute_foreign_h2o(PropmatVector& propmat_clearsky,
                         const Vector& f_grid,
                         const AtmPoint& atm_point,
                         const WaterData& data);

void compute_self_h2o(PropmatVector& propmat_clearsky,
                      const Vector& f_grid,
                      const AtmPoint& atm_point,
                      const WaterData& data);
}  // namespace MT_CKD400

}  // namespace Absorption::PredefinedModel
