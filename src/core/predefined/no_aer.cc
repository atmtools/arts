/*

This file contains all the functions that are necessary to not compile
when AER is disabled

*/

#include "predef.h"

#define NOPE {ARTS_USER_ERROR("Compiled with -DENABLE_ARTS_LGPL=1")}

namespace Absorption::PredefinedModel {
namespace MT_CKD400 {
void compute_foreign_h2o(PropmatVector &,
                         const Vector &,
                         const AtmPoint &,
                         const WaterData &) NOPE;
void compute_self_h2o(PropmatVector &,
                      const Vector &,
                      const AtmPoint &,
                      const WaterData &) NOPE;
}  // namespace MT_CKD400
namespace MT_CKD430 {
void compute_foreign_h2o(PropmatVector &,
                         const Vector &,
                         const AtmPoint &,
                         const WaterData &) NOPE;
void compute_self_h2o(PropmatVector &,
                      const Vector &,
                      const AtmPoint &,
                      const WaterData &) NOPE;
}  // namespace MT_CKD430

namespace MT_CKD100 {
void oxygen_cia(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
void oxygen_v0v0(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
void oxygen_v0v1(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
}  // namespace MT_CKD100

namespace MT_CKD252 {
void carbon_dioxide(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
void oxygen_vis(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
void nitrogen_fun(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
void nitrogen_rot(PropmatVector &, const Vector &, const AtmPoint &) NOPE;
}  // namespace MT_CKD252

namespace CKDMT350 {
void compute_self_h2o(PropmatVector &, const Vector &, const AtmPoint &) NOPE;

void compute_foreign_h2o(PropmatVector &,
                         const Vector &,
                         const AtmPoint &) NOPE;
}  // namespace CKDMT350

namespace CKDMT320 {
void compute_self_h2o(PropmatVector &, const Vector &, const AtmPoint &) NOPE;

void compute_foreign_h2o(PropmatVector &,
                         const Vector &,
                         const AtmPoint &) NOPE;
}  // namespace CKDMT320
}  // namespace Absorption::PredefinedModel
