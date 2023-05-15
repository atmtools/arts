/*

This file contains all the functions that are necessary to not compile
when AER is disabled

*/

#include "predef.h"

#define NOPE                                                                   \
  { ARTS_USER_ERROR("Compiled wiht -DENABLE_ARTS_LGPL=1") }

namespace Absorption::PredefinedModel {
namespace MT_CKD400 {
void compute_foreign_h2o(PropagationMatrix &, const Vector &, const Numeric &,
                         const Numeric &, const Numeric &,
                         const WaterData &) NOPE;
void compute_self_h2o(PropagationMatrix &, const Vector &, const Numeric &,
                      const Numeric &, const Numeric &, const WaterData &) NOPE;
} // namespace MT_CKD400

namespace MT_CKD100 {
void oxygen_cia(PropagationMatrix &, const Vector &, const Numeric,
                const Numeric, const Numeric) NOPE;
void oxygen_v0v0(PropagationMatrix &, const Vector &, const Numeric,
                 const Numeric, const Numeric, const Numeric) NOPE;
void oxygen_v0v1(PropagationMatrix &, const Vector &, const Numeric,
                 const Numeric, const Numeric) NOPE;
} // namespace MT_CKD100

namespace MT_CKD252 {
void carbon_dioxide(PropagationMatrix &, const Vector &, const Numeric,
                    const Numeric, const Numeric) NOPE;
void oxygen_vis(PropagationMatrix &, const Vector &, const Numeric,
                const Numeric, const Numeric) NOPE;
void nitrogen_fun(PropagationMatrix &, const Vector &, const Numeric,
                  const Numeric, const Numeric, const Numeric,
                  const Numeric) NOPE;
void nitrogen_rot(PropagationMatrix &, const Vector &, const Numeric,
                  const Numeric, const Numeric, const Numeric,
                  const Numeric) NOPE;
} // namespace MT_CKD252

namespace CKDMT350 {
void compute_self_h2o(PropagationMatrix &, const Vector &, const Numeric &,
                      const Numeric &, const Numeric &) NOPE;

void compute_foreign_h2o(PropagationMatrix &, const Vector &, const Numeric &,
                         const Numeric &, const Numeric &) NOPE;
} // namespace CKDMT350
} // namespace Absorption::PredefinedModel
