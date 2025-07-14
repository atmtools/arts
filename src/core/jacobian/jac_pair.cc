#include "jac_pair.h"

#include "jac_log.h"
#include "jac_rel.h"

void make_logrelfit(Jacobian::AtmTarget& x, const AtmField& atm) {
  auto sorig = std::make_shared<Vector>(atm[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};
  const logfwd lfwd{.N = sorig->size()};
  const loginv linv{.N = sorig->size()};

  const pairinv<AtmField> inv{
      .vfunc_1 = linv, .vfunc_2 = rinv, .mfunc_1 = linv, .mfunc_2 = rinv};

  x.inverse_state = inv;

  x.inverse_jacobian = inv;

  x.transform_state = pairfwd<AtmField>{.vfunc_1 = rfwd, .vfunc_2 = lfwd};
}

void make_logrelfit(Jacobian::SurfaceTarget& x, const SurfaceField& surf) {
  auto sorig = std::make_shared<Vector>(surf[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};
  const logfwd lfwd{.N = sorig->size()};
  const loginv linv{.N = sorig->size()};

  x.inverse_state = pairinv<SurfaceField>{
      .vfunc_1 = rinv, .vfunc_2 = linv, .mfunc_1 = rinv, .mfunc_2 = linv};
  x.inverse_jacobian = pairinv<SurfaceField>{
      .vfunc_1 = rinv, .vfunc_2 = linv, .mfunc_1 = rinv, .mfunc_2 = linv};
  x.transform_state = pairfwd<SurfaceField>{.vfunc_1 = rfwd, .vfunc_2 = lfwd};
}

void make_logrelfit(Jacobian::SubsurfaceTarget& x,
                    const SubsurfaceField& subsurf) {
  auto sorig = std::make_shared<Vector>(subsurf[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};
  const logfwd lfwd{.N = sorig->size()};
  const loginv linv{.N = sorig->size()};

  x.inverse_state = pairinv<SubsurfaceField>{
      .vfunc_1 = rinv, .vfunc_2 = linv, .mfunc_1 = rinv, .mfunc_2 = linv};
  x.inverse_jacobian = pairinv<SubsurfaceField>{
      .vfunc_1 = rinv, .vfunc_2 = linv, .mfunc_1 = rinv, .mfunc_2 = linv};
  x.transform_state =
      pairfwd<SubsurfaceField>{.vfunc_1 = rfwd, .vfunc_2 = lfwd};
}
