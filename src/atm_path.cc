#include "atm_path.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "auto_md.h"
#include "matpack_view.h"
#include "propagationmatrix.h"
#include "species_tags.h"
#include "transmissionmatrix.h"

#include <algorithm>
#include <array>

ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &atm_path,
                                 const Ppath &ppath) {
  atm_path.resize(ppath.np);
  return atm_path;
}

void forward_atm_path(ArrayOfAtmPoint &atm_path, const Ppath &ppath,
                      const AtmField &atm) {
  std::transform(ppath.pos.begin(), ppath.pos.end(), atm_path.begin(),
                 [&](const auto &path_pos) {
                   return atm.at(path_pos[0], path_pos[1], path_pos[2]);
                 });
}

ArrayOfAtmPoint forward_atm_path(const Ppath &ppath, const AtmField &atm) {
  ArrayOfAtmPoint atm_path(ppath.np);
  forward_atm_path(atm_path, ppath, atm);
  return atm_path;
}

ArrayOfVector &path_freq_resize(ArrayOfVector &path_freq,
                                const Vector &main_freq,
                                const ArrayOfAtmPoint &atm_path) {
  path_freq.resize(atm_path.nelem());
  for (auto &v : path_freq)
    v.resize(main_freq.nelem());
  return path_freq;
}

void forward_ppath_freq(ArrayOfVector &path_freq, const Vector &main_freq,
                        const Ppath &ppath, const ArrayOfAtmPoint &atm_path,
                        const Numeric along_path_atm_speed) {
  auto dot_prod = [&](Index ip) {
    const auto &[u, v, w] = atm_path[ip].wind;
    const auto &za = ppath.los(ip, 0);
    const auto &aa = ppath.los(ip, 1);
    const auto f = std::hypot(u, v, w);
    const auto za_f = std::acos(w / f);
    const auto aa_f = std::atan2(u, v);
    const auto za_p = Conversion::deg2rad(180 - za);
    const auto aa_p = Conversion::deg2rad(aa + 180);
    return f * (std::cos(za_f) * std::cos(za_p) +
                std::sin(za_f) * std::sin(za_p) * std::cos(aa_f - aa_p));
  };

  for (Index ip = 0; ip < atm_path.nelem(); ip++) {
    std::transform(
        main_freq.begin(), main_freq.end(), path_freq[ip].begin(),
        [fac = 1.0 - (along_path_atm_speed +
                      (atm_path[ip].zero_wind() ? 0.0 : dot_prod(ip))) /
                         Constant::speed_of_light](const auto &f) {
          return fac * f;
        });
  }
}

ArrayOfVector forward_ppath_freq(const Vector &main_freq, const Ppath &ppath,
                                 const ArrayOfAtmPoint &atm_path,
                                 const Numeric along_path_atm_speed) {
  ArrayOfVector path_freq(atm_path.nelem(), Vector(main_freq.nelem()));
  forward_ppath_freq(path_freq, main_freq, ppath, atm_path,
                     along_path_atm_speed);
  return path_freq;
}

std::tuple<ArrayOfPropagationMatrix, ArrayOfRadiationVector,
           ArrayOfArrayOfPropagationMatrix, ArrayOfArrayOfRadiationVector>
forward_propmat(Workspace &ws, const Agenda &propmat_clearsky_agenda,
                const Ppath &ppath, const ArrayOfVector &path_freq,
                const ArrayOfAtmPoint &atm_path,
                const ArrayOfRetrievalQuantity &jacobian_quantities) {
  const static ArrayOfSpeciesTag empty_tag_list(0);

  const Index np = atm_path.nelem();
  const Index nq = jacobian_quantities.nelem();
  ArrayOfPropagationMatrix K(np);
  ArrayOfArrayOfPropagationMatrix dK(np);
  StokesVector S;
  ArrayOfStokesVector dS;

  ArrayOfRadiationVector J(np);
  ArrayOfArrayOfRadiationVector dJ(np);

  for (Index ip = 0; ip < np; ip++) {
    propmat_clearsky_agendaExecute(
        ws, K[ip], S, dK[ip], dS, jacobian_quantities, empty_tag_list,
        path_freq[ip], Vector{atm_path[ip].mag}, Vector{ppath.los[ip]},
        atm_path[ip].pressure, atm_path[ip].temperature, {}, {},
        propmat_clearsky_agenda);
  }

  return {std::move(K), std::move(J), std::move(dK), std::move(dJ)};
}