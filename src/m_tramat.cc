#include <workspace.h>

void spectral_tramat_pathFromPath(
    TransmittanceMatrix& spectral_tramat_path,
    const ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfPropmatMatrix& spectral_propmat_jac_path,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& atm_path,
    const SurfaceField& surf_field,
    const JacobianTargets& jac_targets,
    const TransmittanceOption& rte_option,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  const Index it = jac_targets.target_position(AtmKey::t);
  const Vector r = distance(ray_path, surf_field.ellipsoid);
  Tensor3 dr(2, r.size(), jac_targets.target_count(), 0.0);
  if (hse_derivative and it >= 0) {
    for (Size ip = 0; ip < r.size(); ip++) {
      dr[0, ip, it] = r[ip] / (2.0 * atm_path[ip].temperature);
      dr[1, ip, it] = r[ip] / (2.0 * atm_path[ip + 1].temperature);
    }
  }

  spectral_tramat_path.init(
      spectral_propmat_path, spectral_propmat_jac_path, r, dr, rte_option);
}
ARTS_METHOD_ERROR_CATCH
