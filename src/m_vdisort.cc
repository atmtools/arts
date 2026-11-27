#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <geodetic.h>
#include <lagrange_interp.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <subsurface_methods.h>
#include <vdisort_settings.h>
#include <workspace.h>

/* The VDISORT counterparts of the DISORT solver methods (m_disort.cc), one for one, with Stokes vectors. */

void vdisort_spectral_rad_fieldCalc(VDisortRadiance&       vdisort_spectral_rad_field,
                                    ZenGriddedField1&      quadrature,
                                    const VDisortSettings& vdisort_settings,
                                    const AziGrid&         phis) try {
  ARTS_TIME_REPORT

  const Index nv    = vdisort_settings.frequency_count();
  const Index np    = vdisort_settings.layer_count();
  const Index nquad = vdisort_settings.quadrature_dimension;
  const Index nphi  = phis.size();

  vdisort::main_data dis = vdisort_settings.init();

  //! Supplementary outputs
  quadrature = vdisort_quadrature_weights(dis);

  //! Main output
  vdisort_spectral_rad_field.resize(
      vdisort_settings.freq_grid, vdisort_settings.alt_grid, phis, ZenGrid{quadrature.grid<0>()});

  // The propagation azimuths as VDISORT's
  Vector phi(nphi);
  for (Index p = 0; p < nphi; p++) phi[p] = vdisort_azimuth(phis[p]);

  Tensor4 u(np, nphi, nquad, 4);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis, u)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      vdisort_settings.set(dis, iv);
      dis.gridded_u(u, phi);
      for (Index l = 0; l < np; l++)
        for (Index p = 0; p < nphi; p++)
          for (Index i = 0; i < nquad; i++)
            vdisort_spectral_rad_field.data[iv, l, p, i] =
                Stokvec{u[l, p, i, 0], u[l, p, i, 1], u[l, p, i, 2], u[l, p, i, 3]};
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  vdisort_spectral_rad_field.sort(dis.mu());

  if (error.size()) arts::user_error("Error occurred in vdisort-spectral:\n{}", error);
}
ARTS_METHOD_ERROR_CATCH

void spectral_flux_components_profileVDisortCalc(FluxComponentsProfile& spectral_flux_components_profile,
                                                 const VDisortSettings& vdisort_settings) try {
  ARTS_TIME_REPORT

  const Index nv = vdisort_settings.frequency_count();

  spectral_flux_components_profile.resize(vdisort_settings.freq_grid, vdisort_settings.alt_grid);

  vdisort::main_data dis = vdisort_settings.init();

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      vdisort_settings.set(dis, iv);

      dis.gridded_flux(spectral_flux_components_profile.up[iv],
                       spectral_flux_components_profile.down_diffuse[iv],
                       spectral_flux_components_profile.down_direct[iv],
                       spectral_flux_components_profile.dfdt[iv]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (error.size()) arts::user_error("Error occurred in vdisort:\n{}", error);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radFromVDisort(StokvecVector&              spectral_rad,
                             const VDisortRadiance&      vdisort_spectral_rad_field,
                             const PropagationPathPoint& ray_point) {
  using namespace lagrange_interp;
  using aa_cyc_t = cycler<0.0, 360.0>;

  ARTS_TIME_REPORT

  const auto& f_grid   = vdisort_spectral_rad_field.freq_grid;
  const auto& alt_grid = vdisort_spectral_rad_field.alt_grid;
  const auto& azi_grid = vdisort_spectral_rad_field.azi_grid;
  const auto& zen_grid = vdisort_spectral_rad_field.zen_grid;
  const auto& data     = vdisort_spectral_rad_field.data;

  const Size nf = f_grid.size();
  spectral_rad.resize(nf);

  if (nf == 0) return;

  if (alt_grid.size() < 2) arts::user_error("VDISORT altitude grid must have at least two points");

  const Numeric z     = ray_point.altitude();
  const bool    above = alt_grid.size() < 2 or z >= alt_grid[1];
  const bool    below = z <= alt_grid.back();

  // One Stokes component at a time, as spectral_radFromDisort interpolates I
  Tensor4 component(data.shape());
  for (Index s = 0; s < 4; s++) {
    for (Size i = 0; i < data.size(); i++) component.data_handle()[i] = data.data_handle()[i][s];

    if (above or below) {
      std::visit(
          [&, idx = above ? 0 : alt_grid.size() - 2](auto&& aa, auto&& za) {
            for (Size v = 0; v < nf; v++) spectral_rad[v][s] = interp(component[v, idx], aa, za);
          },
          variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
          variant_lag(zen_grid, ray_point.los[0]));
    } else {
      std::visit(
          [&](auto&& alt, auto&& aa, auto&& za) {
            for (Size v = 0; v < nf; v++) spectral_rad[v][s] = interp(component[v], alt, aa, za);
          },
          variant_lag(std::span{alt_grid}.subspan(1), z),
          variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
          variant_lag(zen_grid, ray_point.los[0]));
    }
  }
}

void vdisort_spectral_rad_fieldApplyUnit(VDisortRadiance&                         vdisort_spectral_rad_field,
                                         const PropagationPathPoint&              ray_point,
                                         const SpectralRadianceTransformOperator& spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecVector spectral_rad(vdisort_spectral_rad_field.freq_grid.size());
  StokvecMatrix spectral_rad_jac(0, vdisort_spectral_rad_field.freq_grid.size());

  const Index nv = vdisort_spectral_rad_field.data.nbooks();
  const Index np = vdisort_spectral_rad_field.data.npages();
  const Index nz = vdisort_spectral_rad_field.data.nrows();
  const Index nq = vdisort_spectral_rad_field.data.ncols();

  std::string error{};

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3) firstprivate(spectral_rad, spectral_rad_jac)
  for (Index i = 0; i < np; i++) {
    for (Index j = 0; j < nz; j++) {
      for (Index k = 0; k < nq; k++) {
        try {
          for (Index v = 0; v < nv; v++) spectral_rad[v] = vdisort_spectral_rad_field.data[v, i, j, k];

          spectral_rad_transform_operator(
              spectral_rad, spectral_rad_jac, vdisort_spectral_rad_field.freq_grid, ray_point);

          for (Index v = 0; v < nv; v++) vdisort_spectral_rad_field.data[v, i, j, k] = spectral_rad[v];
        } catch (std::exception& e) {
#pragma omp critical
          if (error.empty()) error = e.what();
        }
      }
    }
  }

  if (not error.empty()) arts::user_error(error);
}
ARTS_METHOD_ERROR_CATCH

namespace {
//! The polarization basis ARTS gives radiation propagating along the unit ENU vector k of propagation azimuth [deg]
rtepack::polarization_basis propagation_basis(const Vector3& k, const Numeric azimuth) {
  return rtepack::meridional_basis(k, Vector3{Conversion::cosd(azimuth), -Conversion::sind(azimuth), 0.0});
}
}  // namespace

void spectral_radSubsurfaceVDisort(const Workspace&               ws,
                                   StokvecVector&                 spectral_rad,
                                   StokvecMatrix&                 spectral_rad_jac,
                                   const AscendingGrid&           freq_grid,
                                   const AtmField&                atm_field,
                                   const SurfaceField&            surf_field,
                                   const SubsurfaceField&         subsurf_field,
                                   const SubsurfaceSpectralField& subsurf_spectral_field,
                                   const JacobianTargets&         jac_targets,
                                   const PropagationPathPoint&    ray_point,
                                   const Index&                   quadrature_dimension,
                                   const Index&                   fourier_mode_dimension,
                                   const Index&                   legendre_polynomial_dimension,
                                   const Agenda&                  spectral_rad_observer_agenda,
                                   const Agenda&                  subsurf_vdisort_settings_agenda,
                                   const DescendingGrid&          depth_profile) try {
  ARTS_TIME_REPORT

  if (not jac_targets.surf.empty() or not jac_targets.subsurf.empty() or not jac_targets.subsurf_spectral.empty())
    arts::user_error(
        "The derivatives of the subsurface radiance by surface, subsurface and spectral subsurface targets are not "
        "implemented");

  const Size nf = freq_grid.size();

  // The outward facet normal and the propagation directions in east-north-up: the outgoing ray, and the sky's ray
  // that the facet reflects into it.  The normal of a surface point is the line of sight into the surface along it.
  const Vector3 normal = -1.0 * los2enu(surf_field.at(ray_point.latitude(), ray_point.longitude()).normal);
  const Vector3 k_out  = los2enu(ray_point.los);
  const Numeric mu     = dot(normal, k_out);
  if (not(mu > 0.0))
    arts::user_error("The ray leaves the facet from below: the cosine between it and the facet normal is {}", mu);
  const Vector3 sky_direction = 2.0 * mu * normal - k_out;
  const Vector3 k_inc         = -1.0 * sky_direction;
  const Vector2 sky_los       = enu2los(sky_direction);

  // The radiation from the sky that the facet reflects, with the sun if it is there
  StokvecVector               sky;
  StokvecMatrix               sky_jac;
  ArrayOfPropagationPathPoint sky_path;
  spectral_rad_observer_agendaExecute(ws,
                                      sky,
                                      sky_jac,
                                      sky_path,
                                      freq_grid,
                                      jac_targets,
                                      ray_point.pos,
                                      sky_los,
                                      atm_field,
                                      surf_field,
                                      subsurf_field,
                                      spectral_rad_observer_agenda);

  // The radiation within the subsurface, at the interface
  ArrayOfPropagationPathPoint subsurf_path;
  ray_pathFromPointAndDepth(subsurf_path, ray_point, depth_profile);
  VDisortSettings vdisort_settings;
  subsurf_vdisort_settings_agendaExecute(ws,
                                         vdisort_settings,
                                         freq_grid,
                                         subsurf_path,
                                         quadrature_dimension,
                                         fourier_mode_dimension,
                                         legendre_polynomial_dimension,
                                         subsurf_vdisort_settings_agenda);
  if (stdr::any_of(vdisort_settings.single_scattering_albedo | by_elem, Cmp::gt(0.0)))
    arts::user_error(
        "Scattering within the subsurface needs its phase matrices at the refracted angle, which are not "
        "implemented");

  const ComplexVector n = subsurface_refractive_index(
      subsurface_spectral_point(subsurf_spectral_field, surf_field, freq_grid, subsurf_path.front()), nf);

  // The facet's bases (p, s) of the incoming and outgoing rays, s across the plane of incidence, and ARTS's bases.  At
  // normal incidence the plane is undefined and any s across the rays will do.
  const rtepack::polarization_basis in  = propagation_basis(k_inc, sky_los[1] + 180.0);
  const rtepack::polarization_basis out = propagation_basis(k_out, ray_point.los[1]);
  Vector3                           s   = cross(k_inc, normal);
  const Numeric                     ns  = hypot(s);
  s                                     = ns > 1e-12 ? s / ns : out.second;
  const rtepack::muelmat to_facet       = rtepack::basis_rotation(in, {cross(s, k_inc), s});
  const rtepack::muelmat from_facet     = rtepack::basis_rotation({cross(s, k_out), s}, out);

  spectral_rad.resize(nf);
  spectral_rad_jac.resize(jac_targets.x_size(), nf);

  vdisort::main_data               dis = vdisort_settings.init();
  const vdisort::phase_matrix_data phase(2,
                                         vdisort_settings.fourier_mode_dimension,
                                         vdisort_settings.layer_count(),
                                         1,
                                         vdisort_settings.quadrature_dimension,
                                         rtepack::muelmat{0.0});
  const Numeric                    sin2 = 1.0 - mu * mu;

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Size iv = 0; iv < nf; iv++) {
    try {
      // The radiance within that refracts into the outgoing ray, at the interface, propagating towards it.  It is
      // VDISORT's downward radiance, in the plane of incidence with Q = I_p - I_s.
      const Numeric nr = n[iv].real();
      if (not(sin2 < nr * nr))
        arts::user_error("The subsurface refractive index {} cannot refract the ray, at {} deg from the normal",
                         nr,
                         Conversion::acosd(mu));
      const Numeric        mu_inside = std::sqrt(1.0 - sin2 / (nr * nr));
      vdisort::user_u_data data;
      vdisort_settings.set(dis, iv);
      dis.u_user(data,
                 vdisort_settings.optical_thicknesses[iv, vdisort_settings.layer_count() - 1],
                 0.0,
                 Vector{-mu_inside},
                 phase);
      const Stokvec& inside = data.intensities[0];

      // The interface transmits each polarization with one minus its reflectance, which is the same from
      // either side
      const auto [Rv, Rh]       = fresnel(1.0, Complex{nr, 0.0}, Conversion::acosd(mu));
      const Numeric          Tp = 1.0 - std::norm(Rv);
      const Numeric          Ts = 1.0 - std::norm(Rh);
      const Numeric          a  = 0.5 * (Tp + Ts);
      const Numeric          b  = 0.5 * (Tp - Ts);
      const Stokvec          transmitted{a * inside.I() + b * inside.Q(), b * inside.I() + a * inside.Q(), 0.0, 0.0};
      const rtepack::muelmat R = from_facet * rtepack::fresnel_reflectance(Rv, Rh) * to_facet;
      spectral_rad[iv]         = R * sky[iv] + from_facet * transmitted;
      for (Index ix = 0; ix < spectral_rad_jac.nrows(); ix++) spectral_rad_jac[ix, iv] = R * sky_jac[ix, iv];
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (error.size()) arts::user_error("Error in the subsurface radiance:\n{}", error);
}
ARTS_METHOD_ERROR_CATCH
