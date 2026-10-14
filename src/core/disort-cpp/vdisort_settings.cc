#include "vdisort_settings.h"

#include <arts_conversions.h>
#include <debug.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "vdisort_arts.h"

void VDisortRadiance::resize(AscendingGrid f_grid, DescendingGrid alt_grid_, AziGrid azi_grid_, ZenGrid zen_grid_) {
  freq_grid = std::move(f_grid);
  alt_grid  = std::move(alt_grid_);
  azi_grid  = std::move(azi_grid_);
  zen_grid  = std::move(zen_grid_);
  data.resize(freq_grid.size(), alt_grid.size() - 1, azi_grid.size(), zen_grid.size());
}

void VDisortRadiance::sort(const Vector& solver_mu) {
  ArrayOfIndex idx(solver_mu.size());
  stdr::iota(idx, 0);
  stdr::sort(idx, [&solver_mu](Index i, Index j) { return solver_mu[i] > solver_mu[j]; });

  const StokvecTensor4 data_copy = data;
  for (Size i = 0; i < solver_mu.size(); i++) data[joker, joker, joker, i] = data_copy[joker, joker, joker, idx[i]];
}

void VDisortSettings::resize(Index          quadrature_dimension_,
                             Index          legendre_polynomial_dimension_,
                             Index          fourier_mode_dimension_,
                             AscendingGrid  f_grid,
                             DescendingGrid alt_grid_) {
  quadrature_dimension          = quadrature_dimension_;
  legendre_polynomial_dimension = legendre_polynomial_dimension_;
  fourier_mode_dimension        = fourier_mode_dimension_;
  const Size nfreq              = f_grid.size();
  const Size nlay               = alt_grid_.size() - 1;

  freq_grid = std::move(f_grid);
  alt_grid  = std::move(alt_grid_);

  solar_source.resize(nfreq);
  solar_zenith_angle.resize(nfreq);
  solar_azimuth_angle.resize(nfreq);
  surface_brdf.resize(nfreq, 2, 0, quadrature_dimension / 2, quadrature_dimension / 2);
  surface_beam_brdf.resize(nfreq, 2, 0, quadrature_dimension / 2);
  surface_beam_mu.resize(nfreq);
  surface_beam_mu = std::numeric_limits<Numeric>::quiet_NaN();
  surface_specular.resize(nfreq, 0);
  optical_thicknesses.resize(nfreq, nlay);
  single_scattering_albedo.resize(nfreq, nlay);
  source_polynomial.resize(nfreq, nlay, 0);
  legendre_coefficients.resize(nfreq, nlay, legendre_polynomial_dimension, 6);
  upward_boundary_condition.resize(nfreq, 2, fourier_mode_dimension, quadrature_dimension / 2);
  downward_boundary_condition.resize(nfreq, 2, fourier_mode_dimension, quadrature_dimension / 2);
}

bool VDisortSettings::ok() const {
  const Index nfreq = freq_grid.size();
  const Index nlay  = alt_grid.size() - 1;
  const Index N     = quadrature_dimension / 2;

  return not(
      solar_source.shape() != std::array{nfreq} or solar_zenith_angle.shape() != std::array{nfreq} or
      solar_azimuth_angle.shape() != std::array{nfreq} or
      (surface_brdf.shape() != std::array{nfreq, Index{2}, surface_brdf.extent(2), N, N}) or
      (surface_beam_brdf.shape() != std::array{nfreq, Index{2}, surface_brdf.extent(2), N}) or
      (surface_beam_mu.shape() != std::array{nfreq}) or
      (surface_specular.shape() != std::array{nfreq, N} and surface_specular.shape() != std::array{nfreq, Index{0}}) or
      (optical_thicknesses.shape() != std::array{nfreq, nlay}) or
      (single_scattering_albedo.shape() != std::array{nfreq, nlay}) or
      (source_polynomial.shape() != std::array{nfreq, nlay, source_polynomial.ncols()}) or
      (legendre_coefficients.shape() != std::array{nfreq, nlay, legendre_coefficients.nrows(), Index{6}}) or
      (upward_boundary_condition.shape() != std::array{nfreq, Index{2}, fourier_mode_dimension, N}) or
      (downward_boundary_condition.shape() != std::array{nfreq, Index{2}, fourier_mode_dimension, N}) or
      legendre_polynomial_dimension > legendre_coefficients.nrows());
}

void VDisortSettings::check() const {
  ARTS_USER_ERROR_IF(not ok(),
                     R"-x-(Input is of incorrect size.

{:Bs,}

Also note that the reduced Legendre polynomial dimension is {}.  It must be at most {}.
)-x-",
                     *this,
                     legendre_polynomial_dimension,
                     legendre_coefficients.nrows());
}

vdisort::main_data VDisortSettings::init() const {
  check();
  return vdisort::main_data(alt_grid.size() - 1,
                            quadrature_dimension,
                            fourier_mode_dimension,
                            source_polynomial.ncols(),
                            surface_brdf.extent(2),
                            surface_specular.ncols() > 0);
}

Numeric vdisort_azimuth(Numeric aa) {
  Numeric phi = std::fmod(-Conversion::deg2rad(aa), Constant::two_pi);
  if (phi < 0.0) phi += Constant::two_pi;
  if (phi >= Constant::two_pi) phi = 0.0;
  return phi;
}

ZenGriddedField1 vdisort_quadrature_weights(const vdisort::main_data& dis) {
  Vector mu = dis.mu();
  stdr::sort(mu);
  std::transform(mu.begin(), mu.end(), mu.begin(), [](const Numeric& m) { return 180.0 - Conversion::acosd(m); });

  const Index N = dis.weights().size();
  ZenGriddedField1 quadrature{
      .data_name = "VDISORT quadrature weights", .data = mu, .grid_names = {"Zenith grid"}, .grids = {std::move(mu)}};
  quadrature.data[Range(0, N)] = dis.weights();
  quadrature.data[Range(N, N)] = dis.weights();
  return quadrature;
}

vdisort::main_data& VDisortSettings::set(vdisort::main_data& dis, Index iv) const {
  using Conversion::cosd;

  const Index nlay = layer_count();
  const Index NQ   = quadrature_dimension;
  const Index NF   = fourier_mode_dimension;
  const Index L    = legendre_polynomial_dimension;

  dis.surface(surface_brdf[iv], surface_beam_brdf[iv], surface_beam_mu[iv], surface_specular[iv]);

  dis.tau(optical_thicknesses[iv]);
  dis.solar_zenith() = cosd(solar_zenith_angle[iv]);
  dis.beam_azimuth() = vdisort_azimuth(solar_azimuth_angle[iv] + 180.0);
  dis.omega()        = single_scattering_albedo[iv];

  // The phase matrices of every layer on VDISORT's streams (and from the beam) from the Legendre moments
  const bool beam = solar_source[iv].I() != 0.0;
  const Vector mu{dis.mu()};
  Vector       mu_in(NQ + (beam ? 1 : 0));
  mu_in[Range(0, NQ)] = mu;
  if (beam) mu_in[NQ] = -dis.solar_zenith();

  rtepack::muelmat_tensor4 C(NF, nlay, NQ, NQ, rtepack::muelmat{0.0}), S = C;
  rtepack::muelmat_tensor3 Cb(NF, nlay, NQ, rtepack::muelmat{0.0}), Sb   = Cb;
  for (Index l = 0; l < nlay; l++) {
    if (single_scattering_albedo[iv, l] == 0.0) continue;
    const auto optics = vdisort::legendre_fourier_optics(legendre_coefficients[iv, l, Range(0, L)], mu, mu_in, NF);
    for (Index m = 0; m < NF; m++) {
      for (Index o = 0; o < NQ; o++) {
        for (Index i = 0; i < NQ; i++) {
          C[m, l, o, i] = optics.cosine[m, o, i];
          S[m, l, o, i] = optics.sine[m, o, i];
        }
        if (beam) {
          Cb[m, l, o] = optics.cosine[m, o, NQ];
          Sb[m, l, o] = optics.sine[m, o, NQ];
        }
      }
    }
  }
  dis.all_phase_matrices()      = vdisort::combine_phase_matrices(C, S);
  dis.all_beam_phase_matrices() = vdisort::combine_beam_phase_matrices(Cb, Sb);

  dis.upward_boundary()   = upward_boundary_condition[iv];
  dis.downward_boundary() = downward_boundary_condition[iv];
  dis.source_poly()       = source_polynomial[iv];
  dis.set_beam_source(solar_source[iv]);

  dis.update_all();

  return dis;
}

void xml_io_stream<VDisortSettings>::read(std::istream& is_xml, VDisortSettings& v, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  xml_read_from_stream(is_xml, v.quadrature_dimension, pbifs);
  xml_read_from_stream(is_xml, v.legendre_polynomial_dimension, pbifs);
  xml_read_from_stream(is_xml, v.fourier_mode_dimension, pbifs);
  xml_read_from_stream(is_xml, v.freq_grid, pbifs);
  xml_read_from_stream(is_xml, v.alt_grid, pbifs);
  xml_read_from_stream(is_xml, v.solar_azimuth_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_zenith_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_source, pbifs);
  xml_read_from_stream(is_xml, v.surface_brdf, pbifs);
  xml_read_from_stream(is_xml, v.surface_beam_brdf, pbifs);
  xml_read_from_stream(is_xml, v.surface_beam_mu, pbifs);
  xml_read_from_stream(is_xml, v.surface_specular, pbifs);
  xml_read_from_stream(is_xml, v.optical_thicknesses, pbifs);
  xml_read_from_stream(is_xml, v.single_scattering_albedo, pbifs);
  xml_read_from_stream(is_xml, v.source_polynomial, pbifs);
  xml_read_from_stream(is_xml, v.legendre_coefficients, pbifs);
  xml_read_from_stream(is_xml, v.upward_boundary_condition, pbifs);
  xml_read_from_stream(is_xml, v.downward_boundary_condition, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
}

void xml_io_stream<VDisortSettings>::write(std::ostream&          os_xml,
                                           const VDisortSettings& v,
                                           bofstream*             pbofs,
                                           std::string_view       name) {
  XMLTag open_tag(type_name, "name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, v.quadrature_dimension, pbofs, "quadrature_dimension");
  xml_write_to_stream(os_xml, v.legendre_polynomial_dimension, pbofs, "legendre_polynomial_dimension");
  xml_write_to_stream(os_xml, v.fourier_mode_dimension, pbofs, "fourier_mode_dimension");
  xml_write_to_stream(os_xml, v.freq_grid, pbofs, "freq_grid");
  xml_write_to_stream(os_xml, v.alt_grid, pbofs, "alt_grid");
  xml_write_to_stream(os_xml, v.solar_azimuth_angle, pbofs, "solar_azimuth_angle");
  xml_write_to_stream(os_xml, v.solar_zenith_angle, pbofs, "solar_zenith_angle");
  xml_write_to_stream(os_xml, v.solar_source, pbofs, "solar_source");
  xml_write_to_stream(os_xml, v.surface_brdf, pbofs, "surface_brdf");
  xml_write_to_stream(os_xml, v.surface_beam_brdf, pbofs, "surface_beam_brdf");
  xml_write_to_stream(os_xml, v.surface_beam_mu, pbofs, "surface_beam_mu");
  xml_write_to_stream(os_xml, v.surface_specular, pbofs, "surface_specular");
  xml_write_to_stream(os_xml, v.optical_thicknesses, pbofs, "optical_thicknesses");
  xml_write_to_stream(os_xml, v.single_scattering_albedo, pbofs, "single_scattering_albedo");
  xml_write_to_stream(os_xml, v.source_polynomial, pbofs, "source_polynomial");
  xml_write_to_stream(os_xml, v.legendre_coefficients, pbofs, "legendre_coefficients");
  xml_write_to_stream(os_xml, v.upward_boundary_condition, pbofs, "upward_boundary_condition");
  xml_write_to_stream(os_xml, v.downward_boundary_condition, pbofs, "downward_boundary_condition");

  open_tag.write_to_end_stream(os_xml);
}
