#include "vdisort_settings.h"

#include <arts_conversions.h>
#include <debug.h>

#include <algorithm>
#include <cmath>
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
  bidirectional_reflectance_distribution_functions.resize(nfreq, 0);
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
      (bidirectional_reflectance_distribution_functions.shape() !=
       std::array{nfreq, bidirectional_reflectance_distribution_functions.ncols()}) or
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
                            bidirectional_reflectance_distribution_functions.ncols());
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

  for (Index i = 0; i < bidirectional_reflectance_distribution_functions.ncols(); i++) {
    dis.brdf_modes()[i] = bidirectional_reflectance_distribution_functions[iv, i];
  }

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

void xml_io_stream<VDisortBDRF>::read(std::istream& is, VDisortBDRF& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.cosine, pbifs);
  xml_read_from_stream(is, x.sine, pbifs);
  xml_read_from_stream(is, x.beam_cosine, pbifs);
  xml_read_from_stream(is, x.beam_sine, pbifs);
  xml_read_from_stream(is, x.specular, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<VDisortBDRF>::write(std::ostream& os, const VDisortBDRF& x, bofstream* pbofs, std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.cosine, pbofs);
  xml_write_to_stream(os, x.sine, pbofs);
  xml_write_to_stream(os, x.beam_cosine, pbofs);
  xml_write_to_stream(os, x.beam_sine, pbofs);
  xml_write_to_stream(os, x.specular, pbofs);

  tag.write_to_end_stream(os);
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
  xml_read_from_stream(is_xml, v.bidirectional_reflectance_distribution_functions, pbifs);
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
  xml_write_to_stream(os_xml,
                      v.bidirectional_reflectance_distribution_functions,
                      pbofs,
                      "bidirectional_reflectance_distribution_functions");
  xml_write_to_stream(os_xml, v.optical_thicknesses, pbofs, "optical_thicknesses");
  xml_write_to_stream(os_xml, v.single_scattering_albedo, pbofs, "single_scattering_albedo");
  xml_write_to_stream(os_xml, v.source_polynomial, pbofs, "source_polynomial");
  xml_write_to_stream(os_xml, v.legendre_coefficients, pbofs, "legendre_coefficients");
  xml_write_to_stream(os_xml, v.upward_boundary_condition, pbofs, "upward_boundary_condition");
  xml_write_to_stream(os_xml, v.downward_boundary_condition, pbofs, "downward_boundary_condition");

  open_tag.write_to_end_stream(os_xml);
}
