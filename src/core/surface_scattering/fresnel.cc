#include "fresnel.h"

#include <physics_funcs.h>
#include <xml_io_base.h>

#include "rtepack.h"

namespace {

//! Fold an ARTS zenith angle in [0, 180] onto [0, 90] — the domain where the
//! Fresnel formulas are physical.  R(180 - za) = R(za) for real indices, so
//! the fold is exact for the downward hemisphere.
Numeric za_fold(Numeric za) { return za <= 90.0 ? za : 180.0 - za; }

//! Compute the Fresnel (specular) BRDF and complementary emissivity from a
//! per-frequency real refractive-index pair (n1 constant, n2 varying).
//!
//! Uses the *pair* overload of fresnel() (physics_funcs.h).  Beyond the
//! critical angle its complex sqrt makes both amplitudes unimodular, so total
//! internal reflection gives a reflectance of 1 only up to a few ulp -- never
//! bit-exactly, and with a libm-dependent residual.  Do not test it for exact
//! equality.
//!
//! The BRDF matrix at a frequency depends only on the incoming zenith angle
//! (degrees); the emissivity tensor (no incidence axis) is the energy
//! complement I4 - M evaluated at each scattering zenith angle.
void fresnel_tensors(const Numeric   n1,
                     const auto&     n2_data,
                     const Vector&   za_inc_grid,
                     const Vector&   za_scat_grid,
                     Index           nf,
                     Index           nzi,
                     Index           nai,
                     Index           nzs,
                     Index           nas,
                     MuelmatTensor5& brdf_specular,
                     MuelmatTensor3& emissivity_specular) {
  for (Index f = 0; f < nf; ++f) {
    const Numeric n2 = n2_data[f];

    for (Index zi = 0; zi < nzi; ++zi) {
      // ARTS zenith angles span [0, 180]; the Fresnel formulas are defined for
      // the angle from the surface normal on the illuminated side.  Fold the
      // downward hemisphere onto the upward one (R(180-za) = R(za)); passing
      // za > 90 raw inverts the reflectance (R -> 1/R) and blows up near the
      // conjugate Brewster angle.
      const Numeric za_eff = za_fold(za_inc_grid[zi]);
      const auto [Rv, Rh]  = fresnel(Complex{n1}, Complex{n2}, za_eff);
      const Muelmat M      = rtepack::fresnel_reflectance(Rv, Rh);

      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as) brdf_specular[f, zi, ai, zs, as] = M;
    }

    for (Index zs = 0; zs < nzs; ++zs) {
      for (Index as = 0; as < nas; ++as) {
        const Numeric za_eff = za_fold(za_scat_grid[zs]);
        const auto [Rv, Rh]  = fresnel(Complex{n1}, Complex{n2}, za_eff);
        const Muelmat M      = rtepack::fresnel_reflectance(Rv, Rh);
        Muelmat       E{};
        for (Index r = 0; r < 4; ++r)
          for (Index c = 0; c < 4; ++c) E[r, c] = (r == c ? 1.0 : 0.0) - M[r, c];
        emissivity_specular[f, zs, as] = E;
      }
    }
  }
}

}  // namespace

using namespace surface_scattering;

void FresnelSurfaceScatterer::check_n1() const {
  if (not(n1 > 0.0)) arts::user_error("FresnelSurfaceScatterer n1 must be positive, but is {}.", n1);
}

FresnelSurfaceScatterer::FresnelSurfaceScatterer(Data refractive_index_, Numeric n1_)
    : refractive_index(std::move(refractive_index_)), n1{n1_} {
  check_n1();
}

SurfaceScatteringModelProperties FresnelSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric       lat,
    Numeric       lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  check_n1();

  const Vector n2_data = refractive_index.at(f_grid, lat, lon);

  const Size nf_size = f_grid.size();
  for (Index f = 0; f < Index(nf_size); ++f) {
    if (not(n2_data[f] > 0.0))
      arts::user_error("Refractive index must be positive, but {} (f={}).", n2_data[f], f_grid[f]);
  }

  MuelmatTensor5 brdf_specular(
      nf_size, za_inc_grid.size(), aa_inc_grid.size(), za_scat_grid.size(), aa_scat_grid.size(), rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_specular(f_grid.size(), za_scat_grid.size(), aa_scat_grid.size(), rtepack::muelmat{0.0});

  fresnel_tensors(n1,
                  n2_data,
                  za_inc_grid,
                  za_scat_grid,
                  f_grid.size(),
                  za_inc_grid.size(),
                  aa_inc_grid.size(),
                  za_scat_grid.size(),
                  aa_scat_grid.size(),
                  brdf_specular,
                  emissivity_specular);

  // Purely specular model: the diffuse components stay missing (std::nullopt)
  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse        = std::nullopt,
      .emissivity_vector_diffuse  = std::nullopt,
      .brdf_matrix_specular       = std::move(brdf_specular),
      .emissivity_vector_specular = std::move(emissivity_specular),
  };
}

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const FresnelSurfaceScatterer& s) {
  return os << std::format("FresnelSurfaceScatterer n1={}", s.n1);
}

void xml_io_stream<FresnelSurfaceScatterer>::write(std::ostream&                  os,
                                                   const FresnelSurfaceScatterer& x,
                                                   bofstream*                     pbofs,
                                                   std::string_view               name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.n1, pbofs, "n1"sv);
  xml_write_to_stream(os, x.refractive_index, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FresnelSurfaceScatterer>::read(std::istream& is, FresnelSurfaceScatterer& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.n1, pbifs);
  xml_read_from_stream(is, x.refractive_index, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
