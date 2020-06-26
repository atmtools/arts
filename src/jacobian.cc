/* Copyright (C) 2004-2012 Mattias Ekstrom  <ekstrom@rss.chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  @file   jacobian.cc
  @author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  @date   2004-09-14

  @brief  Routines for setting up the jacobian.
 */

#include "jacobian.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "global_data.h"
#include "lin_alg.h"
#include "physics_funcs.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric NAT_LOG_TEN;
extern const Numeric PI;

extern const String PROPMAT_SUBSUBTAG;

ostream& operator<<(ostream& os, const RetrievalQuantity& ot) {
  return os << "\n       Main tag = " << ot.MainTag()
            << "\n       Sub  tag = " << ot.Subtag()
            << "\n           Mode = " << ot.Mode();
}

void jac_ranges_indices(ArrayOfArrayOfIndex& jis,
                        bool& any_affine,
                        const ArrayOfRetrievalQuantity& jqs,
                        const bool& before_affine) {
  jis.resize(jqs.nelem());

  any_affine = false;

  // Indices before affine transformation
  if (before_affine) {
    for (Index i = 0; i < jqs.nelem(); ++i) {
      jis[i] = ArrayOfIndex(2);
      if (i > 0) {
        jis[i][0] = jis[i - 1][1] + 1;
      } else {
        jis[i][0] = 0;
      }
      const RetrievalQuantity& jq = jqs[i];
      jis[i][1] = jis[i][0] + jq.nelem() - 1;
      if (jq.HasAffine()) {
        any_affine = true;
      }
    }
  }

  // After affine transformation
  else {
    for (Index i = 0; i < jqs.nelem(); ++i) {
      jis[i] = ArrayOfIndex(2);
      if (i > 0) {
        jis[i][0] = jis[i - 1][1] + 1;
      } else {
        jis[i][0] = 0;
      }
      const RetrievalQuantity& jq = jqs[i];
      if (jq.HasAffine()) {
        jis[i][1] = jis[i][0] + jq.TransformationMatrix().ncols() - 1;
        any_affine = true;
      } else {
        jis[i][1] = jis[i][0] + jq.nelem() - 1;
      }
    }
  }
}

void transform_jacobian(Matrix& jacobian,
                        const Vector x,
                        const ArrayOfRetrievalQuantity& jqs) {
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices(jis, any_affine, jqs, true);

  Vector x_t(x);
  transform_x_back(x_t, jqs, false);

  // Apply functional transformations
  for (Index i = 0; i < jqs.nelem(); ++i) {
    const RetrievalQuantity& jq = jqs[i];
    const String tfun = jq.TransformationFunc();
    // Remember to add new functions also to transform_jacobian and transform_x_back
    if (tfun == "") {
      // Nothing to do
    } else if (tfun == "log") {
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker, c) *= exp(x_t[c]);
      }
    } else if (tfun == "log10") {
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker, c) *= NAT_LOG_TEN * pow(10.0, x_t[c]);
      }
    } else if (tfun == "atanh") {
      const Vector& pars = jq.TFuncParameters();
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker, c) *=
            2 * (pars[1] - pars[0]) / pow(exp(-x_t[c]) + exp(x_t[c]), 2.0);
      }
    } else {
      assert(0);
    }
  }

  // Apply affine transformations
  if (any_affine) {
    ArrayOfArrayOfIndex jis_t;
    jac_ranges_indices(jis_t, any_affine, jqs);

    Matrix jacobian_t(jacobian.nrows(), jis_t.back()[1] + 1);

    for (Index i = 0; i < jqs.nelem(); ++i) {
      const RetrievalQuantity& jq = jqs[i];
      Index col_start = jis[i][0];
      Index col_extent = jis[i][1] - jis[i][0] + 1;
      Range col_range(col_start, col_extent);
      Index col_start_t = jis_t[i][0];
      Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
      Range col_range_t(col_start_t, col_extent_t);
      if (jq.HasAffine()) {
        mult(jacobian_t(joker, col_range_t),
             jacobian(joker, col_range),
             jq.TransformationMatrix());
      } else {
        jacobian_t(joker, col_range_t) = jacobian(joker, col_range);
      }
    }
    swap(jacobian_t, jacobian);
  }
}

void transform_x(Vector& x, const ArrayOfRetrievalQuantity& jqs) {
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices(jis, any_affine, jqs, true);

  // Apply functional transformations
  for (Index i = 0; i < jqs.nelem(); ++i) {
    const RetrievalQuantity& jq = jqs[i];
    const String tfun = jq.TransformationFunc();
    // Remember to add new functions also to transform_jacobian and transform_x_back
    if (tfun == "") {
      // Nothing to do
    } else if (tfun == "log") {
      const Vector& pars = jq.TFuncParameters();
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if (x[r] <= pars[0]) {
          ostringstream os;
          os << "log-transformation selected for retrieval quantity with\n"
             << "index " << i << " (0-based), but at least one value <= z_min\n"
             << "found for this quantity. This is not allowed.";
          throw std::runtime_error(os.str());
        }
        x[r] = log(x[r] - pars[0]);
      }
    } else if (tfun == "log10") {
      const Vector& pars = jq.TFuncParameters();
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if (x[r] <= 0) {
          ostringstream os;
          os << "log10-transformation selected for retrieval quantity with\n"
             << "index " << i << " (0-based), but at least one value <= z_min\n"
             << "found for this quantity. This is not allowed.";
          throw std::runtime_error(os.str());
        }
        x[r] = log10(x[r] - pars[0]);
      }
    } else if (tfun == "atanh") {
      const Vector& pars = jq.TFuncParameters();
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if (x[r] <= pars[0]) {
          ostringstream os;
          os << "atanh-transformation selected for retrieval quantity with\n"
             << "index " << i << " (0-based), but at least one value <= z_min\n"
             << "found for this quantity. This is not allowed.";
          throw std::runtime_error(os.str());
        }
        if (x[r] >= pars[1]) {
          ostringstream os;
          os << "atanh-transformation selected for retrieval quantity with\n"
             << "index " << i << " (0-based), but at least one value is\n"
             << ">= z_max. This is not allowed.";
          throw std::runtime_error(os.str());
        }
        x[r] = atanh(2 * (x[r] - pars[0]) / (pars[1] - pars[0]) - 1);
      }
    } else {
      assert(0);
    }
  }

  // Apply affine transformations
  if (any_affine) {
    ArrayOfArrayOfIndex jis_t;
    jac_ranges_indices(jis_t, any_affine, jqs);

    Vector x_t(jis_t.back()[1] + 1);

    for (Index i = 0; i < jqs.nelem(); ++i) {
      const RetrievalQuantity& jq = jqs[i];
      Index col_start = jis[i][0];
      Index col_extent = jis[i][1] - jis[i][0] + 1;
      Range col_range(col_start, col_extent);
      Index col_start_t = jis_t[i][0];
      Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
      Range col_range_t(col_start_t, col_extent_t);
      if (jq.HasAffine()) {
        Vector t(x[col_range]);
        t -= jq.OffsetVector();
        mult(x_t[col_range_t], transpose(jq.TransformationMatrix()), t);
      } else {
        x_t[col_range_t] = x[col_range];
      }
    }
    swap(x, x_t);
  }
}

void transform_x_back(Vector& x_t,
                      const ArrayOfRetrievalQuantity& jqs,
                      bool revert_functional_transforms) {
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices(jis, any_affine, jqs, true);

  // Revert affine transformations
  // Apply affine transformations
  if (any_affine) {
    ArrayOfArrayOfIndex jis_t;
    jac_ranges_indices(jis_t, any_affine, jqs);

    Vector x(jis.back()[1] + 1);

    for (Index i = 0; i < jqs.nelem(); ++i) {
      const RetrievalQuantity& jq = jqs[i];
      Index col_start = jis[i][0];
      Index col_extent = jis[i][1] - jis[i][0] + 1;
      Range col_range(col_start, col_extent);
      Index col_start_t = jis_t[i][0];
      Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
      Range col_range_t(col_start_t, col_extent_t);
      if (jq.HasAffine()) {
        mult(x[col_range], jq.TransformationMatrix(), x_t[col_range_t]);
        x[col_range] += jq.OffsetVector();
      } else {
        x[col_range] = x_t[col_range_t];
      }
    }
    swap(x_t, x);
  }

  if (revert_functional_transforms) {
    // Revert functional transformations
    for (Index i = 0; i < jqs.nelem(); ++i) {
      const RetrievalQuantity& jq = jqs[i];
      const String tfun = jq.TransformationFunc();
      // Remember to add new functions also to transform_jacobian and transform_x_back
      if (tfun == "") {
        // Nothing to do
      } else if (tfun == "log") {
        const Vector& pars = jq.TFuncParameters();
        for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
          x_t[r] = pars[0] + exp(x_t[r]);
        }
      } else if (tfun == "log10") {
        const Vector& pars = jq.TFuncParameters();
        for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
          x_t[r] = pars[0] + pow(10.0, x_t[r]);
        }
      } else if (tfun == "atanh") {
        const Vector& pars = jq.TFuncParameters();
        for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
          x_t[r] = pars[0] + ((pars[1] - pars[0]) / 2) * (1 + tanh(x_t[r]));
        }
      } else {
        assert(0);
      }
    }
  }
}

/*===========================================================================
  === Help sub-functions to handle analytical jacobians (in alphabetical order)
  ===========================================================================*/

// Small help function, to make the code below cleaner
void from_dpath_to_dx(MatrixView diy_dx,
                      ConstMatrixView diy_dq,
                      const Numeric& w) {
  for (Index irow = 0; irow < diy_dx.nrows(); irow++) {
    for (Index icol = 0; icol < diy_dx.ncols(); icol++) {
      diy_dx(irow, icol) += w * diy_dq(irow, icol);
    }
  }
}

void diy_from_path_to_rgrids(Tensor3View diy_dx,
                             const RetrievalQuantity& jacobian_quantity,
                             ConstTensor3View diy_dpath,
                             const Index& atmosphere_dim,
                             const Ppath& ppath,
                             ConstVectorView ppath_p) {
  assert(jacobian_quantity.Grids().nelem() == atmosphere_dim or jacobian_quantity.Grids().empty());

  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric extpolfac = 1.0e99;

  if (ppath.np > 1)  // Otherwise nothing to do here
  {
    // Pressure
    Index nr1 = jacobian_quantity.Grids().empty() ? 0 : jacobian_quantity.Grids()[0].nelem();
    ArrayOfGridPos gp_p(ppath.np);
    if (nr1 > 1) {
      p2gridpos(gp_p, jacobian_quantity.Grids()[0], ppath_p, extpolfac);
      jacobian_type_extrapol(gp_p);
    } else {
      gp4length1grid(gp_p);
    }

    // Latitude
    Index nr2 = 1;
    ArrayOfGridPos gp_lat;
    if (atmosphere_dim > 1) {
      gp_lat.resize(ppath.np);
      nr2 = jacobian_quantity.Grids().empty() ? 0 : jacobian_quantity.Grids()[1].nelem();
      if (nr2 > 1) {
        gridpos(gp_lat,
                jacobian_quantity.Grids()[1],
                ppath.pos(joker, 1),
                extpolfac);
        jacobian_type_extrapol(gp_lat);
      } else {
        gp4length1grid(gp_lat);
      }
    }

    // Longitude
    ArrayOfGridPos gp_lon;
    if (atmosphere_dim > 2) {
      Index nr3 = jacobian_quantity.Grids().empty() ? 0 : jacobian_quantity.Grids()[2].nelem();
      gp_lon.resize(ppath.np);
      if (nr3 > 1) {
        gridpos(gp_lon,
                jacobian_quantity.Grids()[2],
                ppath.pos(joker, 2),
                extpolfac);
        jacobian_type_extrapol(gp_lon);
      } else {
        gp4length1grid(gp_lon);
      }
    }

    //- 1D
    if (atmosphere_dim == 1) {
      for (Index ip = 0; ip < ppath.np; ip++) {
        if (gp_p[ip].fd[1] > 0) {
          from_dpath_to_dx(diy_dx(gp_p[ip].idx, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_p[ip].fd[1]);
        }
        if (gp_p[ip].fd[0] > 0) {
          from_dpath_to_dx(diy_dx(gp_p[ip].idx + 1, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_p[ip].fd[0]);
        }
      }
    }

    //- 2D
    else if (atmosphere_dim == 2) {
      for (Index ip = 0; ip < ppath.np; ip++) {
        Index ix = nr1 * gp_lat[ip].idx + gp_p[ip].idx;
        // Low lat, low p
        if (gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(diy_dx(ix, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_lat[ip].fd[1] * gp_p[ip].fd[1]);
        // Low lat, high p
        if (gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(diy_dx(ix + 1, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_lat[ip].fd[1] * gp_p[ip].fd[0]);
        // High lat, low p
        if (gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(diy_dx(ix + nr1, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_lat[ip].fd[0] * gp_p[ip].fd[1]);
        // High lat, high p
        if (gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(diy_dx(ix + nr1 + 1, joker, joker),
                           diy_dpath(ip, joker, joker),
                           gp_lat[ip].fd[0] * gp_p[ip].fd[0]);
      }
    }

    //- 3D
    else if (atmosphere_dim == 3) {
      for (Index ip = 0; ip < ppath.np; ip++) {
        Index ix =
            nr2 * nr1 * gp_lon[ip].idx + nr1 * gp_lat[ip].idx + gp_p[ip].idx;
        // Low lon, low lat, low p
        if (gp_lon[ip].fd[1] > 0 && gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(
              diy_dx(ix, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[1] * gp_lat[ip].fd[1] * gp_p[ip].fd[1]);
        // Low lon, low lat, high p
        if (gp_lon[ip].fd[1] > 0 && gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(
              diy_dx(ix + 1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[1] * gp_lat[ip].fd[1] * gp_p[ip].fd[0]);
        // Low lon, high lat, low p
        if (gp_lon[ip].fd[1] > 0 && gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(
              diy_dx(ix + nr1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[1] * gp_lat[ip].fd[0] * gp_p[ip].fd[1]);
        // Low lon, high lat, high p
        if (gp_lon[ip].fd[1] > 0 && gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(
              diy_dx(ix + nr1 + 1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[1] * gp_lat[ip].fd[0] * gp_p[ip].fd[0]);

        // Increase *ix* (to be valid for high lon level)
        ix += nr2 * nr1;

        // High lon, low lat, low p
        if (gp_lon[ip].fd[0] > 0 && gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(
              diy_dx(ix, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[0] * gp_lat[ip].fd[1] * gp_p[ip].fd[1]);
        // High lon, low lat, high p
        if (gp_lon[ip].fd[0] > 0 && gp_lat[ip].fd[1] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(
              diy_dx(ix + 1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[0] * gp_lat[ip].fd[1] * gp_p[ip].fd[0]);
        // High lon, high lat, low p
        if (gp_lon[ip].fd[0] > 0 && gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[1] > 0)
          from_dpath_to_dx(
              diy_dx(ix + nr1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[0] * gp_lat[ip].fd[0] * gp_p[ip].fd[1]);
        // High lon, high lat, high p
        if (gp_lon[ip].fd[0] > 0 && gp_lat[ip].fd[0] > 0 && gp_p[ip].fd[0] > 0)
          from_dpath_to_dx(
              diy_dx(ix + nr1 + 1, joker, joker),
              diy_dpath(ip, joker, joker),
              gp_lon[ip].fd[0] * gp_lat[ip].fd[0] * gp_p[ip].fd[0]);
      }
    }
  }
}

void diy_from_pos_to_rgrids(Tensor3View diy_dx,
                            const RetrievalQuantity& jacobian_quantity,
                            ConstMatrixView diy_dpos,
                            const Index& atmosphere_dim,
                            ConstVectorView rtp_pos) {
  assert(jacobian_quantity.Grids().nelem() ==
         max(atmosphere_dim - 1, Index(1)));
  assert(rtp_pos.nelem() == atmosphere_dim);

  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric extpolfac = 1.0e99;

  // Handle 1D separately
  if (atmosphere_dim == 1) {
    diy_dx(0, joker, joker) = diy_dpos;
    return;
  }

  // Latitude
  Index nr1 = 1;
  ArrayOfGridPos gp_lat;
  {
    gp_lat.resize(1);
    nr1 = jacobian_quantity.Grids()[0].nelem();
    if (nr1 > 1) {
      gridpos(gp_lat,
              jacobian_quantity.Grids()[0],
              Vector(1, rtp_pos[1]),
              extpolfac);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
  }

  // Longitude
  ArrayOfGridPos gp_lon;
  if (atmosphere_dim > 2) {
    gp_lon.resize(1);
    if (jacobian_quantity.Grids()[1].nelem() > 1) {
      gridpos(gp_lon,
              jacobian_quantity.Grids()[1],
              Vector(1, rtp_pos[2]),
              extpolfac);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
  }

  //- 2D
  if (atmosphere_dim == 2) {
    if (gp_lat[0].fd[1] > 0) {
      from_dpath_to_dx(diy_dx(gp_lat[0].idx, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lat[0].fd[1]);
    }
    if (gp_lat[0].fd[0] > 0) {
      from_dpath_to_dx(diy_dx(gp_lat[0].idx + 1, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lat[0].fd[0]);
    }
  }
  //- 3D
  else {
    Index ix = nr1 * gp_lon[0].idx + gp_lat[0].idx;
    // Low lon, low lat
    if (gp_lon[0].fd[1] > 0 && gp_lat[0].fd[1] > 0)
      from_dpath_to_dx(diy_dx(ix, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lon[0].fd[1] * gp_lat[0].fd[1]);
    // Low lon, high lat
    if (gp_lon[0].fd[1] > 0 && gp_lat[0].fd[0] > 0)
      from_dpath_to_dx(diy_dx(ix + 1, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lon[0].fd[1] * gp_lat[0].fd[0]);
    // High lon, low lat
    if (gp_lon[0].fd[0] > 0 && gp_lat[0].fd[1] > 0)
      from_dpath_to_dx(diy_dx(ix + nr1, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lon[0].fd[0] * gp_lat[0].fd[1]);
    // High lon, high lat
    if (gp_lon[0].fd[0] > 0 && gp_lat[0].fd[0] > 0)
      from_dpath_to_dx(diy_dx(ix + nr1 + 1, joker, joker),
                       diy_dpos(joker, joker),
                       gp_lon[0].fd[0] * gp_lat[0].fd[0]);
  }
}

void get_pointers_for_analytical_jacobians(
    ArrayOfIndex& abs_species_i,
    ArrayOfIndex& scat_species_i,
    ArrayOfIndex& is_t,
    ArrayOfIndex& wind_i,
    ArrayOfIndex& magfield_i,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Index& cloudbox_on,
    const ArrayOfString& scat_species) {
  FOR_ANALYTICAL_JACOBIANS_DO(
      //
      if (jacobian_quantities[iq] == Jacobian::Atm::Temperature &&
          jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG) {
        is_t[iq] = Index(JacobianType::Temperature);
      } else { is_t[iq] = Index(JacobianType::None); }
      //
      if (jacobian_quantities[iq].Target().isSpeciesVMR()) {
        if (jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG) {
          bool test_available = false;
          for (Index ii = 0; ii < abs_species.nelem(); ii++) {
            if (abs_species[ii][0].Species() ==
                SpeciesTag(jacobian_quantities[iq].Subtag()).Species()) {
              test_available = true;
              abs_species_i[iq] = ii;
              break;
            }
          }
          if (!test_available) {
            ostringstream os;
            os << "Could not find " << jacobian_quantities[iq].Subtag()
               << " in species of abs_species.\n";
            throw std::runtime_error(os.str());
          }
        } else {
          ArrayOfSpeciesTag atag;
          array_species_tag_from_string(atag, jacobian_quantities[iq].Subtag());
          abs_species_i[iq] = chk_contains("abs_species", abs_species, atag);
        }
      } else if (jacobian_quantities[iq] == Jacobian::Atm::Particulates ||
                 jacobian_quantities[iq] == Jacobian::Atm::Electrons) {
        abs_species_i[iq] = -9999;
      } else { abs_species_i[iq] = -1; }
      //
      if (cloudbox_on &&
          jacobian_quantities[iq] == Jacobian::Special::ScatteringString) {
        scat_species_i[iq] =
            find_first(scat_species, jacobian_quantities[iq].Subtag());
        if (scat_species_i[iq] < 0) {
          ostringstream os;
          os << "Jacobian quantity with index " << iq << " refers to\n"
             << "  " << jacobian_quantities[iq].Subtag()
             << "\nbut this species could not be found in *scat_species*.";
          throw runtime_error(os.str());
        }
      } else { scat_species_i[iq] = -1; }
      //
      if (jacobian_quantities[iq].Target().isWind() &&
          jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG) {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index(c) - 116;
        if (test == 1)
          wind_i[iq] = Index(JacobianType::WindFieldU);
        else if (test == 2)
          wind_i[iq] = Index(JacobianType::WindFieldV);
        else if (test == 3)
          wind_i[iq] = Index(JacobianType::WindFieldW);
        else if (test == (Index('s') - 116))
          wind_i[iq] = Index(JacobianType::AbsWind);
      } else { wind_i[iq] = Index(JacobianType::None); }
      //
      if (jacobian_quantities[iq].Target().isMagnetic() &&
          jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG) {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index(c) - 116;
        if (test == 1)
          magfield_i[iq] = Index(JacobianType::MagFieldU);
        else if (test == 2)
          magfield_i[iq] = Index(JacobianType::MagFieldV);
        else if (test == 3)
          magfield_i[iq] = Index(JacobianType::MagFieldW);
        else if (test == (Index('s') - 116))
          magfield_i[iq] = Index(JacobianType::AbsMag);
      } else { magfield_i[iq] = Index(JacobianType::None); }

  )
}

/*===========================================================================
  === Other functions, in alphabetical order
  ===========================================================================*/

bool check_retrieval_grids(ArrayOfVector& grids,
                           ostringstream& os,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Vector& p_retr,
                           const Vector& lat_retr,
                           const Vector& lon_retr,
                           const String& p_retr_name,
                           const String& lat_retr_name,
                           const String& lon_retr_name,
                           const Index& dim) {
  assert(grids.nelem() == dim);

  if (p_retr.nelem() == 0) {
    os << "The grid vector *" << p_retr_name << "* is empty,"
       << " at least one pressure level\n"
       << "should be specified.";
    return false;
  } else if (!is_decreasing(p_retr)) {
    os << "The pressure grid vector *" << p_retr_name << "* is not a\n"
       << "strictly decreasing vector, which is required.";
    return false;
  } else if (p_grid.nelem() == 1 and p_grid.nelem() == p_retr.nelem()) {
    if (p_grid[0] not_eq p_retr[0]) {
      os << "Mismatching 1-long grids for " << p_retr_name;
      return false;
    }

    // Necessary repeat but grids are OK
    grids[0] = p_retr;
  } else if (log(p_retr[0]) > 1.5 * log(p_grid[0]) - 0.5 * log(p_grid[1]) ||
             log(p_retr[p_retr.nelem() - 1]) <
                 1.5 * log(p_grid[p_grid.nelem() - 1]) -
                     0.5 * log(p_grid[p_grid.nelem() - 2])) {
    os << "The grid vector *" << p_retr_name << "* is not covered by the\n"
       << "corresponding atmospheric grid.";
    return false;
  } else {
    // Pressure grid ok, add it to grids
    grids[0] = p_retr;
  }

  if (dim >= 2) {
    // If 2D and 3D atmosphere, check latitude grid
    if (lat_retr.nelem() == 0) {
      os << "The grid vector *" << lat_retr_name << "* is empty,"
         << " at least one latitude\n"
         << "should be specified for a 2D/3D atmosphere.";
      return false;
    } else if (!is_increasing(lat_retr)) {
      os << "The latitude grid vector *" << lat_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;
    } else if (lat_grid.nelem() == 1 and lat_grid.nelem() == lat_retr.nelem()) {
      if (lat_grid[0] not_eq lat_retr[0]) {
        os << "Mismatching 1-long grids for " << lat_retr_name;
        return false;
      }

      // Necessary repeat but grids are OK
      grids[1] = lat_retr;
    } else if (lat_retr[0] < 1.5 * lat_grid[0] - 0.5 * lat_grid[1] ||
               lat_retr[lat_retr.nelem() - 1] >
                   1.5 * lat_grid[lat_grid.nelem() - 1] -
                       0.5 * lat_grid[lat_grid.nelem() - 2]) {
      os << "The grid vector *" << lat_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    } else {
      // Latitude grid ok, add it to grids
      grids[1] = lat_retr;
    }
    if (dim == 3) {
      // For 3D atmosphere check longitude grid
      if (lon_retr.nelem() == 0) {
        os << "The grid vector *" << lon_retr_name << "* is empty,"
           << " at least one longitude\n"
           << "should be specified for a 3D atmosphere.";
        return false;
      } else if (!is_increasing(lon_retr)) {
        os << "The longitude grid vector *" << lon_retr_name << "* is not a\n"
           << "strictly increasing vector, which is required.";
        return false;
      } else if (lon_grid.nelem() == 1 and
                 lon_grid.nelem() == lon_retr.nelem()) {
        if (lon_grid[0] not_eq lon_retr[0]) {
          os << "Mismatching 1-long grids for " << lon_retr_name;
          return false;
        }

        // Necessary repeat but grids are OK
        grids[2] = lon_retr;
      } else if (lon_retr[0] < 1.5 * lon_grid[0] - 0.5 * lon_grid[1] ||
                 lon_retr[lon_retr.nelem() - 1] >
                     1.5 * lon_grid[lon_grid.nelem() - 1] -
                         0.5 * lon_grid[lon_grid.nelem() - 2]) {
        os << "The grid vector *" << lon_retr_name
           << "* is not covered by the\n"
           << "corresponding atmospheric grid.";
        return false;
      } else {
        // Longitude grid ok, add it to grids
        grids[2] = lon_retr;
      }
    }
  }
  return true;
}

bool check_retrieval_grids(ArrayOfVector& grids,
                           ostringstream& os,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Vector& lat_retr,
                           const Vector& lon_retr,
                           const String& lat_retr_name,
                           const String& lon_retr_name,
                           const Index& dim) {
  assert(grids.nelem() == max(dim - 1, Index(1)));

  if (dim == 1) {
    // Here we only need to create a length 1 dummy grid
    grids[0].resize(1);
    grids[0][0] = 0;
  }

  if (dim >= 2) {
    // If 2D and 3D atmosphere, check latitude grid
    if (lat_retr.nelem() == 0) {
      os << "The grid vector *" << lat_retr_name << "* is empty,"
         << " at least one latitude\n"
         << "should be specified for a 2D/3D atmosphere.";
      return false;
    } else if (!is_increasing(lat_retr)) {
      os << "The latitude grid vector *" << lat_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;
    } else if (lat_grid.nelem() == 1 and lat_grid.nelem() == lat_retr.nelem()) {
      if (lat_grid[0] not_eq lat_retr[0]) {
        os << "Mismatching 1-long grids for " << lat_retr_name;
        return false;
      }

      // Necessary repeat but grids are OK
      grids[0] = lat_retr;
    } else if (lat_retr[0] < 1.5 * lat_grid[0] - 0.5 * lat_grid[1] ||
               lat_retr[lat_retr.nelem() - 1] >
                   1.5 * lat_grid[lat_grid.nelem() - 1] -
                       0.5 * lat_grid[lat_grid.nelem() - 2]) {
      os << "The grid vector *" << lat_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    } else {
      // Latitude grid ok, add it to grids
      grids[0] = lat_retr;
    }

    if (dim == 3) {
      // For 3D atmosphere check longitude grid
      if (lon_retr.nelem() == 0) {
        os << "The grid vector *" << lon_retr_name << "* is empty,"
           << " at least one longitude\n"
           << "should be specified for a 3D atmosphere.";
        return false;
      } else if (!is_increasing(lon_retr)) {
        os << "The longitude grid vector *" << lon_retr_name << "* is not a\n"
           << "strictly increasing vector, which is required.";
        return false;
      } else if (lon_grid.nelem() == 1 and
                 lon_grid.nelem() == lon_retr.nelem()) {
        if (lon_grid[0] not_eq lon_retr[0]) {
          os << "Mismatching 1-long grids for " << lon_retr_name;
          return false;
        }

        // Necessary repeat but grids are OK
        grids[1] = lon_retr;
      } else if (lon_retr[0] < 1.5 * lon_grid[0] - 0.5 * lon_grid[1] ||
                 lon_retr[lon_retr.nelem() - 1] >
                     1.5 * lon_grid[lon_grid.nelem() - 1] -
                         0.5 * lon_grid[lon_grid.nelem() - 2]) {
        os << "The grid vector *" << lon_retr_name
           << "* is not covered by the\n"
           << "corresponding atmospheric grid.";
        return false;
      } else {
        // Longitude grid ok, add it to grids
        grids[1] = lon_retr;
      }
    }
  }
  return true;
}

void jacobian_type_extrapol(ArrayOfGridPos& gp) {
  for (Index i = 0; i < gp.nelem(); i++) {
    if (gp[i].fd[0] < 0) {
      gp[i].fd[0] = 0;
      gp[i].fd[1] = 1;
    } else if (gp[i].fd[0] > 1) {
      gp[i].fd[0] = 1;
      gp[i].fd[1] = 0;
    }
  }
}

void polynomial_basis_func(Vector& b,
                           const Vector& x,
                           const Index& poly_coeff) {
  const Index l = x.nelem();

  assert(l > poly_coeff);

  if (b.nelem() != l) b.resize(l);

  if (poly_coeff == 0) {
    b = 1.0;
  } else {
    const Numeric xmin = min(x);
    const Numeric dx = 0.5 * (max(x) - xmin);
    //
    for (Index i = 0; i < l; i++) {
      b[i] = (x[i] - xmin) / dx - 1.0;
      b[i] = pow(b[i], int(poly_coeff));
    }
    //
    b -= mean(b);
  }
}

void calcBaselineFit(Vector& y_baseline,
                     const Vector& x,
                     const Index& mblock_index,
                     const Sparse& sensor_response,
                     const ArrayOfIndex& sensor_response_pol_grid,
                     const Vector& sensor_response_f_grid,
                     const Matrix& sensor_response_dlos_grid,
                     const RetrievalQuantity& rq,
                     const Index rq_index,
                     const ArrayOfArrayOfIndex& jacobian_indices) {
  bool is_sine_fit = false;
  if (rq == Jacobian::Sensor::Polyfit) {
    is_sine_fit = false;
  } else if (rq == Jacobian::Sensor::Sinefit) {
    is_sine_fit = true;
  } else {
    throw runtime_error(
        "Retrieval quantity is neither a polynomial or a sine "
        " baseline fit.");
  }

  // Size and check of sensor_response
  //
  const Index nf = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nlos = sensor_response_dlos_grid.nrows();

  // Evaluate basis functions for fits.
  Vector w, s, c;
  if (is_sine_fit) {
    s.resize(nf);
    c.resize(nf);
    Numeric period = rq.Grids()[0][0];
    for (Index f = 0; f < nf; f++) {
      Numeric a = (sensor_response_f_grid[f] - sensor_response_f_grid[0]) * 2 *
                  PI / period;
      s[f] = sin(a);
      c[f] = cos(a);
    }
  } else {
    Numeric poly_coeff = rq.Grids()[0][0];
    polynomial_basis_func(
        w, sensor_response_f_grid, static_cast<Index>(poly_coeff));
  }

  // Compute baseline
  ArrayOfVector jg = rq.Grids();
  const Index n1 = jg[1].nelem();
  const Index n2 = jg[2].nelem();
  const Index n3 = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row4 = rowind.get_start();
  Index col4 = jacobian_indices[rq_index][0];

  if (n3 > 1) {
    col4 += mblock_index * n2 * n1;
  }

  for (Index l = 0; l < nlos; l++) {
    const Index row3 = row4 + l * nf * npol;
    const Index col3 = col4 + l * n1 * (is_sine_fit ? 2 : 1);

    for (Index f = 0; f < nf; f++) {
      const Index row2 = row3 + f * npol;

      for (Index p = 0; p < npol; p++) {
        Index col1 = col3;
        if (n1 > 1) {
          col1 += p;
        }
        if (is_sine_fit) {
          y_baseline[row2 + p] += x[col1] * s[f] + x[col1 + 1] * c[f];
        } else {
          y_baseline[row2 + p] += w[f] * x[col1];
        }
      }
    }
  }
}

void vmrunitscf(Numeric& x,
                const String& unit,
                const Numeric& vmr,
                const Numeric& p,
                const Numeric& t) {
  if (unit == "rel" || unit == "logrel") {
    x = 1;
  } else if (unit == "vmr") {
    if (vmr == 0) {
      x = 0;
      return;
    }
    x = 1 / vmr;
  } else if (unit == "nd") {
    if (vmr == 0) {
      x = 0;
      return;
    }
    x = 1 / (vmr * number_density(p, t));
  } else {
    ostringstream os;
    os << "Allowed options for gas species jacobians are "
          "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "
       << unit << std::endl;
    throw std::runtime_error(os.str());
  }
}

void dxdvmrscf(Numeric& x,
               const String& unit,
               const Numeric& vmr,
               const Numeric& p,
               const Numeric& t) {
  if (unit == "rel" || unit == "logrel") {
    x = vmr;
  } else if (unit == "vmr") {
    x = 1;
  } else if (unit == "nd") {
    x = 1 / number_density(p, t);
  } else {
    ostringstream os;
    os << "Allowed options for gas species jacobians are "
          "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "
       << unit << std::endl;
    throw std::runtime_error(os.str());
  }
}

void get_diydx(VectorView diy1,
               VectorView diy2,
               ConstMatrixView ImT,
               ConstMatrixView cumulative_transmission,
               ConstMatrixView dT1,
               ConstMatrixView dT2,
               ConstVectorView iYmJ,
               ConstVectorView dJ1,
               ConstVectorView dJ2,
               const Index stokes_dim,
               const bool transmission_only) {
  /*
   * Solves 
   * 
   * diy1 = PiT [ dT1 iYmJ + (1-T) dJ1 ],
   * 
   * and
   * 
   * diy2 += PiT [ dT2 iYmJ + (1-T) dJ2 ],
   * 
   * where diy2 is diy1 from a prior layer
   * 
   * FIXME:  Needs HSE
  */

  // Computation vectors
  Vector a(stokes_dim), b(stokes_dim);

  // The first time a level is involved in a layer
  mult(a, dT1, iYmJ);
  if (not transmission_only) {
    mult(b, ImT, dJ1);
    a += b;
  }
  mult(diy1, cumulative_transmission, a);

  // The second time a level is involved in a layer
  mult(a, dT2, iYmJ);
  if (not transmission_only) {
    mult(b, ImT, dJ2);
    a += b;
  }
  mult(b, cumulative_transmission, a);
  diy2 += b;
}

//======================================================================
//             Propmat partials descriptions
//======================================================================

ArrayOfIndex equivalent_propmattype_indexes(const ArrayOfRetrievalQuantity& js) {
  ArrayOfIndex pos;
  pos.reserve(js.nelem());
  for (Index i = 0; i < js.nelem(); i++)
    if (not (js[i] == Jacobian::Type::Special or js[i] == Jacobian::Type::Sensor))
      pos.push_back(i);
  return pos;
}

Index equivalent_propmattype_index(const ArrayOfRetrievalQuantity& js,
                                   const Index i) noexcept {
  return equivalent_propmattype_indexes(js)[i];
}

bool is_wind_parameter(const RetrievalQuantity& t) noexcept {
  return t.Target().isWind();
}

bool is_frequency_parameter(const RetrievalQuantity& t) noexcept {
  return t.Target().isWind() or t.Target().isFrequency();
}

bool is_derived_magnetic_parameter(const RetrievalQuantity& t) noexcept {
  return t == Jacobian::Atm::MagneticMagnitude;
}

bool is_magnetic_parameter(const RetrievalQuantity& t) noexcept {
  return t.Target().isMagnetic();
}

bool is_nlte_parameter(const RetrievalQuantity& t) noexcept {
  return t == Jacobian::Line::NLTE;
}

#define ISLINESHAPETYPE(X)                                               \
  bool is_pressure_broadening_##X(const RetrievalQuantity& t) noexcept { \
    return t == Jacobian::Line::Shape##X##X0 or                          \
           t == Jacobian::Line::Shape##X##X1 or                          \
           t == Jacobian::Line::Shape##X##X2 or                          \
           t == Jacobian::Line::Shape##X##X3;                            \
  }
ISLINESHAPETYPE(G0)
ISLINESHAPETYPE(D0)
ISLINESHAPETYPE(G2)
ISLINESHAPETYPE(D2)
ISLINESHAPETYPE(FVC)
ISLINESHAPETYPE(ETA)
ISLINESHAPETYPE(Y)
ISLINESHAPETYPE(G)
ISLINESHAPETYPE(DV)
#undef ISLINESHAPETYPE

#define VARISLINESHAPEPARAM(X, Y) (t == Jacobian::Line::Shape##X##Y)
bool is_lineshape_parameter_X0(const RetrievalQuantity& t) noexcept {
  return VARISLINESHAPEPARAM(G0, X0) or VARISLINESHAPEPARAM(D0, X0) or
         VARISLINESHAPEPARAM(G2, X0) or VARISLINESHAPEPARAM(D2, X0) or
         VARISLINESHAPEPARAM(FVC, X0) or VARISLINESHAPEPARAM(ETA, X0) or
         VARISLINESHAPEPARAM(Y, X0) or VARISLINESHAPEPARAM(G, X0) or
         VARISLINESHAPEPARAM(DV, X0);
}

bool is_lineshape_parameter_X1(const RetrievalQuantity& t) noexcept {
  return VARISLINESHAPEPARAM(G0, X1) or VARISLINESHAPEPARAM(D0, X1) or
         VARISLINESHAPEPARAM(G2, X1) or VARISLINESHAPEPARAM(D2, X1) or
         VARISLINESHAPEPARAM(FVC, X1) or VARISLINESHAPEPARAM(ETA, X1) or
         VARISLINESHAPEPARAM(Y, X1) or VARISLINESHAPEPARAM(G, X1) or
         VARISLINESHAPEPARAM(DV, X1);
}

bool is_lineshape_parameter_X2(const RetrievalQuantity& t) noexcept {
  return VARISLINESHAPEPARAM(G0, X2) or VARISLINESHAPEPARAM(D0, X2) or
         VARISLINESHAPEPARAM(G2, X2) or VARISLINESHAPEPARAM(D2, X2) or
         VARISLINESHAPEPARAM(FVC, X2) or VARISLINESHAPEPARAM(ETA, X2) or
         VARISLINESHAPEPARAM(Y, X2) or VARISLINESHAPEPARAM(G, X2) or
         VARISLINESHAPEPARAM(DV, X2);
}

bool is_lineshape_parameter_X3(const RetrievalQuantity& t) noexcept {
  return VARISLINESHAPEPARAM(G0, X3) or VARISLINESHAPEPARAM(D0, X3) or
         VARISLINESHAPEPARAM(G2, X3) or VARISLINESHAPEPARAM(D2, X3) or
         VARISLINESHAPEPARAM(FVC, X3) or VARISLINESHAPEPARAM(ETA, X3) or
         VARISLINESHAPEPARAM(Y, X3) or VARISLINESHAPEPARAM(G, X3) or
         VARISLINESHAPEPARAM(DV, X3);
}
#undef VARISLINESHAPEPARAM

bool is_lineshape_parameter_bar_linemixing(
    const RetrievalQuantity& t) noexcept {
  return is_pressure_broadening_G0(t) or is_pressure_broadening_D0(t) or
         is_pressure_broadening_G2(t) or is_pressure_broadening_D2(t) or
         is_pressure_broadening_FVC(t) or is_pressure_broadening_ETA(t);
}

bool is_lineshape_parameter(const RetrievalQuantity& t) noexcept {
  return is_pressure_broadening_G0(t) or is_pressure_broadening_D0(t) or
         is_pressure_broadening_G2(t) or is_pressure_broadening_D2(t) or
         is_pressure_broadening_FVC(t) or is_pressure_broadening_ETA(t) or
         is_pressure_broadening_Y(t) or is_pressure_broadening_G(t) or
         is_pressure_broadening_DV(t);
}

bool is_line_parameter(const RetrievalQuantity& t) noexcept {
  return t == Jacobian::Type::Line;
}

bool supports_CIA(const ArrayOfRetrievalQuantity& js) {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return (j == Jacobian::Atm::Temperature or j == Jacobian::Line::VMR or is_frequency_parameter(j));});
}

bool supports_hitran_xsec(const ArrayOfRetrievalQuantity& js) {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return (j == Jacobian::Atm::Temperature or j == Jacobian::Line::VMR or is_frequency_parameter(j));});
}

bool supports_continuum(const ArrayOfRetrievalQuantity& js) {
  if (std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_line_parameter(j);}))
    throw std::runtime_error("Line specific parameters are not supported while using continuum tags.\nWe do not track what lines are in the continuum.\n");
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return (j == Jacobian::Atm::Temperature or is_frequency_parameter(j));});
}

bool supports_relaxation_matrix(const ArrayOfRetrievalQuantity& js) {
  if (std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_line_parameter(j);}))
    throw std::runtime_error("Line specific parameters are not supported while\n using the relaxation matrix line mixing routine.\n We do not yet track individual lines in the relaxation matrix calculations.\n");
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return (j == Jacobian::Atm::Temperature or is_frequency_parameter(j));});
}

bool supports_lookup(const ArrayOfRetrievalQuantity& js) {
  if (std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_line_parameter(j);}))
    throw std::runtime_error("Line specific parameters are not supported while using Lookup table.\nWe do not track lines in the Lookup.\n");
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return (j == Jacobian::Atm::Temperature or j == Jacobian::Special::ArrayOfSpeciesTagVMR or is_frequency_parameter(j));});
}

bool supports_faraday(const ArrayOfRetrievalQuantity& js) {
  if (std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_derived_magnetic_parameter(j);}))
    throw std::runtime_error("This method does not yet support Zeeman-style magnetic Jacobian calculations.\n Please use u, v, and w Jacobians instead.\n");
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return j == Jacobian::Atm::MagneticU or j == Jacobian::Atm::MagneticV or j == Jacobian::Atm::MagneticW or j == Jacobian::Atm::Electrons or is_frequency_parameter(j);});
}

bool supports_particles(const ArrayOfRetrievalQuantity& js) {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return j == Jacobian::Atm::Temperature or j == Jacobian::Atm::Particulates;});
}

bool supports_propmat_clearsky(const ArrayOfRetrievalQuantity& js) {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return not (j == Jacobian::Type::Sensor);});
}

bool species_match(const RetrievalQuantity& rq, const ArrayOfSpeciesTag& ast) {
  if (rq.QuantumIdentity().Type() == QuantumIdentifier::ALL) {  // Single tag
    for (const auto& s : ast) {
      if (rq.QuantumIdentity().Species() not_eq s.Species())
        return false;  // Species must match
      if (rq.QuantumIdentity().Isotopologue() >= 0 and s.Isotopologue() >= 0 and
          rq.QuantumIdentity().Isotopologue() not_eq s.Isotopologue())
        return false;  // Isotopologue must match or be undefined
    }
  } else {
    ArrayOfSpeciesTag test;
    array_species_tag_from_string(test, rq.Subtag());
    if (ast not_eq test)
      return false;  // Match single tag perfectly or throw out
  }
  return true;
}

bool species_match(const RetrievalQuantity& rq, const Index species) {
  if (rq == Jacobian::Line::VMR and rq.QuantumIdentity().Species() == species)
    return true;
  return false;
}

bool species_iso_match(const RetrievalQuantity& rq,
                       const Index species,
                       const Index iso) {
  const SpeciesRecord& spr = global_data::species_data[species];
  if (species_match(rq, species))
    if (rq.QuantumIdentity().Isotopologue() == iso or
        spr.Isotopologue().nelem() == iso or iso < 0)
      return true;
  return false;
}

bool do_temperature_jacobian(const ArrayOfRetrievalQuantity& js) noexcept {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return j == Jacobian::Atm::Temperature;});
}

jacobianVMRcheck do_vmr_jacobian(const ArrayOfRetrievalQuantity& js,
                                 const QuantumIdentifier& line_qid) noexcept {
  auto p = std::find_if(js.cbegin(), js.cend(), [&line_qid](auto& j){return j == Jacobian::Line::VMR and j.QuantumIdentity().In(line_qid);});
  if (p not_eq js.cend())
    return {true, p -> QuantumIdentity()};
  else
    return {false, line_qid};
}

bool do_line_center_jacobian(const ArrayOfRetrievalQuantity& js) noexcept {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return j == Jacobian::Line::Center;});
}

bool do_frequency_jacobian(const ArrayOfRetrievalQuantity& js) noexcept {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_frequency_parameter(j);});
}

bool do_magnetic_jacobian(const ArrayOfRetrievalQuantity& js) noexcept {
  return std::any_of(js.cbegin(), js.cend(), [](auto& j){return is_magnetic_parameter(j);});
}

Numeric temperature_perturbation(const ArrayOfRetrievalQuantity& js) noexcept {
  auto p = std::find_if(js.cbegin(), js.cend(), [](auto& j){return j == Jacobian::Atm::Temperature;});
  if (p not_eq js.cend())
    return p -> Target().Perturbation();
  else
    return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric frequency_perturbation(const ArrayOfRetrievalQuantity& js) noexcept {
  auto p = std::find_if(js.cbegin(), js.cend(), [](auto& j){return is_frequency_parameter(j);});
  if (p not_eq js.cend())
    return p -> Target().Perturbation();
  else
    return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric magnetic_field_perturbation(const ArrayOfRetrievalQuantity& js) noexcept {
  auto p = std::find_if(js.cbegin(), js.cend(), [](auto& j){return is_magnetic_parameter(j);});
  if (p not_eq js.cend())
    return p -> Target().Perturbation();
  else
    return std::numeric_limits<Numeric>::quiet_NaN();
}
