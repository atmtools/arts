/* Copyright (C) 2004-2012 Mattias Ekstrom <ekstrom@rss.chalmers.se>

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
   USA. 
*/

/**
  @file   m_jacobian.cc
  @author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  @date   2004-09-14

  @brief  Workspace functions related to the jacobian.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <string>
#include "absorption.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "interpolation_poly.h"
#include "jacobian.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"
#include "rte.h"

extern const Numeric PI;

extern const String ABSSPECIES_MAINTAG;
extern const String FREQUENCY_MAINTAG;
extern const String FREQUENCY_SUBTAG_0;
extern const String FREQUENCY_SUBTAG_1;
extern const String POINTING_MAINTAG;
extern const String POINTING_SUBTAG_A;
extern const String POINTING_CALCMODE_A;
extern const String POINTING_CALCMODE_B;
extern const String POLYFIT_MAINTAG;
extern const String SCATSPECIES_MAINTAG;
extern const String SINEFIT_MAINTAG;
extern const String TEMPERATURE_MAINTAG;
extern const String NLTE_MAINTAG;
extern const String WIND_MAINTAG;
extern const String MAGFIELD_MAINTAG;
extern const String FLUX_MAINTAG;
extern const String PROPMAT_SUBSUBTAG;
extern const String ELECTRONS_MAINTAG;
extern const String PARTICULATES_MAINTAG;
extern const String CATALOGPARAMETER_MAINTAG;

extern const String SURFACE_MAINTAG;

// Generic modes
extern const String PRESSUREBROADENINGGAMMA_MODE;
extern const String LINESTRENGTH_MODE;
extern const String LINECENTER_MODE;
extern const String LINEMIXINGY_MODE;
extern const String LINEMIXINGG_MODE;
extern const String LINEMIXINGDF_MODE;

// Modes for "some" catalogs
//  Pressure Broadening
extern const String SELFBROADENING_MODE;
extern const String FOREIGNBROADENING_MODE;
extern const String WATERBROADENING_MODE;
extern const String SELFBROADENINGEXPONENT_MODE;
extern const String FOREIGNBROADENINGEXPONENT_MODE;
extern const String WATERBROADENINGEXPONENT_MODE;
extern const String SELFPRESSURESHIFT_MODE;
extern const String FOREIGNPRESSURESHIFT_MODE;
extern const String WATERPRESSURESHIFT_MODE;

//  Line Mixing
extern const String LINEMIXINGY0_MODE;
extern const String LINEMIXINGG0_MODE;
extern const String LINEMIXINGDF0_MODE;
extern const String LINEMIXINGY1_MODE;
extern const String LINEMIXINGG1_MODE;
extern const String LINEMIXINGDF1_MODE;
extern const String LINEMIXINGYEXPONENT_MODE;
extern const String LINEMIXINGGEXPONENT_MODE;
extern const String LINEMIXINGDFEXPONENT_MODE;

/*===========================================================================
  === The methods, with general methods first followed by the Add/Calc method
  === pairs for each retrieval quantity.
  ===========================================================================*/

//----------------------------------------------------------------------------
// Basic methods:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcDoNothing(Matrix& jacobian _U_,
                           const Index& mblock_index _U_,
                           const Vector& iyb _U_,
                           const Vector& yb _U_,
                           const Verbosity&) {
  /* Nothing to do here for the analytical case, this function just exists
   to satisfy the required inputs and outputs of the jacobian_agenda */
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianClose(Workspace& ws,
                   Index& jacobian_do,
                   Agenda& jacobian_agenda,
                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                   const Verbosity& verbosity) {
  // Make sure that the array is not empty
  if (jacobian_quantities.empty())
    throw runtime_error(
        "No retrieval quantities has been added to *jacobian_quantities*.");

  jacobian_agenda.check(ws, verbosity);
  jacobian_do = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianInit(ArrayOfRetrievalQuantity& jacobian_quantities,
                  Agenda& jacobian_agenda,
                  const Verbosity&) {
  jacobian_quantities.resize(0);
  jacobian_agenda = Agenda();
  jacobian_agenda.set_name("jacobian_agenda");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianOff(Index& jacobian_do,
                 Agenda& jacobian_agenda,
                 ArrayOfRetrievalQuantity& jacobian_quantities,
                 const Verbosity& verbosity) {
  jacobian_do = 0;
  jacobianInit(jacobian_quantities, jacobian_agenda, verbosity);
}

//----------------------------------------------------------------------------
// Absorption species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddAbsSpecies(Workspace&,
                           ArrayOfRetrievalQuantity& jq,
                           Agenda& jacobian_agenda,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Vector& rq_p_grid,
                           const Vector& rq_lat_grid,
                           const Vector& rq_lon_grid,
                           const String& species,
                           const String& mode,
                           const Index& for_species_tag,
                           const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  QuantumIdentifier qi;
  if (not for_species_tag) {
    ArrayOfSpeciesTag test;
    array_species_tag_from_string(test, species);
    if (test.nelem() not_eq 1)
      throw std::runtime_error(
          "Trying to add a species as a species tag of multiple species.\n"
          "This is not supported.  Please give just a single species instead.\n"
          "Otherwise consider if you intended for_species_tag to be evaluated true.\n");
    qi.SetAll();
    qi.SetIsotopologue(test[0].Isotopologue());
    qi.SetSpecies(test[0].Species());
  }

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == ABSSPECIES_MAINTAG &&
        jq[it].SubSubtag() != PROPMAT_SUBSUBTAG && jq[it].Subtag() == species) {
      ostringstream os;
      os << "The gas species:\n"
         << species << "\nis already included in "
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    } else if (jq[it].MainTag() == ABSSPECIES_MAINTAG &&
               jq[it].SubSubtag() == PROPMAT_SUBSUBTAG) {
      if (SpeciesTag(jq[it].Subtag()) == SpeciesTag(species)) {
        ostringstream os;
        os << "The atmospheric species of:\n"
           << species << "\nis already included in "
           << "*jacobian_quantities*.";
        throw runtime_error(os.str());
      }
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Check that mode is correct
  if (mode != "vmr" && mode != "nd" && mode != "rel" && mode != "rh" &&
      mode != "q") {
    throw runtime_error(
        "The retrieval mode can only be \"vmr\", \"nd\", "
        "\"rel\", \"rh\" or \"q\".");
  }
  if ((mode == "rh" || mode == "q") && species.substr(0, 3) != "H2O") {
    throw runtime_error(
        "Retrieval modes \"rh\" and \"q\" can only be applied "
        "on species starting with H2O.");
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(ABSSPECIES_MAINTAG);
  rq.Subtag(species);
  rq.Mode(mode);
  rq.Analytical(1);
  rq.Perturbation(0.001);
  rq.Grids(grids);
  if (not for_species_tag) {
    rq.SubSubtag(PROPMAT_SUBSUBTAG);
    rq.PropType(JacPropMatType::VMR);
  } else
    rq.PropType(JacPropMatType::NotPropagationMatrixType);

  rq.QuantumIdentity(qi);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Frequency shift
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqShift(Workspace& ws _U_,
                          ArrayOfRetrievalQuantity& jacobian_quantities,
                          Agenda& jacobian_agenda,
                          const Vector& f_grid,
                          const Numeric& df,
                          const Verbosity&) {
  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it].MainTag() == FREQUENCY_MAINTAG &&
        jacobian_quantities[it].Subtag() == FREQUENCY_SUBTAG_0) {
      ostringstream os;
      os << "Fit of frequency shift is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Checks of frequencies
  if (df <= 0) throw runtime_error("The argument *df* must be > 0.");
  if (df > 1e6)
    throw runtime_error("The argument *df* is not allowed to exceed 1 MHz.");
  const Index nf = f_grid.nelem();
  if (nf < 2)
    throw runtime_error(
        "Frequency shifts and *f_grid* of length 1 can "
        "not be combined.");
  const Numeric maxdf = f_grid[nf - 1] - f_grid[nf - 2];
  if (df > maxdf) {
    ostringstream os;
    os << "The value of *df* is too big with respect to spacing of "
       << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
       << "between the two last elements of *f_grid*.\n"
       << "This spacing is   : " << maxdf / 1e3 << " kHz\n"
       << "The value of df is: " << df / 1e3 << " kHz";
    throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(FREQUENCY_MAINTAG);
  rq.Subtag(FREQUENCY_SUBTAG_0);
  rq.Mode("");
  rq.Analytical(0);
  rq.Perturbation(df);

  // Dummy vector of length 1
  Vector grid(1, 0);
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.append("jacobianCalcFreqShift", "");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqShift(Matrix& jacobian,
                           const Index& mblock_index,
                           const Vector& iyb,
                           const Vector& yb,
                           const Index& stokes_dim,
                           const Vector& f_grid,
                           const Matrix& mblock_dlos_grid,
                           const Sparse& sensor_response,
                           const ArrayOfRetrievalQuantity& jacobian_quantities,
                           const Verbosity&) {
  // Set some useful (and needed) variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n].MainTag() == FREQUENCY_MAINTAG &&
        jacobian_quantities[n].Subtag() == FREQUENCY_SUBTAG_0) {
      bool any_affine;
      ArrayOfArrayOfIndex jacobian_indices;
      jac_ranges_indices(
          jacobian_indices, any_affine, jacobian_quantities, true);
      //
      found = true;
      rq = jacobian_quantities[n];
      ji = jacobian_indices[n];
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no such frequency retrieval quantity defined.\n");
  }

  // Check that sensor_response is consistent with yb and iyb
  //
  if (sensor_response.nrows() != yb.nelem())
    throw runtime_error("Mismatch in size between *sensor_response* and *yb*.");
  if (sensor_response.ncols() != iyb.nelem())
    throw runtime_error(
        "Mismatch in size between *sensor_response* and *iyb*.");

  // Get disturbed (part of) y
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    const Index nf2 = f_grid.nelem();
    const Index nlos2 = mblock_dlos_grid.nrows();
    const Index niyb = nf2 * nlos2 * stokes_dim;

    // Interpolation weights
    //
    const Index porder = 3;
    //
    ArrayOfGridPosPoly gp(nf2);
    Matrix itw(nf2, porder + 1);
    Vector fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Perturbation();
    gridpos_poly(gp, f_grid, fg_new, porder, 1.0);
    interpweights(itw, gp);

    // Do interpolation
    for (Index ilos = 0; ilos < nlos2; ilos++) {
      const Index row0 = ilos * nf2 * stokes_dim;

      for (Index is = 0; is < stokes_dim; is++) {
        interp(iyb2[Range(row0 + is, nf2, stokes_dim)],
               itw,
               iyb[Range(row0 + is, nf2, stokes_dim)],
               gp);
      }
    }

    // Determine difference
    //
    mult(dy, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (dy[i] - yb[i]) / rq.Perturbation();
    }
  }

  //--- Set jacobian ---
  assert(rq.Grids()[0].nelem() == 1);
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  jacobian(rowind, ji[0]) = dy;
}

//----------------------------------------------------------------------------
// Frequency stretch
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqStretch(Workspace& ws _U_,
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            Agenda& jacobian_agenda,
                            const Vector& f_grid,
                            const Numeric& df,
                            const Verbosity&) {
  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it].MainTag() == FREQUENCY_MAINTAG &&
        jacobian_quantities[it].Subtag() == FREQUENCY_SUBTAG_1) {
      ostringstream os;
      os << "Fit of frequency stretch is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Checks of df
  if (df <= 0) throw runtime_error("The argument *df* must be > 0.");
  if (df > 1e6)
    throw runtime_error("The argument *df* is not allowed to exceed 1 MHz.");
  const Index nf = f_grid.nelem();
  const Numeric maxdf = f_grid[nf - 1] - f_grid[nf - 2];
  if (df > maxdf) {
    ostringstream os;
    os << "The value of *df* is too big with respect to spacing of "
       << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
       << "between the two last elements of *f_grid*.\n"
       << "This spacing is   : " << maxdf / 1e3 << " kHz\n"
       << "The value of df is: " << df / 1e3 << " kHz";
    throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(FREQUENCY_MAINTAG);
  rq.Subtag(FREQUENCY_SUBTAG_1);
  rq.Mode("");
  rq.Analytical(0);
  rq.Perturbation(df);

  // Dummy vector of length 1
  Vector grid(1, 0);
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.append("jacobianCalcFreqStretch", "");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqStretch(
    Matrix& jacobian,
    const Index& mblock_index,
    const Vector& iyb,
    const Vector& yb,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Matrix& mblock_dlos_grid,
    const Sparse& sensor_response,
    const ArrayOfIndex& sensor_response_pol_grid,
    const Vector& sensor_response_f_grid,
    const Matrix& sensor_response_dlos_grid,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Verbosity&) {
  // The code here is close to identical to the one for Shift. The main
  // difference is that dy is weighted with poly_order 1 basis function.

  // Set some useful (and needed) variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n].MainTag() == FREQUENCY_MAINTAG &&
        jacobian_quantities[n].Subtag() == FREQUENCY_SUBTAG_1) {
      bool any_affine;
      ArrayOfArrayOfIndex jacobian_indices;
      jac_ranges_indices(
          jacobian_indices, any_affine, jacobian_quantities, true);
      //
      found = true;
      rq = jacobian_quantities[n];
      ji = jacobian_indices[n];
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no such frequency retrieval quantity defined.\n");
  }

  // Check that sensor_response is consistent with yb and iyb
  //
  if (sensor_response.nrows() != yb.nelem())
    throw runtime_error("Mismatch in size between *sensor_response* and *yb*.");
  if (sensor_response.ncols() != iyb.nelem())
    throw runtime_error(
        "Mismatch in size between *sensor_response* and *iyb*.");

  // Get disturbed (part of) y
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    const Index nf2 = f_grid.nelem();
    const Index nlos2 = mblock_dlos_grid.nrows();
    const Index niyb = nf2 * nlos2 * stokes_dim;

    // Interpolation weights
    //
    const Index porder = 3;
    //
    ArrayOfGridPosPoly gp(nf2);
    Matrix itw(nf2, porder + 1);
    Vector fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Perturbation();
    gridpos_poly(gp, f_grid, fg_new, porder, 1.0);
    interpweights(itw, gp);

    // Do interpolation
    for (Index ilos = 0; ilos < nlos2; ilos++) {
      const Index row0 = ilos * nf2 * stokes_dim;

      for (Index is = 0; is < stokes_dim; is++) {
        interp(iyb2[Range(row0 + is, nf2, stokes_dim)],
               itw,
               iyb[Range(row0 + is, nf2, stokes_dim)],
               gp);
      }
    }

    // Determine difference
    //
    mult(dy, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (dy[i] - yb[i]) / rq.Perturbation();
    }

    // dy above corresponds now to shift. Convert to stretch:
    //
    Vector w;
    polynomial_basis_func(w, sensor_response_f_grid, 1);
    //
    const Index nf = sensor_response_f_grid.nelem();
    const Index npol = sensor_response_pol_grid.nelem();
    const Index nlos = sensor_response_dlos_grid.nrows();
    //
    for (Index l = 0; l < nlos; l++) {
      for (Index f = 0; f < nf; f++) {
        const Index row1 = (l * nf + f) * npol;
        for (Index p = 0; p < npol; p++) {
          dy[row1 + p] *= w[f];
        }
      }
    }
  }

  //--- Set jacobians ---
  assert(rq.Grids()[0].nelem() == 1);
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  jacobian(rowind, ji[0]) = dy;
}

//----------------------------------------------------------------------------
// Pointing:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPointingZa(Workspace& ws _U_,
                           ArrayOfRetrievalQuantity& jacobian_quantities,
                           Agenda& jacobian_agenda,
                           const Matrix& sensor_pos,
                           const Vector& sensor_time,
                           const Index& poly_order,
                           const String& calcmode,
                           const Numeric& dza,
                           const Verbosity&) {
  // Check that poly_order is -1 or positive
  if (poly_order < -1)
    throw runtime_error(
        "The polynomial order has to be positive or -1 for gitter.");

  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it].MainTag() == POINTING_MAINTAG &&
        jacobian_quantities[it].Subtag() == POINTING_SUBTAG_A) {
      ostringstream os;
      os << "Fit of zenith angle pointing off-set is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Checks of dza
  if (dza <= 0) throw runtime_error("The argument *dza* must be > 0.");
  if (dza > 0.1)
    throw runtime_error("The argument *dza* is not allowed to exceed 0.1 deg.");

  // Check that sensor_time is consistent with sensor_pos
  if (sensor_time.nelem() != sensor_pos.nrows()) {
    ostringstream os;
    os << "The WSV *sensor_time* must be defined for every "
       << "measurement block.\n";
    throw runtime_error(os.str());
  }

  // Do not allow that *poly_order* is not too large compared to *sensor_time*
  if (poly_order > sensor_time.nelem() - 1) {
    throw runtime_error(
        "The polynomial order can not be >= length of *sensor_time*.");
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(POINTING_MAINTAG);
  rq.Subtag(POINTING_SUBTAG_A);
  rq.Analytical(0);
  rq.Perturbation(dza);

  // To store the value or the polynomial order, create a vector with length
  // poly_order+1, in case of gitter set the size of the grid vector to be the
  // number of measurement blocks, all elements set to -1.
  Vector grid(0, poly_order + 1, 1);
  if (poly_order == -1) {
    grid.resize(sensor_pos.nrows());
    grid = -1.0;
  }
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  if (calcmode == "recalc") {
    rq.Mode(POINTING_CALCMODE_A);
    jacobian_agenda.append("jacobianCalcPointingZaRecalc", "");
  } else if (calcmode == "interp") {
    rq.Mode(POINTING_CALCMODE_B);
    jacobian_agenda.append("jacobianCalcPointingZaInterp", "");
  } else
    throw runtime_error(
        "Possible choices for *calcmode* are \"recalc\" and \"interp\".");

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaInterp(
    Matrix& jacobian,
    const Index& mblock_index,
    const Vector& iyb,
    const Vector& yb _U_,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Matrix& DEBUG_ONLY(sensor_los),
    const Matrix& mblock_dlos_grid,
    const Sparse& sensor_response,
    const Vector& sensor_time,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Verbosity&) {
  if (mblock_dlos_grid.nrows() < 2)
    throw runtime_error(
        "The method demands that *mblock_dlos_grid* has "
        "more than one row.");

  if (!(is_increasing(mblock_dlos_grid(joker, 0)) ||
        is_decreasing(mblock_dlos_grid(joker, 0))))
    throw runtime_error(
        "The method demands that the zenith angles in "
        "*mblock_dlos_grid* are sorted (increasing or decreasing).");

  // Set some useful variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n].MainTag() == POINTING_MAINTAG &&
        jacobian_quantities[n].Subtag() == POINTING_SUBTAG_A &&
        jacobian_quantities[n].Mode() == POINTING_CALCMODE_B) {
      bool any_affine;
      ArrayOfArrayOfIndex jacobian_indices;
      jac_ranges_indices(
          jacobian_indices, any_affine, jacobian_quantities, true);
      //
      found = true;
      rq = jacobian_quantities[n];
      ji = jacobian_indices[n];
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no such pointing retrieval quantity defined.\n");
  }

  // Get "dy", by inter/extra-polation of existing iyb
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    // Sizes
    const Index nf = f_grid.nelem();
    const Index nza = mblock_dlos_grid.nrows();

    // Shifted zenith angles
    Vector za1 = mblock_dlos_grid(joker, 0);
    za1 -= rq.Perturbation();
    Vector za2 = mblock_dlos_grid(joker, 0);
    za2 += rq.Perturbation();

    // Find interpolation weights
    ArrayOfGridPos gp1(nza), gp2(nza);
    gridpos(
        gp1, mblock_dlos_grid(joker, 0), za1, 1e6);  // Note huge extrapolation!
    gridpos(
        gp2, mblock_dlos_grid(joker, 0), za2, 1e6);  // Note huge extrapolation!
    Matrix itw1(nza, 2), itw2(nza, 2);
    interpweights(itw1, gp1);
    interpweights(itw2, gp2);

    // Make interpolation (for all azimuth angles, frequencies and Stokes)
    //
    Vector iyb1(iyb.nelem()), iyb2(iyb.nelem());
    //
    for (Index iza = 0; iza < nza; iza++) {
      for (Index iv = 0; iv < nf; iv++) {
        for (Index is = 0; is < stokes_dim; is++) {
          const Range r(iv * stokes_dim + is, nza, nf * stokes_dim);
          interp(iyb1[r], itw1, iyb[r], gp1);
          interp(iyb2[r], itw2, iyb[r], gp2);
        }
      }
    }

    // Apply sensor and take difference
    //
    Vector y1(n1y), y2(n1y);
    mult(y1, sensor_response, iyb1);
    mult(y2, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (y2[i] - y1[i]) / (2 * rq.Perturbation());
    }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row0 = rowind.get_start();

  // Handle pointing "jitter" seperately
  if (rq.Grids()[0][0] == -1)          // Not all values are set here,
  {                                    // but should already have been
    assert(lg == sensor_los.nrows());  // set to 0
    assert(rq.Grids()[0][mblock_index] == -1);
    jacobian(rowind, it + mblock_index) = dy;
  }

  // Polynomial representation
  else {
    Vector w;
    for (Index c = 0; c < lg; c++) {
      assert(Numeric(c) == rq.Grids()[0][c]);
      //
      polynomial_basis_func(w, sensor_time, c);
      //
      for (Index i = 0; i < n1y; i++) {
        jacobian(row0 + i, it + c) = w[mblock_index] * dy[i];
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaRecalc(
    Workspace& ws,
    Matrix& jacobian,
    const Index& mblock_index,
    const Vector& iyb _U_,
    const Vector& yb,
    const Index& atmosphere_dim,
    const EnergyLevelMap& nlte_field,              
    const Index& cloudbox_on,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Matrix& sensor_pos,
    const Matrix& sensor_los,
    const Matrix& transmitter_pos,
    const Matrix& mblock_dlos_grid,
    const Sparse& sensor_response,
    const Vector& sensor_time,
    const String& iy_unit,
    const Agenda& iy_main_agenda,
    const Agenda& geo_pos_agenda,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Verbosity& verbosity) {
  // Set some useful variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n].MainTag() == POINTING_MAINTAG &&
        jacobian_quantities[n].Subtag() == POINTING_SUBTAG_A &&
        jacobian_quantities[n].Mode() == POINTING_CALCMODE_A) {
      bool any_affine;
      ArrayOfArrayOfIndex jacobian_indices;
      jac_ranges_indices(
          jacobian_indices, any_affine, jacobian_quantities, true);
      //
      found = true;
      rq = jacobian_quantities[n];
      ji = jacobian_indices[n];
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no such pointing retrieval quantity defined.\n");
  }

  // Get "dy", by calling iyb_calc with shifted sensor_los.
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    Vector iyb2;
    Matrix los = sensor_los;
    Matrix geo_pos;
    ArrayOfVector iyb_aux;
    ArrayOfMatrix diyb_dx;

    los(joker, 0) += rq.Perturbation();

    iyb_calc(ws,
             iyb2,
             iyb_aux,
             diyb_dx,
             geo_pos,
             mblock_index,
             atmosphere_dim,
             nlte_field,
             cloudbox_on,
             stokes_dim,
             f_grid,
             sensor_pos,
             los,
             transmitter_pos,
             mblock_dlos_grid,
             iy_unit,
             iy_main_agenda,
             geo_pos_agenda,
             0,
             ArrayOfRetrievalQuantity(),
             ArrayOfArrayOfIndex(),
             ArrayOfString(),
             verbosity);

    // Apply sensor and take difference
    //
    mult(dy, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (dy[i] - yb[i]) / rq.Perturbation();
    }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row0 = rowind.get_start();

  // Handle "jitter" seperately
  if (rq.Grids()[0][0] == -1)          // Not all values are set here,
  {                                    // but should already have been
    assert(lg == sensor_los.nrows());  // set to 0
    assert(rq.Grids()[0][mblock_index] == -1);
    jacobian(rowind, it + mblock_index) = dy;
  }

  // Polynomial representation
  else {
    Vector w;
    for (Index c = 0; c < lg; c++) {
      assert(Numeric(c) == rq.Grids()[0][c]);
      //
      polynomial_basis_func(w, sensor_time, c);
      //
      for (Index i = 0; i < n1y; i++) {
        jacobian(row0 + i, it + c) = w[mblock_index] * dy[i];
      }
    }
  }
}

//----------------------------------------------------------------------------
// Polynomial baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPolyfit(Workspace& ws _U_,
                        ArrayOfRetrievalQuantity& jq,
                        Agenda& jacobian_agenda,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        const Matrix& sensor_response_dlos_grid,
                        const Matrix& sensor_pos,
                        const Index& poly_order,
                        const Index& no_pol_variation,
                        const Index& no_los_variation,
                        const Index& no_mblock_variation,
                        const Verbosity&) {
  // Check that poly_order is >= 0
  if (poly_order < 0)
    throw runtime_error("The polynomial order has to be >= 0.");

  // Check that polyfit is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == POLYFIT_MAINTAG) {
      ostringstream os;
      os << "Polynomial baseline fit is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // "Grids"
  //
  // Grid dimensions correspond here to
  //   1: polynomial order
  //   2: polarisation
  //   3: viewing direction
  //   4: measurement block
  //
  ArrayOfVector grids(4);
  //
  if (no_pol_variation)
    grids[1] = Vector(1, 1);
  else
    grids[1] = Vector(0, sensor_response_pol_grid.nelem(), 1);
  if (no_los_variation)
    grids[2] = Vector(1, 1);
  else
    grids[2] = Vector(0, sensor_response_dlos_grid.nrows(), 1);
  if (no_mblock_variation)
    grids[3] = Vector(1, 1);
  else
    grids[3] = Vector(0, sensor_pos.nrows(), 1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(POLYFIT_MAINTAG);
  rq.Mode("");
  rq.Analytical(0);
  rq.Perturbation(0);

  // Each polynomial coeff. is treated as a retrieval quantity
  //
  for (Index i = 0; i <= poly_order; i++) {
    ostringstream sstr;
    sstr << "Coefficient " << i;
    rq.Subtag(sstr.str());

    // Grid is a scalar, use polynomial coeff.
    grids[0] = Vector(1, (Numeric)i);
    rq.Grids(grids);

    // Add it to the *jacobian_quantities*
    jq.push_back(rq);

    // Add pointing method to the jacobian agenda
    jacobian_agenda.append("jacobianCalcPolyfit", i);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPolyfit(Matrix& jacobian,
                         const Index& mblock_index,
                         const Vector& iyb _U_,
                         const Vector& yb _U_,
                         const Sparse& sensor_response,
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Vector& sensor_response_f_grid,
                         const Matrix& sensor_response_dlos_grid,
                         const ArrayOfRetrievalQuantity& jacobian_quantities,
                         const Index& poly_coeff,
                         const Verbosity&) {
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  ostringstream sstr;
  sstr << "Coefficient " << poly_coeff;
  for (iq = 0; iq < jacobian_quantities.nelem() && !found; iq++) {
    if (jacobian_quantities[iq].MainTag() == POLYFIT_MAINTAG &&
        jacobian_quantities[iq].Subtag() == sstr.str()) {
      found = true;
      break;
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no Polyfit jacobian defined, in general "
        "or for the selected polynomial coefficient.\n");
  }

  // Size and check of sensor_response
  //
  const Index nf = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nlos = sensor_response_dlos_grid.nrows();

  // Make a vector with values to distribute over *jacobian*
  //
  Vector w;
  //
  polynomial_basis_func(w, sensor_response_f_grid, poly_coeff);

  // Fill J
  //
  ArrayOfArrayOfIndex jacobian_indices;
  {
    bool any_affine;
    jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities, true);
  }
  //
  ArrayOfVector jg = jacobian_quantities[iq].Grids();
  const Index n1 = jg[1].nelem();
  const Index n2 = jg[2].nelem();
  const Index n3 = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row4 = rowind.get_start();
  Index col4 = jacobian_indices[iq][0];

  if (n3 > 1) {
    col4 += mblock_index * n2 * n1;
  }

  for (Index l = 0; l < nlos; l++) {
    const Index row3 = row4 + l * nf * npol;
    const Index col3 = col4 + l * n1;

    for (Index f = 0; f < nf; f++) {
      const Index row2 = row3 + f * npol;

      for (Index p = 0; p < npol; p++) {
        Index col1 = col3;
        if (n1 > 1) {
          col1 += p;
        }

        jacobian(row2 + p, col1) = w[f];
      }
    }
  }
}

//----------------------------------------------------------------------------
// Scattering species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddScatSpecies(Workspace&,
                            ArrayOfRetrievalQuantity& jq,
                            Agenda& jacobian_agenda,
                            const Index& atmosphere_dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const Vector& rq_p_grid,
                            const Vector& rq_lat_grid,
                            const Vector& rq_lon_grid,
                            const String& species,
                            const String& quantity,
                            const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Check that this species+quantity combination is not already included in
  // the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == SCATSPECIES_MAINTAG && jq[it].Subtag() == species &&
        jq[it].SubSubtag() == quantity) {
      ostringstream os;
      os << "The combintaion of\n   scattering species: " << species
         << "\n   retrieval quantity: " << quantity
         << "\nis already included in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(SCATSPECIES_MAINTAG);
  rq.Subtag(species);
  rq.SubSubtag(quantity);
  rq.Analytical(1);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Sinusoidal baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSinefit(Workspace& ws _U_,
                        ArrayOfRetrievalQuantity& jq,
                        Agenda& jacobian_agenda,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        const Matrix& sensor_response_dlos_grid,
                        const Matrix& sensor_pos,
                        const Vector& period_lengths,
                        const Index& no_pol_variation,
                        const Index& no_los_variation,
                        const Index& no_mblock_variation,
                        const Verbosity&) {
  const Index np = period_lengths.nelem();

  // Check that poly_order is >= 0
  if (np == 0) throw runtime_error("No sinusoidal periods has benn given.");

  // Check that polyfit is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == SINEFIT_MAINTAG) {
      ostringstream os;
      os << "Polynomial baseline fit is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // "Grids"
  //
  // Grid dimensions correspond here to
  //   1: polynomial order
  //   2: polarisation
  //   3: viewing direction
  //   4: measurement block
  //
  ArrayOfVector grids(4);
  //
  if (no_pol_variation)
    grids[1] = Vector(1, 1);
  else
    grids[1] = Vector(0, sensor_response_pol_grid.nelem(), 1);
  if (no_los_variation)
    grids[2] = Vector(1, 1);
  else
    grids[2] = Vector(0, sensor_response_dlos_grid.nrows(), 1);
  if (no_mblock_variation)
    grids[3] = Vector(1, 1);
  else
    grids[3] = Vector(0, sensor_pos.nrows(), 1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(SINEFIT_MAINTAG);
  rq.Mode("");
  rq.Analytical(0);
  rq.Perturbation(0);

  // Each sinefit coeff. pair is treated as a retrieval quantity
  //
  for (Index i = 0; i < np; i++) {
    ostringstream sstr;
    sstr << "Period " << i;
    rq.Subtag(sstr.str());

    // "Grid" has length 2, set to period length
    grids[0] = Vector(2, period_lengths[i]);
    rq.Grids(grids);

    // Add it to the *jacobian_quantities*
    jq.push_back(rq);

    // Add pointing method to the jacobian agenda
    jacobian_agenda.append("jacobianCalcSinefit", i);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcSinefit(Matrix& jacobian,
                         const Index& mblock_index,
                         const Vector& iyb _U_,
                         const Vector& yb _U_,
                         const Sparse& sensor_response,
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Vector& sensor_response_f_grid,
                         const Matrix& sensor_response_dlos_grid,
                         const ArrayOfRetrievalQuantity& jacobian_quantities,
                         const Index& period_index,
                         const Verbosity&) {
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  ostringstream sstr;
  sstr << "Period " << period_index;
  for (iq = 0; iq < jacobian_quantities.nelem() && !found; iq++) {
    if (jacobian_quantities[iq].MainTag() == SINEFIT_MAINTAG &&
        jacobian_quantities[iq].Subtag() == sstr.str()) {
      found = true;
      break;
    }
  }
  if (!found) {
    throw runtime_error(
        "There is no Sinefit jacobian defined, in general "
        "or for the selected period length.\n");
  }

  // Size and check of sensor_response
  //
  const Index nf = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nlos = sensor_response_dlos_grid.nrows();

  // Make vectors with values to distribute over *jacobian*
  //
  // (period length stored in grid 0)
  ArrayOfVector jg = jacobian_quantities[iq].Grids();
  //
  Vector s(nf), c(nf);
  //
  for (Index f = 0; f < nf; f++) {
    Numeric a = (sensor_response_f_grid[f] - sensor_response_f_grid[0]) * 2 *
                PI / jg[0][0];
    s[f] = sin(a);
    c[f] = cos(a);
  }

  // Fill J
  //
  ArrayOfArrayOfIndex jacobian_indices;
  {
    bool any_affine;
    jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities, true);
  }
  //
  const Index n1 = jg[1].nelem();
  const Index n2 = jg[2].nelem();
  const Index n3 = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row4 = rowind.get_start();
  Index col4 = jacobian_indices[iq][0];

  if (n3 > 1) {
    col4 += mblock_index * n2 * n1 * 2;
  }

  for (Index l = 0; l < nlos; l++) {
    const Index row3 = row4 + l * nf * npol;
    const Index col3 = col4 + l * n1 * 2;

    for (Index f = 0; f < nf; f++) {
      const Index row2 = row3 + f * npol;

      for (Index p = 0; p < npol; p++) {
        Index col1 = col3;
        if (n1 > 1) {
          col1 += p * 2;
        }

        jacobian(row2 + p, col1) = s[f];
        jacobian(row2 + p, col1 + 1) = c[f];
      }
    }
  }
}

//----------------------------------------------------------------------------
// Surface quantities
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSurfaceQuantity(Workspace&,
                                ArrayOfRetrievalQuantity& jq,
                                Agenda& jacobian_agenda,
                                const Index& atmosphere_dim,
                                const Vector& lat_grid,
                                const Vector& lon_grid,
                                const Vector& rq_lat_grid,
                                const Vector& rq_lon_grid,
                                const String& quantity,
                                const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == SURFACE_MAINTAG && jq[it].Subtag() == quantity) {
      ostringstream os;
      os << quantity << " is already included as a surface variable "
         << "in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(max(atmosphere_dim - 1, Index(1)));
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               lat_grid,
                               lon_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(SURFACE_MAINTAG);
  rq.Subtag(quantity);
  rq.Analytical(0);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  // Add dummy
  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Temperatures (atmospheric):
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddTemperature(Workspace&,
                            ArrayOfRetrievalQuantity& jq,
                            Agenda& jacobian_agenda,
                            const Index& atmosphere_dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const Vector& rq_p_grid,
                            const Vector& rq_lat_grid,
                            const Vector& rq_lon_grid,
                            const String& hse,
                            const Verbosity& verbosity) {
  CREATE_OUT3;

  // Check that temperature is not already included in the jacobian.
  // We only check the main tag.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == TEMPERATURE_MAINTAG) {
      ostringstream os;
      os << "Temperature is already included in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Set subtag
  String subtag;
  if (hse == "on") {
    subtag = "HSE on";
  } else if (hse == "off") {
    subtag = "HSE off";
  } else {
    ostringstream os;
    os << "The keyword for hydrostatic equilibrium can only be set to\n"
       << "\"on\" or \"off\"\n";
    throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag(TEMPERATURE_MAINTAG);
  rq.Subtag(subtag);
  rq.Mode("abs");
  rq.Analytical(1);
  rq.Perturbation(0.1);
  rq.Grids(grids);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.PropType(JacPropMatType::Temperature);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Winds:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddWind(Workspace&,
                     ArrayOfRetrievalQuantity& jq,
                     Agenda& jacobian_agenda,
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Vector& rq_p_grid,
                     const Vector& rq_lat_grid,
                     const Vector& rq_lon_grid,
                     const String& component,
                     const Numeric& dfrequency,
                     const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == WIND_MAINTAG && jq[it].Subtag() == component) {
      ostringstream os;
      os << "The wind component:\n"
         << component << "\nis already included "
         << "in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;

  if (component == "u")
    rq.PropType(JacPropMatType::WindU);
  else if (component == "v")
    rq.PropType(JacPropMatType::WindV);
  else if (component == "w")
    rq.PropType(JacPropMatType::WindW);
  else if (component == "strength")
    rq.PropType(JacPropMatType::WindMagnitude);
  else
    throw std::runtime_error(
        "The selection for *component* can only be \"u\", \"v\", \"w\" or \"strength\".");

  rq.MainTag(WIND_MAINTAG);
  rq.Subtag(component);
  rq.Analytical(1);
  rq.Grids(grids);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.Perturbation(dfrequency);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  out3 << "  Calculations done by propagation matrix expression.\n";
  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Magnetic field:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddMagField(Workspace&,
                         ArrayOfRetrievalQuantity& jq,
                         Agenda& jacobian_agenda,
                         const Index& atmosphere_dim,
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Vector& rq_p_grid,
                         const Vector& rq_lat_grid,
                         const Vector& rq_lon_grid,
                         const String& component,
                         const Numeric& dB,
                         const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == MAGFIELD_MAINTAG && jq[it].Subtag() == component) {
      ostringstream os;
      os << "The magnetic field component:\n"
         << component << "\nis already "
         << "included in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  if (component == "u")
    rq.PropType(JacPropMatType::MagneticU);
  else if (component == "v")
    rq.PropType(JacPropMatType::MagneticV);
  else if (component == "w")
    rq.PropType(JacPropMatType::MagneticW);
  else if (component == "strength")
    rq.PropType(JacPropMatType::MagneticMagnitude);
  else
    throw runtime_error(
        "The selection for *component* can only be \"u\", \"v\", \"w\", or \"strength\".");

  rq.MainTag(MAGFIELD_MAINTAG);
  rq.Subtag(component);
  rq.Analytical(1);
  rq.Grids(grids);

  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.Perturbation(dB);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  // Add gas species method to the jacobian agenda
  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Catalog parameters:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddShapeCatalogParameter(Workspace&,
                                      ArrayOfRetrievalQuantity& jq,
                                      Agenda& jacobian_agenda,
                                      const QuantumIdentifier& line_identity,
                                      const String& species,
                                      const String& variable,
                                      const String& coefficient,
                                      const Verbosity& verbosity) {
  CREATE_OUT3;

  if (line_identity.Type() not_eq QuantumIdentifier::TRANSITION)
    throw std::runtime_error("Identity has to identify a line");

  const JacPropMatType jpt = select_derivativeLineShape(variable, coefficient);

  out3 << "Attempting to create RT tag for " << line_identity << " " << variable
       << " " << coefficient << " for ";
  if (species not_eq LineShape::self_broadening and
      species not_eq LineShape::bath_broadening)
    out3 << SpeciesTag(species).SpeciesNameMain() << "\n";
  else
    out3 << species << "\n";

  // Create the quantity
  RetrievalQuantity rq;
  rq.MainTag(CATALOGPARAMETER_MAINTAG);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.Mode(species);
  rq.Analytical(1);
  rq.Grids(ArrayOfVector(0, Vector(0)));
  rq.QuantumIdentity(line_identity);
  rq.PropType(jpt);
  rq.IntegrationOn();

  // Test this is not a copy
  for (auto& q : jq)
    if (q.HasSameInternalsAs(rq))
      throw std::runtime_error("Error with copies of the quantities");

  // Append and do housekeeping
  jq.push_back(rq);
  out3 << "Creation was successful!\n";
  jacobian_agenda.append("jacobianCalcDoNothing",
                         TokVal());  // old code activation
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddShapeCatalogParameters(
    Workspace& ws,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    const ArrayOfQuantumIdentifier& line_identities,
    const ArrayOfString& species,
    const ArrayOfString& variables,
    const ArrayOfString& coefficients,
    const Verbosity& verbosity) {
  if (not(line_identities.nelem() or species.nelem() or variables.nelem() or
          coefficients.nelem()))
    throw std::runtime_error("Must have at least 1-long lists for all GINs");

  ArrayOfString vars;
  if (variables[0] == "ALL")
    vars = AllLineShapeVars();
  else
    vars = variables;

  ArrayOfString coeffs;
  if (coefficients[0] == "ALL")
    coeffs = AllLineShapeCoeffs();
  else
    coeffs = coefficients;

  for (auto& l : line_identities)
    for (auto& s : species)
      for (auto& v : vars)
        for (auto& c : coeffs)
          jacobianAddShapeCatalogParameter(
              ws, jq, jacobian_agenda, l, s, v, c, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddBasicCatalogParameter(Workspace&,
                                      ArrayOfRetrievalQuantity& jq,
                                      Agenda& jacobian_agenda,
                                      const QuantumIdentifier& catalog_identity,
                                      const String& catalog_parameter,
                                      const Verbosity& verbosity) {
  CREATE_OUT3;

  // Check that this is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == CATALOGPARAMETER_MAINTAG &&
        jq[it].QuantumIdentity() == catalog_identity &&
        jq[it].Mode() == catalog_parameter) {
      ostringstream os;
      os << "The catalog identifier:\n"
         << catalog_identity << "\nis already included in "
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;

  // Check catalog_parameter here
  if (LINESTRENGTH_MODE == catalog_parameter)
    rq.PropType(JacPropMatType::LineStrength);
  else if (LINECENTER_MODE == catalog_parameter)
    rq.PropType(JacPropMatType::LineCenter);
  else {
    ostringstream os;
    os << "You have selected:\n"
       << catalog_parameter
       << "\nas your catalog parameter. This is not supported.\n"
       << "Please see user guide for supported parameters.\n";
    throw std::runtime_error(os.str());
  }

  rq.MainTag(CATALOGPARAMETER_MAINTAG);
  rq.Mode(catalog_parameter);
  rq.QuantumIdentity(catalog_identity);
  rq.Analytical(1);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.IntegrationOn();

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  out3 << "  Calculations done by propagation matrix expressions.\n";

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddBasicCatalogParameters(
    Workspace& ws,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    const ArrayOfQuantumIdentifier& catalog_identities,
    const ArrayOfString& catalog_parameters,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  out2 << " Adding " << catalog_identities.nelem() * catalog_parameters.nelem()
       << " expression(s) to the Jacobian calculations.\n";

  for (auto& qi : catalog_identities)
    for (auto& param : catalog_parameters)
      jacobianAddBasicCatalogParameter(
          ws, jq, jacobian_agenda, qi, param, verbosity);
}

//----------------------------------------------------------------------------
// NLTE temperature:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddNLTE(Workspace&,
                     ArrayOfRetrievalQuantity& jq,
                     Agenda& jacobian_agenda,
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Vector& rq_p_grid,
                     const Vector& rq_lat_grid,
                     const Vector& rq_lon_grid,
                     const QuantumIdentifier& energy_level_identity,
                     const Numeric& dx,
                     const Verbosity& verbosity) {
  CREATE_OUT3;

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it].MainTag() == NLTE_MAINTAG and
        jq[it].QuantumIdentity() == energy_level_identity) {
      ostringstream os;
      os << "The NLTE identifier:\n"
         << energy_level_identity << "\nis already included in "
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (not check_retrieval_grids(grids,
                                  os,
                                  p_grid,
                                  lat_grid,
                                  lon_grid,
                                  rq_p_grid,
                                  rq_lat_grid,
                                  rq_lon_grid,
                                  "retrieval pressure grid",
                                  "retrieval latitude grid",
                                  "retrievallongitude_grid",
                                  atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;

  rq.MainTag(NLTE_MAINTAG);
  rq.QuantumIdentity(energy_level_identity);
  rq.Perturbation(dx);
  rq.Grids(grids);
  rq.Analytical(1);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  out3 << "  Calculations done by propagation matrix expressions.\n";

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

void jacobianAddNLTEs(Workspace& ws,
                      ArrayOfRetrievalQuantity& jq,
                      Agenda& jacobian_agenda,
                      const Index& atmosphere_dim,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Vector& rq_p_grid,
                      const Vector& rq_lat_grid,
                      const Vector& rq_lon_grid,
                      const ArrayOfQuantumIdentifier& energy_level_identities,
                      const Numeric& dx,
                      const Verbosity& verbosity) {
  for (const auto& qi : energy_level_identities)
    jacobianAddNLTE(ws,
                    jq,
                    jacobian_agenda,
                    atmosphere_dim,
                    p_grid,
                    lat_grid,
                    lon_grid,
                    rq_p_grid,
                    rq_lat_grid,
                    rq_lon_grid,
                    qi,
                    dx,
                    verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSpecialSpecies(Workspace&,
                               ArrayOfRetrievalQuantity& jq,
                               Agenda& jacobian_agenda,
                               const Index& atmosphere_dim,
                               const Vector& p_grid,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const Vector& rq_p_grid,
                               const Vector& rq_lat_grid,
                               const Vector& rq_lon_grid,
                               const String& species,
                               const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if (!check_retrieval_grids(grids,
                               os,
                               p_grid,
                               lat_grid,
                               lon_grid,
                               rq_p_grid,
                               rq_lat_grid,
                               rq_lon_grid,
                               "retrieval pressure grid",
                               "retrieval latitude grid",
                               "retrievallongitude_grid",
                               atmosphere_dim))
      throw runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Grids(grids);
  rq.Analytical(1);
  rq.SubSubtag(PROPMAT_SUBSUBTAG);

  // Make sure modes are valid and complain if they are repeated
  if (species == "electrons") {
    for (Index it = 0; it < jq.nelem(); it++) {
      if (jq[it].MainTag() == ELECTRONS_MAINTAG) {
        ostringstream os;
        os << "Electrons are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.MainTag(ELECTRONS_MAINTAG);
    rq.PropType(JacPropMatType::Electrons);
  } else if (species == "particulates") {
    for (Index it = 0; it < jq.nelem(); it++) {
      if (jq[it].MainTag() == PARTICULATES_MAINTAG) {
        ostringstream os;
        os << "Particulates are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.MainTag(PARTICULATES_MAINTAG);
    rq.PropType(JacPropMatType::Particulates);
  } else {
    ostringstream os;
    os << "Unknown special species jacobian: \"" << species
       << "\"\nPlease see *jacobianAddSpecialSpecies* for viable options.";
    throw std::runtime_error(os.str());
  }

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.append("jacobianCalcDoNothing", TokVal());
}

//----------------------------------------------------------------------------
// Adjustments and transformations
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAdjustAndTransform(
    Matrix& jacobian,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& x,
    const Verbosity&) {
  // For flexibility inside inversion_iteration_agenda, we should accept empty
  // Jacobian
  if (jacobian.empty()) {
    return;
  }

  // Adjustments
  //
  // Unfortunately the adjustment requires both range indices and the
  // untransformed x, which makes things a bit messy
  bool vars_init = false;
  ArrayOfArrayOfIndex jis0;
  Vector x0;
  //
  for (Index q = 0; q < jacobian_quantities.nelem(); q++) {
    if (jacobian_quantities[q].MainTag() == ABSSPECIES_MAINTAG &&
        jacobian_quantities[q].Mode() == "rel") {
      if (!vars_init) {
        bool any_affine;
        jac_ranges_indices(jis0, any_affine, jacobian_quantities, true);
        x0 = x;
        transform_x_back(x0, jacobian_quantities);
        vars_init = true;
      }
      for (Index i = jis0[q][0]; i <= jis0[q][1]; i++) {
        if (x[i] != 1) {
          jacobian(joker, i) /= x[i];
        }
      }
    }
  }

  // Transformations
  transform_jacobian(jacobian, x, jacobian_quantities);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianSetAffineTransformation(ArrayOfRetrievalQuantity& jqs,
                                     const Matrix& transformation_matrix,
                                     const Vector& offset_vector,
                                     const Verbosity& /*v*/
) {
  if (jqs.empty()) {
    runtime_error(
        "Jacobian quantities is empty, so there is nothing to add the "
        "transformation to.");
  }

  Index nelem = jqs.back().Grids().nelem();

  if (!(nelem == transformation_matrix.nrows())) {
    runtime_error(
        "Dimension of transformation matrix incompatible with retrieval grids.");
  }
  if (!(nelem == offset_vector.nelem())) {
    runtime_error(
        "Dimension of offset vector incompatible with retrieval grids.");
  }

  jqs.back().SetTransformationMatrix(transpose(transformation_matrix));
  jqs.back().SetOffsetVector(offset_vector);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianSetFuncTransformation(ArrayOfRetrievalQuantity& jqs,
                                   const String& transformation_func,
                                   const Numeric& z_min,
                                   const Numeric& z_max,
                                   const Verbosity& /*v*/
) {
  if (jqs.empty())
    throw runtime_error(
        "Jacobian quantities is empty, so there is nothing to add the "
        "transformation to.");

  if (transformation_func == "none") {
    jqs.back().SetTransformationFunc("");
    return;
  }

  Vector pars;

  if (transformation_func == "atanh") {
    if (z_max <= z_min)
      throw runtime_error(
          "For option atanh, the GIN *z_max* must be set and be > z_min.");
    pars.resize(2);
    pars[0] = z_min;
    pars[1] = z_max;
  } else if (transformation_func == "log" || transformation_func == "log10") {
    pars.resize(1);
    pars[0] = z_min;
  } else {
    ostringstream os;
    os << "Valid options for *transformation_func* are:\n"
       << "\"none\", \"log\", \"log10\" and \"atanh\"\n"
       << "But found: \"" << transformation_func << "\"";
    throw runtime_error(os.str());
  }

  jqs.back().SetTransformationFunc(transformation_func);
  jqs.back().SetTFuncParameters(pars);
}

//----------------------------------------------------------------------------
// Methods for doing perturbations
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPerturb(Tensor3& perturbed_field,
                    const Index& atmosphere_dim,
                    const Vector& p_grid,
                    const Vector& lat_grid,
                    const Vector& lon_grid,
                    const Tensor3& original_field,
                    const Vector& p_ret_grid,
                    const Vector& lat_ret_grid,
                    const Vector& lon_ret_grid,
                    const Index& pert_index,
                    const Numeric& pert_size,
                    const String& pert_mode,
                    const Verbosity&) {
  // Input checks (more below)
  chk_atm_field("original_field",
                original_field,
                atmosphere_dim,
                p_grid,
                lat_grid,
                lon_grid,
                false );

  // Pack retrieval grids into an ArrayOfVector
  ArrayOfVector ret_grids(atmosphere_dim);
  ret_grids[0] = p_ret_grid;
  if (atmosphere_dim>1){
    ret_grids[1] = lat_ret_grid;
    if (atmosphere_dim>2){
      ret_grids[2] = lon_ret_grid;
    }
  }

  // Find mapping from retrieval grids to atmospheric grids
  ArrayOfGridPos gp_p, gp_lat, gp_lon;
  Index n_p, n_lat, n_lon;
  get_gp_rq_to_atmgrids(gp_p,
                        gp_lat,
                        gp_lon,
                        n_p,
                        n_lat,
                        n_lon,
                        ret_grids,
                        atmosphere_dim,
                        p_grid,
                        lat_grid,
                        lon_grid);

  // Now we can chec *pert_index*
  if (pert_index<0){
    throw runtime_error("Bad *pert_index*. It is negative.");
  }
  const Index n_tot = n_p * n_lat * n_lon;
  if (pert_index >= n_tot){
    throw runtime_error("Bad *pert_index*. It is too high with respect "
                        "to length of retrieval grids.");
  }    
  
  // Create x-vector that matches perturbation
  Vector x(n_tot);
  if (pert_mode == "absolute" ){
    x = 0;
    x[pert_index] = pert_size;
  }
  else if (pert_mode == "relative" ){
    x = 1;
    x[pert_index] += pert_size;
  }
  else{
    throw runtime_error("Bad *pert_mode*. Allowed choices are: "
                        """absolute"" and ""relative"".");
  }
  
  // Map x to a perturbation defined at atmospheric grids
  Tensor3 x3d(n_p, n_lat, n_lon), pert(n_p, n_lat, n_lon);
  reshape(x3d, x);
  regrid_atmfield_by_gp_oem(pert, atmosphere_dim, x3d, gp_p, gp_lat, gp_lon);
  
  // Init perturbed_field, if not equal to original_field
  if (&perturbed_field != &original_field) {
    perturbed_field = original_field;
  }

  // Apply perturbation
  if (pert_mode == "absolute" ){
    perturbed_field += pert;
  }
  else{
    perturbed_field *= pert;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPerturbAtmGrids(Tensor3& perturbed_field,
                             const Index& atmosphere_dim,
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Tensor3& original_field,
                             const Index& pert_index,
                             const Numeric& pert_size,
                             const String& pert_mode,
                             const Verbosity&) {
  // Some sizes
  const Index n_p = p_grid.nelem();
  const Index n_lat = atmosphere_dim<2 ? 1 : lat_grid.nelem();
  const Index n_lon = atmosphere_dim<3 ? 1 : lon_grid.nelem();
  
  // Check input
  chk_atm_field("original_field",
                original_field,
                atmosphere_dim,
                p_grid,
                lat_grid,
                lon_grid,
                false );
  if (pert_index<0){
    throw runtime_error("Bad *pert_index*. It is negative.");
  }
  if (pert_index >= n_p * n_lat * n_lon){
    throw runtime_error("Bad *pert_index*. It is too high with respect "
                        "to length of atmospheric grids.");
  }    

  // Determine indexes with respect to atmospheric grids
  Index tot_index = pert_index;
  const Index lon_index = atmosphere_dim<3 ? 0 : tot_index / (n_lat * n_p);
  tot_index -= lon_index * n_lat * n_p;
  const Index lat_index = atmosphere_dim<2 ? 0 : tot_index / n_p;
  tot_index -= lat_index * n_p;
  const Index p_index = tot_index;
  
  // Init perturbed_field, if not equal to original_field
  if (&perturbed_field != &original_field) {
    perturbed_field = original_field;
  }

  // Perturb
  if (pert_mode == "absolute" ){
    perturbed_field(p_index,
                    atmosphere_dim>1 ? lat_index : 0,
                    atmosphere_dim>2 ? lon_index : 0) += pert_size;
  }
  else if (pert_mode == "relative"){
    perturbed_field(p_index,
                    atmosphere_dim>1 ? lat_index : 0,
                    atmosphere_dim>2 ? lon_index : 0) *= 1 + pert_size;
  }
  else{
    throw runtime_error("Bad *pert_mode*. Allowed choices are: "
                        """absolute"" and ""relative"".");
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexNumberOfAtmosphericPoints(Index& n,
                                    const Index& atmosphere_dim,
                                    const Vector& p_grid,
                                    const Vector& lat_grid,
                                    const Vector& lon_grid,
                                    const Verbosity&) {
  const Index n_p = p_grid.nelem();
  const Index n_lat = atmosphere_dim<2 ? 1 : lat_grid.nelem();
  const Index n_lon = atmosphere_dim<3 ? 1 : lon_grid.nelem();

  n = n_p * n_lat * n_lon;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianFromTwoY(Matrix& jacobian,
                    const Vector& y_pert,
                    const Vector& y,
                    const Numeric& pert_size,
                    const Verbosity&) {
  const Index n = y.nelem();
  if( y_pert.nelem() != n ){
    throw runtime_error("Inconsistency in length of *y_pert* and *y*.");
  }
  jacobian = y_pert;
  jacobian -= y;
  jacobian /= pert_size;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianFromYbatch(Matrix& jacobian,
                    const ArrayOfVector& ybatch,
                    const Vector& y,
                    const Numeric& pert_size,
                    const Verbosity&) {
  const Index n = y.nelem();
  const Index l = ybatch.nelem();
  if (l>0){
      if( ybatch[0].nelem() != n )
        throw runtime_error("Inconsistency in length of y and ybatch[0].");
    }
  jacobian.resize(n,l);
  for (Index i=0; i<l; i++) {
    jacobian(joker,i) = ybatch[i];
    jacobian(joker,i) -= y;
  }
  jacobian /= pert_size;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_bulkprop_fieldPerturb(Tensor4& particle_bulkprop_field,
                                    const Index& atmosphere_dim,
                                    const Vector& p_grid,
                                    const Vector& lat_grid,
                                    const Vector& lon_grid,
                                    const ArrayOfString& particle_bulkprop_names,
                                    const String& particle_type,
                                    const Vector& p_ret_grid,
                                    const Vector& lat_ret_grid,
                                    const Vector& lon_ret_grid,
                                    const Index& pert_index,
                                    const Numeric& pert_size,
                                    const String& pert_mode,
                                    const Verbosity& verbosity) {
  // Locate particle_type among particle_bulkprop_names
  Index iq = find_first(particle_bulkprop_names, particle_type);
  if (iq < 0) {
    ostringstream os;
    os << "Could not find " << particle_type << " in *particle_bulkprop_names*.\n";
    throw std::runtime_error(os.str());
  }

  Tensor3 original_field, perturbed_field;
  original_field = particle_bulkprop_field(iq,joker,joker,joker);
  AtmFieldPerturb(perturbed_field,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid,
                  original_field,
                  p_ret_grid,
                  lat_ret_grid,
                  lon_ret_grid,
                  pert_index,
                  pert_size,
                  pert_mode,
                  verbosity);
  particle_bulkprop_field(iq,joker,joker,joker) = perturbed_field;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_bulkprop_fieldPerturbAtmGrids(Tensor4& particle_bulkprop_field,
                                            const Index& atmosphere_dim,
                                            const Vector& p_grid,
                                            const Vector& lat_grid,
                                            const Vector& lon_grid,
                                            const ArrayOfString& particle_bulkprop_names,
                                            const String& particle_type,
                                            const Index& pert_index,
                                            const Numeric& pert_size,
                                            const String& pert_mode,
                                            const Verbosity& verbosity) {
  // Locate particle_type among particle_bulkprop_names
  Index iq = find_first(particle_bulkprop_names, particle_type);
  if (iq < 0) {
    ostringstream os;
    os << "Could not find " << particle_type << " in *particle_bulkprop_names*.\n";
    throw std::runtime_error(os.str());
  }

  Tensor3 original_field, perturbed_field;
  original_field = particle_bulkprop_field(iq,joker,joker,joker);
  AtmFieldPerturbAtmGrids(perturbed_field,
                          atmosphere_dim,
                          p_grid,
                          lat_grid,
                          lon_grid,
                          original_field,
                          pert_index,
                          pert_size,
                          pert_mode,
                          verbosity);
  particle_bulkprop_field(iq,joker,joker,joker) = perturbed_field;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldPerturb(Tensor4& vmr_field,
                      const Index& atmosphere_dim,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const String& species,
                      const Vector& p_ret_grid,
                      const Vector& lat_ret_grid,
                      const Vector& lon_ret_grid,
                      const Index& pert_index,
                      const Numeric& pert_size,
                      const String& pert_mode,
                      const Verbosity& verbosity) {
  // Locate vmr_species among abs_species
  Index iq = -1;
  for (Index i = 0; i < abs_species.nelem(); i++) {
    if (abs_species[i][0].Species() == SpeciesTag(species).Species()) {
      iq = i;
      break;
    }
  }
  if (iq < 0) {
    ostringstream os;
    os << "Could not find " << species << " in *abs_species*.\n";
    throw std::runtime_error(os.str());
  }

  Tensor3 original_field, perturbed_field;
  original_field = vmr_field(iq,joker,joker,joker);
  AtmFieldPerturb(perturbed_field,
                  atmosphere_dim,
                  p_grid,
                  lat_grid,
                  lon_grid,
                  original_field,
                  p_ret_grid,
                  lat_ret_grid,
                  lon_ret_grid,
                  pert_index,
                  pert_size,
                  pert_mode,
                  verbosity);
  vmr_field(iq,joker,joker,joker) = perturbed_field;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldPerturbAtmGrids(Tensor4& vmr_field,
                              const Index& atmosphere_dim,
                              const Vector& p_grid,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const ArrayOfArrayOfSpeciesTag& abs_species,
                              const String& species,
                              const Index& pert_index,
                              const Numeric& pert_size,
                              const String& pert_mode,
                              const Verbosity& verbosity) {
  // Locate vmr_species among abs_species
  Index iq = -1;
  for (Index i = 0; i < abs_species.nelem(); i++) {
    if (abs_species[i][0].Species() == SpeciesTag(species).Species()) {
      iq = i;
      break;
    }
  }
  if (iq < 0) {
    ostringstream os;
    os << "Could not find " << species << " in *abs_species*.\n";
    throw std::runtime_error(os.str());
  }

  Tensor3 original_field, perturbed_field;
  original_field = vmr_field(iq,joker,joker,joker);
  AtmFieldPerturbAtmGrids(perturbed_field,
                          atmosphere_dim,
                          p_grid,
                          lat_grid,
                          lon_grid,
                          original_field,
                          pert_index,
                          pert_size,
                          pert_mode,
                          verbosity);
  vmr_field(iq,joker,joker,joker) = perturbed_field;
}


