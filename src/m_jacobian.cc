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
#include <workspace.h>
#include "check_input.h"
#include "cloudbox.h"
#include "arts_conversions.h"
#include "interp.h"
#include "jacobian.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "matpack_math.h"
#include "physics_funcs.h"
#include "rte.h"
#include "special_interp.h"

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
                           const Vector& yb _U_) {
  /* Nothing to do here for the analytical case, this function just exists
   to satisfy the required inputs and outputs of the jacobian_agenda */
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianClose(Index& jacobian_do,
                   Agenda& jacobian_agenda,
                   const ArrayOfRetrievalQuantity& jacobian_quantities) {
  // Make sure that the array is not empty
  if (jacobian_quantities.empty())
    throw std::runtime_error(
        "No retrieval quantities has been added to *jacobian_quantities*.");

  jacobian_agenda.finalize();
  jacobian_do = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianInit(ArrayOfRetrievalQuantity& jacobian_quantities,
                  Agenda& jacobian_agenda) {
  jacobian_quantities.resize(0);
  jacobian_agenda = Agenda("jacobian_agenda");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianOff(Index& jacobian_do,
                 Agenda& jacobian_agenda,
                 ArrayOfRetrievalQuantity& jacobian_quantities) {
  jacobian_do = 0;
  jacobianInit(jacobian_quantities, jacobian_agenda);
}

//----------------------------------------------------------------------------
// Absorption species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddAbsSpecies(const Workspace&,
                           ArrayOfRetrievalQuantity& jq,
                           Agenda& jacobian_agenda,
                           const Vector& rq_p_grid,
                           const Vector& rq_lat_grid,
                           const Vector& rq_lon_grid,
                           const String& species,
                           const String& mode,
                           const Index& for_species_tag) {
  QuantumIdentifier qi;
  if (not for_species_tag) {
    ArrayOfSpeciesTag test(species);
    if (test.nelem() not_eq 1)
      throw std::runtime_error(
          "Trying to add a species as a species tag of multiple species.\n"
          "This is not supported.  Please give just a single species instead.\n"
          "Otherwise consider if you intended for_species_tag to be evaluated true.\n");
    qi = QuantumIdentifier(test[0].Isotopologue());
  }

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    ARTS_USER_ERROR_IF (jq[it] == Jacobian::Special::ArrayOfSpeciesTagVMR && jq[it].Subtag() == species,
      "The gas species:\n", species, "\nis already included in *jacobian_quantities*.")
    ARTS_USER_ERROR_IF (jq[it] == Jacobian::Line::VMR and jq[it].QuantumIdentity() == qi,
      "The atmospheric species of:\n", species, "\nis already included in *jacobian_quantities*\n"
      "as: ", jq[it])
  }

  // Check that mode is correct
  if (mode != "vmr" && mode != "nd" && mode != "rel" && mode != "rh" &&
      mode != "q") {
    throw std::runtime_error(
        "The retrieval mode can only be \"vmr\", \"nd\", "
        "\"rel\", \"rh\" or \"q\".");
  }
  if ((mode == "rh" || mode == "q") && species.substr(0, 3) != "H2O") {
    throw std::runtime_error(
        "Retrieval modes \"rh\" and \"q\" can only be applied "
        "on species starting with H2O.");
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Subtag(species);
  rq.Mode(mode);
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});
  if (for_species_tag == 0) {
    rq.Target(Jacobian::Target(Jacobian::Line::VMR, qi, qi.Species()));
  } else {
    ArrayOfSpeciesTag test(species);
    rq.Target(Jacobian::Target(Jacobian::Special::ArrayOfSpeciesTagVMR, test));
  }
  rq.Target().perturbation = 0.001;

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add(Method{"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Frequency shift
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqShift(const Workspace& ws _U_,
                          ArrayOfRetrievalQuantity& jacobian_quantities,
                          Agenda& jacobian_agenda,
                          const Vector& f_grid,
                          const Numeric& df) {
  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it] == Jacobian::Sensor::FrequencyShift) {
      std::ostringstream os;
      os << "Fit of frequency shift is already included in\n"
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Checks of frequencies
  if (df <= 0) throw std::runtime_error("The argument *df* must be > 0.");
  if (df > 1e6)
    throw std::runtime_error("The argument *df* is not allowed to exceed 1 MHz.");
  const Index nf = f_grid.nelem();
  if (nf < 2)
    throw std::runtime_error(
        "Frequency shifts and *f_grid* of length 1 can "
        "not be combined.");
  const Numeric maxdf = f_grid[nf - 1] - f_grid[nf - 2];
  if (df > maxdf) {
    std::ostringstream os;
    os << "The value of *df* is too big with respect to spacing of "
       << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
       << "between the two last elements of *f_grid*.\n"
       << "This spacing is   : " << maxdf / 1e3 << " kHz\n"
       << "The value of df is: " << df / 1e3 << " kHz";
    throw std::runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Sensor::FrequencyShift);
  rq.Mode("");
  rq.Target().perturbation = df;

  // Dummy vector of length 1
  Vector grid(1, 0);
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.add(Method{"jacobianCalcFreqShift"});
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqShift(Matrix& jacobian,
                           const Index& mblock_index,
                           const Vector& iyb,
                           const Vector& yb,
                           const Vector& f_grid,
                           const Matrix& mblock_dlos,
                           const Sparse& sensor_response,
                           const ArrayOfRetrievalQuantity& jacobian_quantities) {
  // Set some useful (and needed) variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n] == Jacobian::Sensor::FrequencyShift) {
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
    throw std::runtime_error(
        "There is no such frequency retrieval quantity defined.\n");
  }

  // Check that sensor_response is consistent with yb and iyb
  //
  if (sensor_response.nrows() != yb.nelem())
    throw std::runtime_error("Mismatch in size between *sensor_response* and *yb*.");
  if (sensor_response.ncols() != iyb.nelem())
    throw std::runtime_error(
        "Mismatch in size between *sensor_response* and *iyb*.");

  // Get disturbed (part of) y
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    const Index nf2 = f_grid.nelem();
    const Index nlos2 = mblock_dlos.nrows();
    const Index niyb = nf2 * nlos2 * 4;

    // Interpolation weights
    //
    constexpr Index porder = 3;
    //
    Vector fg_new = f_grid, iyb2(niyb);
    fg_new += rq.Target().perturbation;
    //
    const auto lag = my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<porder>>(fg_new, f_grid);
    const auto itw = interpweights(lag);

    // Do interpolation
    for (Index ilos = 0; ilos < nlos2; ilos++) {
      const Index row0 = ilos * nf2 * 4;

      for (Index is = 0; is < 4; is++) {
        iyb2[Range(row0 + is, nf2, 4)] = reinterp(iyb[Range(row0 + is, nf2, 4)], itw, lag);
      }
    }

    // Determine difference
    //
    mult(dy, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (dy[i] - yb[i]) / rq.Target().perturbation;
    }
  }

  //--- Set jacobian ---
  ARTS_ASSERT(rq.Grids()[0].nelem() == 1);
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  jacobian(rowind, ji[0]) = dy;
}

//----------------------------------------------------------------------------
// Frequency stretch
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqStretch(const Workspace& ws _U_,
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            Agenda& jacobian_agenda,
                            const Vector& f_grid,
                            const Numeric& df) {
  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it] == Jacobian::Sensor::FrequencyStretch) {
      std::ostringstream os;
      os << "Fit of frequency stretch is already included in\n"
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Checks of df
  if (df <= 0) throw std::runtime_error("The argument *df* must be > 0.");
  if (df > 1e6)
    throw std::runtime_error("The argument *df* is not allowed to exceed 1 MHz.");
  const Index nf = f_grid.nelem();
  const Numeric maxdf = f_grid[nf - 1] - f_grid[nf - 2];
  if (df > maxdf) {
    std::ostringstream os;
    os << "The value of *df* is too big with respect to spacing of "
       << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
       << "between the two last elements of *f_grid*.\n"
       << "This spacing is   : " << maxdf / 1e3 << " kHz\n"
       << "The value of df is: " << df / 1e3 << " kHz";
    throw std::runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Sensor::FrequencyStretch);
  rq.Mode("");
  rq.Target().perturbation=df;

  // Dummy vector of length 1
  Vector grid(1, 0);
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.add(Method{"jacobianCalcFreqStretch"});
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqStretch(
    Matrix& jacobian,
    const Index& mblock_index,
    const Vector& iyb,
    const Vector& yb,
    const Vector& f_grid,
    const Matrix& mblock_dlos,
    const Sparse& sensor_response,
    const ArrayOfIndex& sensor_response_pol_grid,
    const Vector& sensor_response_f_grid,
    const Matrix& sensor_response_dlos_grid,
    const ArrayOfRetrievalQuantity& jacobian_quantities) {
  // The code here is close to identical to the one for Shift. The main
  // difference is that dy is weighted with poly_order 1 basis function.

  // Set some useful (and needed) variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n] == Jacobian::Sensor::FrequencyStretch) {
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
    throw std::runtime_error(
        "There is no such frequency retrieval quantity defined.\n");
  }

  // Check that sensor_response is consistent with yb and iyb
  //
  if (sensor_response.nrows() != yb.nelem())
    throw std::runtime_error("Mismatch in size between *sensor_response* and *yb*.");
  if (sensor_response.ncols() != iyb.nelem())
    throw std::runtime_error(
        "Mismatch in size between *sensor_response* and *iyb*.");

  // Get disturbed (part of) y
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    const Index nf2 = f_grid.nelem();
    const Index nlos2 = mblock_dlos.nrows();
    const Index niyb = nf2 * nlos2 * 4;

    // Interpolation weights
    //
    constexpr Index porder = 3;
    //
    Vector fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Target().perturbation;
    //
    const auto lag = my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<porder>>(fg_new, f_grid);
    const auto itw = interpweights(lag);

    // Do interpolation
    for (Index ilos = 0; ilos < nlos2; ilos++) {
      const Index row0 = ilos * nf2 * 4;

      for (Index is = 0; is < 4; is++) {
        iyb2[Range(row0 + is, nf2, 4)] = reinterp(iyb[Range(row0 + is, nf2, 4)], itw, lag);
      }
    }

    // Determine difference
    //
    mult(dy, sensor_response, iyb2);
    //
    for (Index i = 0; i < n1y; i++) {
      dy[i] = (dy[i] - yb[i]) / rq.Target().perturbation;
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
  ARTS_ASSERT(rq.Grids()[0].nelem() == 1);
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  jacobian(rowind, ji[0]) = dy;
}

//----------------------------------------------------------------------------
// Pointing:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPointingZa(const Workspace& ws _U_,
                           ArrayOfRetrievalQuantity& jacobian_quantities,
                           Agenda& jacobian_agenda,
                           const Matrix& sensor_pos,
                           const ArrayOfTime& sensor_time,
                           const Index& poly_order,
                           const String& calcmode,
                           const Numeric& dza) {
  // Check that poly_order is -1 or positive
  if (poly_order < -1)
    throw std::runtime_error(
        "The polynomial order has to be positive or -1 for gitter.");

  // Check that this jacobian type is not already included.
  for (Index it = 0; it < jacobian_quantities.nelem(); it++) {
    if (jacobian_quantities[it].Target().isPointing()) {
      std::ostringstream os;
      os << "Fit of zenith angle pointing off-set is already included in\n"
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Checks of dza
  if (dza <= 0) throw std::runtime_error("The argument *dza* must be > 0.");
  if (dza > 0.1)
    throw std::runtime_error("The argument *dza* is not allowed to exceed 0.1 deg.");

  // Check that sensor_time is consistent with sensor_pos
  if (sensor_time.nelem() != sensor_pos.nrows()) {
    std::ostringstream os;
    os << "The WSV *sensor_time* must be defined for every "
       << "measurement block.\n";
    throw std::runtime_error(os.str());
  }

  // Do not allow that *poly_order* is not too large compared to *sensor_time*
  if (poly_order > sensor_time.nelem() - 1) {
    throw std::runtime_error(
        "The polynomial order can not be >= length of *sensor_time*.");
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  if (calcmode == "recalc") {
    rq.Target() = Jacobian::Target(Jacobian::Sensor::PointingZenithRecalc);
    jacobian_agenda.add({"jacobianCalcPointingZaRecalc"});
  } else if (calcmode == "interp") {
    rq.Target() = Jacobian::Target(Jacobian::Sensor::PointingZenithInterp);
    jacobian_agenda.add({"jacobianCalcPointingZaInterp"});
  } else
    throw std::runtime_error(
      R"(Possible choices for *calcmode* are "recalc" and "interp".)");
  rq.Target().perturbation = dza;

  // To store the value or the polynomial order, create a vector with length
  // poly_order+1, in case of gitter set the size of the grid vector to be the
  // number of measurement blocks, all elements set to -1.
  Vector grid=uniform_grid(0, poly_order + 1, 1);
  if (poly_order == -1) {
    grid.resize(sensor_pos.nrows());
    grid = -1.0;
  }
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back(rq);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaInterp(
    Matrix& jacobian,
    const Index& mblock_index,
    const Vector& iyb,
    const Vector& yb _U_,
    const Vector& f_grid,
    const Matrix& DEBUG_ONLY(sensor_los),
    const Matrix& mblock_dlos,
    const Sparse& sensor_response,
    const ArrayOfTime& sensor_time,
    const ArrayOfRetrievalQuantity& jacobian_quantities) {
  if (mblock_dlos.nrows() < 2)
    throw std::runtime_error(
        "The method demands that *mblock_dlos* has "
        "more than one row.");

  if (!(is_increasing(mblock_dlos(joker, 0)) ||
        is_decreasing(mblock_dlos(joker, 0))))
    throw std::runtime_error(
        "The method demands that the zenith angles in "
        "*mblock_dlos* are sorted (increasing or decreasing).");

  // Set some useful variables.
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for (Index n = 0; n < jacobian_quantities.nelem() && !found; n++) {
    if (jacobian_quantities[n] == Jacobian::Sensor::PointingZenithInterp) {
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
    throw std::runtime_error(
        "There is no such pointing retrieval quantity defined.\n");
  }

  // Get "dy", by inter/extra-polation of existing iyb
  //
  const Index n1y = sensor_response.nrows();
  Vector dy(n1y);
  {
    // Sizes
    const Index nf = f_grid.nelem();
    const Index nza = mblock_dlos.nrows();

    // Shifted zenith angles
    Vector za1{mblock_dlos(joker, 0)};
    za1 -= rq.Target().perturbation;
    Vector za2{mblock_dlos(joker, 0)};
    za2 += rq.Target().perturbation;

    // Find interpolation weights
    ArrayOfGridPos gp1(nza), gp2(nza);
    gridpos(
        gp1, mblock_dlos(joker, 0), za1, 1e6);  // Note huge extrapolation!
    gridpos(
        gp2, mblock_dlos(joker, 0), za2, 1e6);  // Note huge extrapolation!
    Matrix itw1(nza, 2), itw2(nza, 2);
    interpweights(itw1, gp1);
    interpweights(itw2, gp2);

    // Make interpolation (for all azimuth angles, frequencies and Stokes)
    //
    Vector iyb1(iyb.nelem()), iyb2(iyb.nelem());
    //
    for (Index iza = 0; iza < nza; iza++) {
      for (Index iv = 0; iv < nf; iv++) {
        for (Index is = 0; is < 4; is++) {
          const Range r(iv * 4 + is, nza, nf * 4);
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
      dy[i] = (y2[i] - y1[i]) / (2 * rq.Target().perturbation);
    }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
  const Index row0 = rowind.offset;

  // Handle pointing "jitter" seperately
  if (rq.Grids()[0][0] == -1)          // Not all values are set here,
  {                                    // but should already have been
    ARTS_ASSERT(lg == sensor_los.nrows());  // set to 0
    ARTS_ASSERT(rq.Grids()[0][mblock_index] == -1);
    jacobian(rowind, it + mblock_index) = dy;
  }

  // Polynomial representation
  else {
    Vector w;
    for (Index c = 0; c < lg; c++) {
      ARTS_ASSERT(Numeric(c) == rq.Grids()[0][c]);
      //
      polynomial_basis_func(w, time_vector(sensor_time), c);
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
void jacobianAddPolyfit(const Workspace& ws _U_,
                        ArrayOfRetrievalQuantity& jq,
                        Agenda& jacobian_agenda,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        const Matrix& sensor_response_dlos_grid,
                        const Matrix& sensor_pos,
                        const Index& poly_order,
                        const Index& no_pol_variation,
                        const Index& no_los_variation,
                        const Index& no_mblock_variation) {
  // Check that poly_order is >= 0
  if (poly_order < 0)
    throw std::runtime_error("The polynomial order has to be >= 0.");

  // Check that polyfit is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Sensor::Polyfit) {
      std::ostringstream os;
      os << "Polynomial baseline fit is already included in\n"
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
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
    grids[1] = uniform_grid(0, sensor_response_pol_grid.nelem(), 1);
  if (no_los_variation)
    grids[2] = Vector(1, 1);
  else
    grids[2] = uniform_grid(0, sensor_response_dlos_grid.nrows(), 1);
  if (no_mblock_variation)
    grids[3] = Vector(1, 1);
  else
    grids[3] = uniform_grid(0, sensor_pos.nrows(), 1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Sensor::Polyfit);
  rq.Mode("");
  rq.Target().perturbation = 0;

  // Each polynomial coeff. is treated as a retrieval quantity
  //
  for (Index i = 0; i <= poly_order; i++) {
    std::ostringstream sstr;
    sstr << "Coefficient " << i;
    rq.Subtag(sstr.str());

    // Grid is a scalar, use polynomial coeff.
    grids[0] = Vector(1, (Numeric)i);
    rq.Grids(grids);

    // Add it to the *jacobian_quantities*
    jq.push_back(rq);

    // Add pointing method to the jacobian agenda
    jacobian_agenda.add(Method{"my_poly_coeff", i});
    jacobian_agenda.add(Method{"jacobianCalcPolyfit", {}, {{"poly_coeff", "my_poly_coeff"}}});
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
                         const Index& poly_coeff) {
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  std::ostringstream sstr;
  sstr << "Coefficient " << poly_coeff;
  for (iq = 0; iq < jacobian_quantities.nelem() && !found; iq++) {
    if (jacobian_quantities[iq] == Jacobian::Sensor::Polyfit &&
        jacobian_quantities[iq].Subtag() == sstr.str()) {
      found = true;
      break;
    }
  }
  if (!found) {
    throw std::runtime_error(
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
  const Index row4 = rowind.offset;
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
void jacobianAddScatSpecies(const Workspace&,
                            ArrayOfRetrievalQuantity& jq,
                            Agenda& jacobian_agenda,
                            const Vector& rq_p_grid,
                            const Vector& rq_lat_grid,
                            const Vector& rq_lon_grid,
                            const String& species,
                            const String& quantity) {
  // Check that this species+quantity combination is not already included in
  // the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Special::ScatteringString && jq[it].Subtag() == species &&
        jq[it].SubSubtag() == quantity) {
      std::ostringstream os;
      os << "The combintaion of\n   scattering species: " << species
         << "\n   retrieval quantity: " << quantity
         << "\nis already included in *jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Special::ScatteringString, species);
  rq.Subtag(species);
  rq.SubSubtag(quantity);
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Sinusoidal baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSinefit(const Workspace& ws _U_,
                        ArrayOfRetrievalQuantity& jq,
                        Agenda& jacobian_agenda,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        const Matrix& sensor_response_dlos_grid,
                        const Matrix& sensor_pos,
                        const Vector& period_lengths,
                        const Index& no_pol_variation,
                        const Index& no_los_variation,
                        const Index& no_mblock_variation) {
  const Index np = period_lengths.nelem();

  // Check that poly_order is >= 0
  if (np == 0) throw std::runtime_error("No sinusoidal periods has benn given.");

  // Check that polyfit is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Sensor::Sinefit) {
      std::ostringstream os;
      os << "Polynomial baseline fit is already included in\n"
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
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
    grids[1] = uniform_grid(0, sensor_response_pol_grid.nelem(), 1);
  if (no_los_variation)
    grids[2] = Vector(1, 1);
  else
    grids[2] = uniform_grid(0, sensor_response_dlos_grid.nrows(), 1);
  if (no_mblock_variation)
    grids[3] = Vector(1, 1);
  else
    grids[3] = uniform_grid(0, sensor_pos.nrows(), 1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Sensor::Sinefit);
  rq.Mode("");
  rq.Target().perturbation = 0;

  // Each sinefit coeff. pair is treated as a retrieval quantity
  //
  for (Index i = 0; i < np; i++) {
    std::ostringstream sstr;
    sstr << "Period " << i;
    rq.Subtag(sstr.str());

    // "Grid" has length 2, set to period length
    grids[0] = Vector(2, period_lengths[i]);
    rq.Grids(grids);

    // Add it to the *jacobian_quantities*
    jq.push_back(rq);

    // Add pointing method to the jacobian agenda
    jacobian_agenda.add(Method{"my_period_index", i});
    jacobian_agenda.add(Method{"jacobianCalcSinefit", {}, {{"period_index", "my_period_index"}}});
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
                         const Index& period_index) {
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  std::ostringstream sstr;
  sstr << "Period " << period_index;
  for (iq = 0; iq < jacobian_quantities.nelem() && !found; iq++) {
    if (jacobian_quantities[iq] == Jacobian::Sensor::Sinefit &&
        jacobian_quantities[iq].Subtag() == sstr.str()) {
      found = true;
      break;
    }
  }
  if (!found) {
    throw std::runtime_error(
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
    Numeric a = (sensor_response_f_grid[f] - sensor_response_f_grid[0]) *
                Constant::two_pi / jg[0][0];
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
  const Index row4 = rowind.offset;
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
void jacobianAddSurfaceQuantity(const Workspace&,
                                ArrayOfRetrievalQuantity& jq,
                                Agenda& jacobian_agenda,
                                const Vector& rq_lat_grid,
                                const Vector& rq_lon_grid,
                                const String& quantity) {
  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Special::SurfaceString && jq[it].Subtag() == quantity) {
      std::ostringstream os;
      os << quantity << " is already included as a surface variable "
         << "in *jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target() = Jacobian::Target(Jacobian::Special::SurfaceString, quantity);
  rq.Subtag(quantity);
  rq.Grids({rq_lat_grid, rq_lon_grid});

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  // Add dummy
  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Temperatures (atmospheric):
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddTemperature(const Workspace&,
                            ArrayOfRetrievalQuantity& jq,
                            Agenda& jacobian_agenda,
                            const Vector& rq_p_grid,
                            const Vector& rq_lat_grid,
                            const Vector& rq_lon_grid,
                            const String& hse) {
  // Check that temperature is not already included in the jacobian.
  // We only check the main tag.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Atm::Temperature) {
      std::ostringstream os;
      os << "Temperature is already included in *jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Set subtag
  String subtag;
  if (hse == "on") {
    subtag = "HSE on";
  } else if (hse == "off") {
    subtag = "HSE off";
  } else {
    std::ostringstream os;
    os << "The keyword for hydrostatic equilibrium can only be set to\n"
       << "\"on\" or \"off\"\n";
    throw std::runtime_error(os.str());
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Subtag(subtag);
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});
  rq.Target(Jacobian::Target(Jacobian::Atm::Temperature));
  rq.Target().perturbation = 0.1;

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Winds:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddWind(const Workspace&,
                     ArrayOfRetrievalQuantity& jq,
                     Agenda& jacobian_agenda,
                     const Vector& rq_p_grid,
                     const Vector& rq_lat_grid,
                     const Vector& rq_lon_grid,
                     const String& component,
                     const Numeric& dfrequency) {  
  // Create the new retrieval quantity
  const auto opt = Options::toWindMagJacobianOrThrow(component);
  RetrievalQuantity rq;
  switch (opt) {
    case Options::WindMagJacobian::u:
      rq.Target(Jacobian::Target(Jacobian::Atm::WindU));
      break;
    case Options::WindMagJacobian::v:
      rq.Target(Jacobian::Target(Jacobian::Atm::WindV));
      break;
    case Options::WindMagJacobian::w:
      rq.Target(Jacobian::Target(Jacobian::Atm::WindW));
      break;
    case Options::WindMagJacobian::strength:
      rq.Target(Jacobian::Target(Jacobian::Atm::WindMagnitude));
      break;
    case Options::WindMagJacobian::FINAL:
      ARTS_ASSERT(false, "This error should be caught earlier")
  }

  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    ARTS_USER_ERROR_IF (jq[it].Target().sameTargetType(rq.Target()),
      "The wind component:\n", component, "\n"
      "is already included in *jacobian_quantities*.")
  }

  rq.Subtag(component);  // nb.  This should be possible to remove...
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});
  rq.Target().perturbation = dfrequency;

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Magnetic field:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddMagField(const Workspace&,
                         ArrayOfRetrievalQuantity& jq,
                         Agenda& jacobian_agenda,
                         const Vector& rq_p_grid,
                         const Vector& rq_lat_grid,
                         const Vector& rq_lon_grid,
                         const String& component,
                         const Numeric& dB) {  
  // Create the new retrieval quantity
  const auto opt = Options::toWindMagJacobianOrThrow(component);
  RetrievalQuantity rq;
  switch (opt) {
    case Options::WindMagJacobian::u:
      rq.Target(Jacobian::Target(Jacobian::Atm::MagneticU));
      break;
    case Options::WindMagJacobian::v:
      rq.Target(Jacobian::Target(Jacobian::Atm::MagneticV));
      break;
    case Options::WindMagJacobian::w:
      rq.Target(Jacobian::Target(Jacobian::Atm::MagneticW));
      break;
    case Options::WindMagJacobian::strength:
      rq.Target(Jacobian::Target(Jacobian::Atm::MagneticMagnitude));
      break;
    case Options::WindMagJacobian::FINAL:
      ARTS_ASSERT(false, "This error should be caught earlier")
  }
  
  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    ARTS_USER_ERROR_IF (jq[it].Target().sameTargetType(rq.Target()),
                        "The magnetic component:\n", component, "\n"
                        "is already included in *jacobian_quantities*.")
  }

  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});
  rq.Target().perturbation = dB;

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  // Add gas species method to the jacobian agenda
  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Catalog parameters:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddShapeCatalogParameter(const Workspace&,
                                      ArrayOfRetrievalQuantity& jq,
                                      Agenda& jacobian_agenda,
                                      const QuantumIdentifier& line_identity,
                                      const String& species,
                                      const String& variable,
                                      const String& coefficient) {
  const auto jpt = select_derivativeLineShape(variable, coefficient);

  // Create the quantity
  RetrievalQuantity rq;
  rq.Mode(species);
  rq.Grids(ArrayOfVector(0, Vector(0)));
  
  // Map the species
  if (species == LineShape::self_broadening) {
    rq.Target(Jacobian::Target(jpt, line_identity, line_identity.Species()));
  } else if (species == LineShape::bath_broadening) {
    rq.Target(Jacobian::Target(jpt, line_identity, Species::Species::Bath));
  } else {
    rq.Target(Jacobian::Target(jpt, line_identity, SpeciesTag(species).Spec()));
  }

  // Test this is not a copy
  for (auto& q : jq)
    if (q.HasSameInternalsAs(rq))
      throw std::runtime_error("Error with copies of the quantities");

  // Append and do housekeeping
  jq.push_back(rq);
  jacobian_agenda.add({"jacobianCalcDoNothing"});  // old code activation
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddShapeCatalogParameters(
    const Workspace& ws,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    const ArrayOfQuantumIdentifier& line_identities,
    const ArrayOfString& species,
    const ArrayOfString& variables,
    const ArrayOfString& coefficients) {
  if (not(line_identities.nelem() or species.nelem() or variables.nelem() or
          coefficients.nelem()))
    throw std::runtime_error("Must have at least 1-long lists for all GINs");

  ArrayOfString vars;
  if (variables[0] == "ALL")
    vars = ArrayOfString(LineShape::enumstrs::VariableNames);
  else
    vars = variables;

  ArrayOfString coeffs;
  if (coefficients[0] == "ALL")
    coeffs = ArrayOfString(Options::enumstrs::LineShapeCoeffNames);
  else
    coeffs = coefficients;

  for (auto& l : line_identities)
    for (auto& s : species)
      for (auto& v : vars)
        for (auto& c : coeffs)
          jacobianAddShapeCatalogParameter(
              ws, jq, jacobian_agenda, l, s, v, c);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddBasicCatalogParameter(const Workspace&,
                                      ArrayOfRetrievalQuantity& jq,
                                      Agenda& jacobian_agenda,
                                      const QuantumIdentifier& catalog_identity,
                                      const String& catalog_parameter) {
  // Create the new retrieval quantity
  const auto opt = Options::toBasicCatParamJacobianOrThrow(catalog_parameter);
  RetrievalQuantity rq;
  switch (opt) {
    case Options::BasicCatParamJacobian::LineCenter:
      rq.Target(Jacobian::Target(Jacobian::Line::Center, catalog_identity, catalog_identity.Species()));
      break;
    case Options::BasicCatParamJacobian::LineStrength:
      rq.Target(Jacobian::Target(Jacobian::Line::Strength, catalog_identity, catalog_identity.Species()));
      break;
    case Options::BasicCatParamJacobian::FINAL:
      ARTS_ASSERT(false, "This error should be caught earlier")
  }
  
  // Check that this is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    ARTS_USER_ERROR_IF (rq.Target().sameTargetType(jq[it].Target()),
      "The catalog identifier:\n",
      catalog_identity, " for ID: ", catalog_identity, "\n"
      "is already included in jacobian_quantities*.")
  }

  rq.Grids(ArrayOfVector(0, Vector(0)));

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddBasicCatalogParameters(
    const Workspace& ws,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    const ArrayOfQuantumIdentifier& catalog_identities,
    const ArrayOfString& catalog_parameters) {

  for (auto& qi : catalog_identities)
    for (auto& param : catalog_parameters)
      jacobianAddBasicCatalogParameter(
          ws, jq, jacobian_agenda, qi, param);
}

//----------------------------------------------------------------------------
// NLTE temperature:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddNLTE(const Workspace&,
                     ArrayOfRetrievalQuantity& jq,
                     Agenda& jacobian_agenda,
                     const Vector& rq_p_grid,
                     const Vector& rq_lat_grid,
                     const Vector& rq_lon_grid,
                     const QuantumIdentifier& energy_level_identity,
                     const Numeric& dx) {
  // Check that this species is not already included in the jacobian.
  for (Index it = 0; it < jq.nelem(); it++) {
    if (jq[it] == Jacobian::Line::NLTE and
        jq[it].QuantumIdentity() == energy_level_identity) {
      std::ostringstream os;
      os << "The NLTE identifier:\n"
         << energy_level_identity << "\nis already included in "
         << "*jacobian_quantities*.";
      throw std::runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Target(Jacobian::Target(Jacobian::Line::NLTE, energy_level_identity, energy_level_identity.Species()));
  rq.Target().perturbation = dx;
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});
  
  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

void jacobianAddNLTEs(const Workspace& ws,
                      ArrayOfRetrievalQuantity& jq,
                      Agenda& jacobian_agenda,
                      const Vector& rq_p_grid,
                      const Vector& rq_lat_grid,
                      const Vector& rq_lon_grid,
                      const ArrayOfQuantumIdentifier& energy_level_identities,
                      const Numeric& dx) {
  for (const auto& qi : energy_level_identities)
    jacobianAddNLTE(ws,
                    jq,
                    jacobian_agenda,
                    rq_p_grid,
                    rq_lat_grid,
                    rq_lon_grid,
                    qi,
                    dx);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSpecialSpecies(const Workspace&,
                               ArrayOfRetrievalQuantity& jq,
                               Agenda& jacobian_agenda,
                               const Vector& rq_p_grid,
                               const Vector& rq_lat_grid,
                               const Vector& rq_lon_grid,
                               const String& species) {
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Grids({rq_p_grid, rq_lat_grid, rq_lon_grid});

  // Make sure modes are valid and complain if they are repeated
  if (species == "electrons") {
    for (Index it = 0; it < jq.nelem(); it++) {
      if (jq[it] == Jacobian::Atm::Electrons) {
        std::ostringstream os;
        os << "Electrons are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.Target(Jacobian::Target(Jacobian::Atm::Electrons));

  } else if (species == "particulates") {
    for (Index it = 0; it < jq.nelem(); it++) {
      if (jq[it] == Jacobian::Atm::Particulates) {
        std::ostringstream os;
        os << "Particulates are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.Target(Jacobian::Target(Jacobian::Atm::Particulates));
    
  } else {
    std::ostringstream os;
    os << "Unknown special species jacobian: \"" << species
       << "\"\nPlease see *jacobianAddSpecialSpecies* for viable options.";
    throw std::runtime_error(os.str());
  }

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  jacobian_agenda.add({"jacobianCalcDoNothing"});
}

//----------------------------------------------------------------------------
// Adjustments and transformations
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAdjustAndTransform(
    Matrix& jacobian,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& x) {
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
    if (jacobian_quantities[q].Target().isSpeciesVMR() &&
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
                                     const Vector& offset_vector) {
  if (jqs.empty()) {
    std::runtime_error(
        "Jacobian quantities is empty, so there is nothing to add the "
        "transformation to.");
  }

  Index nelem = jqs.back().Grids().nelem();

  if (!(nelem == transformation_matrix.nrows())) {
    std::runtime_error(
        "Dimension of transformation matrix incompatible with retrieval grids.");
  }
  if (!(nelem == offset_vector.nelem())) {
    std::runtime_error(
        "Dimension of offset vector incompatible with retrieval grids.");
  }

  jqs.back().SetTransformationMatrix(Matrix{transpose(transformation_matrix)});
  jqs.back().SetOffsetVector(offset_vector);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianSetFuncTransformation(ArrayOfRetrievalQuantity& jqs,
                                   const String& transformation_func,
                                   const Numeric& z_min,
                                   const Numeric& z_max) {
  if (jqs.empty())
    throw std::runtime_error(
        "Jacobian quantities is empty, so there is nothing to add the "
        "transformation to.");

  if (transformation_func == "none") {
    jqs.back().SetTransformationFunc("");
    return;
  }

  Vector pars;

  if (transformation_func == "atanh") {
    if (z_max <= z_min)
      throw std::runtime_error(
          "For option atanh, the GIN *z_max* must be set and be > z_min.");
    pars.resize(2);
    pars[0] = z_min;
    pars[1] = z_max;
  } else if (transformation_func == "log" || transformation_func == "log10") {
    pars.resize(1);
    pars[0] = z_min;
  } else {
    std::ostringstream os;
    os << "Valid options for *transformation_func* are:\n"
       << "\"none\", \"log\", \"log10\" and \"atanh\"\n"
       << "But found: \"" << transformation_func << "\"";
    throw std::runtime_error(os.str());
  }

  jqs.back().SetTransformationFunc(transformation_func);
  jqs.back().SetTFuncParameters(pars);
}

//----------------------------------------------------------------------------
// Methods for doing perturbations
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPerturb(Tensor3& perturbed_field,
                    const Tensor3& original_field,
                    const Vector& p_ret_grid,
                    const Vector& lat_ret_grid,
                    const Vector& lon_ret_grid,
                    const Index& pert_index,
                    const Numeric& pert_size,
                    const String& pert_mode) {

  // Pack retrieval grids into an ArrayOfVector
  ArrayOfVector ret_grids(3);
  ret_grids[0] = p_ret_grid;
    ret_grids[1] = lat_ret_grid;
      ret_grids[2] = lon_ret_grid;

  // Find mapping from retrieval grids to atmospheric grids
  ArrayOfGridPos gp_p, gp_lat, gp_lon;
  Index n_p, n_lat, n_lon;
  ARTS_ASSERT(false)
    /*get_gp_rq_to_atmgrids(gp_p,
                        gp_lat,
                        gp_lon,
                        n_p,
                        n_lat,
                        n_lon,
                        ret_grids,
                        3,
                        p_grid,
                        lat_grid,
                        lon_grid);*/

  // Now we can chec *pert_index*
  if (pert_index<0){
    throw std::runtime_error("Bad *pert_index*. It is negative.");
  }
  const Index n_tot = n_p * n_lat * n_lon;
  if (pert_index >= n_tot){
    throw std::runtime_error("Bad *pert_index*. It is too high with respect "
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
    throw std::runtime_error("Bad *pert_mode*. Allowed choices are: "
                        """absolute"" and ""relative"".");
  }
  
  // Map x to a perturbation defined at atmospheric grids
  Tensor3 x3d(n_p, n_lat, n_lon), pert(n_p, n_lat, n_lon);
  reshape(x3d, x);
  regrid_atmfield_by_gp_oem(pert, x3d, gp_p, gp_lat, gp_lon);
  
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
                             const Vector& p_grid,
                             const Vector& lat_grid,
                             const Vector& lon_grid,
                             const Tensor3& original_field,
                             const Index& pert_index,
                             const Numeric& pert_size,
                             const String& pert_mode) {
  // Some sizes
  const Index n_p = p_grid.nelem();
  const Index n_lat = lat_grid.nelem();
  const Index n_lon = lon_grid.nelem();
  
  // Check input
  chk_atm_field("original_field",
                original_field,
                3,
                p_grid,
                lat_grid,
                lon_grid,
                false );
  if (pert_index<0){
    throw std::runtime_error("Bad *pert_index*. It is negative.");
  }
  if (pert_index >= n_p * n_lat * n_lon){
    throw std::runtime_error("Bad *pert_index*. It is too high with respect "
                        "to length of atmospheric grids.");
  }    

  // Determine indexes with respect to atmospheric grids
  Index tot_index = pert_index;
  const Index lon_index = tot_index / (n_lat * n_p);
  tot_index -= lon_index * n_lat * n_p;
  const Index lat_index = tot_index / n_p;
  tot_index -= lat_index * n_p;
  const Index p_index = tot_index;
  
  // Init perturbed_field, if not equal to original_field
  if (&perturbed_field != &original_field) {
    perturbed_field = original_field;
  }

  // Perturb
  if (pert_mode == "absolute") {
    perturbed_field(p_index, lat_index, lon_index) += pert_size;
  } else if (pert_mode == "relative") {
    perturbed_field(p_index, lat_index, lon_index) *= 1 + pert_size;
  } else {
    ARTS_USER_ERROR(R"(Bad *pert_mode*. Allowed choices are:  "absolute" and "relative".)");
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexNumberOfAtmosphericPoints(Index& n,
                                    const Vector& p_grid,
                                    const Vector& lat_grid,
                                    const Vector& lon_grid) {
  const Index n_p = p_grid.nelem();
  const Index n_lat = lat_grid.nelem();
  const Index n_lon = lon_grid.nelem();

  n = n_p * n_lat * n_lon;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianFromTwoY(Matrix& jacobian,
                    const Vector& y_pert,
                    const Vector& y,
                    const Numeric& pert_size) {
  const Index n = y.nelem();
  if( y_pert.nelem() != n ){
    throw std::runtime_error("Inconsistency in length of *y_pert* and *y*.");
  }
  jacobian = ExhaustiveConstMatrixView{y_pert};
  jacobian -= ExhaustiveConstMatrixView{y};
  jacobian /= pert_size;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianFromYbatch(Matrix& jacobian,
                    const ArrayOfVector& ybatch,
                    const Vector& y,
                    const Numeric& pert_size) {
  const Index n = y.nelem();
  const Index l = ybatch.nelem();
  if (l>0){
      if( ybatch[0].nelem() != n )
        throw std::runtime_error("Inconsistency in length of y and ybatch[0].");
    }
  jacobian.resize(n,l);
  for (Index i=0; i<l; i++) {
    jacobian(joker,i) = ybatch[i];
    jacobian(joker,i) -= y;
  }
  jacobian /= pert_size;
}
