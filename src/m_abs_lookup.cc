/*!
  \file   m_abs_lookup.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "absorption.h"
#include <workspace.h>
#include "arts_omp.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "gas_abs_lookup.h"
#include "gridded_fields.h"
#include "interp.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "matpack_view.h"
#include "physics_funcs.h"
#include "rng.h"
#include <rtepack.h>
#include "species_tags.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupInit(GasAbsLookup& x) {
  x = GasAbsLookup();
}

//! Find continuum species in abs_species.
/*! 
  Returns an index array with indexes of those species in abs_species
  that have continuum tags that require h2o_abs, and hence require
  nonlinear treatment in the absorption lookup table.

  H2O itself is ignored here since that is treated separately.

  We are a bit conservative here, it is possible that some of the
  continua do not really require H2O. Check yourself, if you want, and
  improve the guessing here.

  \retval   cont         indices of those species with continua.
  \param    abs_species  list of absorption species.
  
  \author Stefan Buehler
  \date   2007-11-16
*/
void find_nonlinear_continua(ArrayOfIndex& cont,
                             const ArrayOfArrayOfSpeciesTag& abs_species) {
  cont.resize(0);

  // This is quite complicated, unfortunately. The approach here
  // is copied from abs_xsec_per_speciesAddConts. For explanation,
  // see there.

  // Loop tag groups:
  for (Index i = 0; i < abs_species.size(); ++i) {
    // Loop tags in tag group
    for (Index s = 0; s < abs_species[i].size(); ++s) {
      // Check for continuum tags
      if (abs_species[i][s].type == Species::TagType::Predefined ||
          abs_species[i][s].type == Species::TagType::Cia) {
        const String thisname = abs_species[i][s].Name();
        // Ok, now we know this is a continuum tag.

        // Check whether we want nonlinear treatment for
        // this or not. We have three classes of continua:
        // 1. Those that we know do not require it
        // 2. Those that we know require h2o_abs
        // 3. Those for which we do not know

        // The list here is not at all perfect, just a first
        // guess. If you need particular models, you better
        // check that the right thing is done with your model.

        // 1. Continua known to not use h2o_abs
        // We take also H2O itself here, since this is
        // handled separately
        if (Species::fromShortName("H2O") == abs_species[i][s].Spec() ||
            "N2-" == thisname.substr(0, 3) || "CO2-" == thisname.substr(0, 4) ||
            "O2-CIA" == thisname.substr(0, 6) ||
            "O2-v0v" == thisname.substr(0, 6) ||
            "O2-v1v" == thisname.substr(0, 6) ||
            "H2-CIA" == thisname.substr(0, 6) ||
            "He-CIA" == thisname.substr(0, 6) ||
            "CH4-CIA" == thisname.substr(0, 7) ||
            "liquidcloud-MPM93" == thisname.substr(0, 17) ||
            "liquidcloud-ELL07" == thisname.substr(0, 17)) {
          break;
        }

        // 2. Continua known to use h2o_abs
        if ("O2-" == thisname.substr(0, 3)) {
          cont.push_back(i);
          break;
        }

        // 3. abs_species tags that are NOT allowed in LUT
        // calculations
        if ("icecloud-" == thisname.substr(0, 9) ||
            "rain-" == thisname.substr(0, 5)) {
          std::ostringstream os;
          os << "Tag " << thisname << " not allowed in absorption "
             << "lookup tables.";
          throw std::runtime_error(os.str());
        }

        // If we get here, then the tag was neither in the
        // posivitive nor in the negative list. We through a
        // runtime error.
        std::ostringstream os;
        os << "Unknown whether tag " << thisname
           << " is a nonlinear species (i.e. uses h2o_abs) or not.\n"
           << "Cannot set abs_nls automatically.";
        throw std::runtime_error(os.str());
      }
    }
  }
}

//! Choose species for abs_nls
/*!
  Make an intelligent choice for abs_nls, based on abs_species.

  \author Stefan Buehler

  \param[out] abs_nls     The list of nonlinear species.
  \param[in]  abs_species Absorption species.
*/
void choose_abs_nls(ArrayOfArrayOfSpeciesTag& abs_nls,
                    const ArrayOfArrayOfSpeciesTag& abs_species) {
  abs_nls.resize(0);

  // Add all H2O species as non-linear:
  Index next_h2o = 0;
  while (-1 !=
         (next_h2o = find_next_species(
              abs_species, Species::fromShortName("H2O"), next_h2o))) {
    abs_nls.push_back(abs_species[next_h2o]);
    ++next_h2o;
  }

  // Certain continuum models also depend on abs_h2o. There is a
  // helper function that contains a list of these.
  ArrayOfIndex cont;
  find_nonlinear_continua(cont, abs_species);

  // Add these to abs_nls:
  for (Index i = 0; i < cont.size(); ++i) {
    abs_nls.push_back(abs_species[cont[i]]);
  }
}

//! Chose the temperature perturbations abs_t_pert
/*!  
  This simple function creates a vector of temperature
  perturbations, relative to the reference temperature profile, that
  covers the minimum and maximum temperature profile. 
  
  \author Stefan Buehler

  \param[out] abs_t_pert Temperature perturbations
  \param[in] abs_t       Reference temperature profile
  \param[in] tmin        Minimum temperature profile
  \param[in] tmax        Maximum temperature profile
  \param[in] t_step      Temperature perturbation step
*/
void choose_abs_t_pert(Vector& abs_t_pert,
                       ConstVectorView abs_t,
                       ConstVectorView tmin,
                       ConstVectorView tmax,
                       const Numeric& step,
                       const Index& p_interp_order,
                       const Index& t_interp_order) {
  // The code to find out the range for perturbation is a bit
  // complicated. The problem is that, since we use higher order
  // interpolation for p, we may require temperatures well outside the
  // local min/max range at the reference level. We solve this by
  // really looking at the min/max range for those levels required by
  // the p_interp_order.

  Numeric mindev = 1e9;
  Numeric maxdev = -1e9;

  Vector the_grid=uniform_grid(0, abs_t.nelem(), 1);
  for (Index i = 0; i < the_grid.nelem(); ++i) {
    const Index idx0 = my_interp::pos_finder<true>(i, Numeric(i), the_grid, p_interp_order);

    for (Index j = 0; j < p_interp_order+1; ++j) {
      // Our pressure grid for the lookup table may be coarser than
      // the original one for the batch cases. This may lead to max/min
      // values in the original data that exceed those we assumed
      // above. We add some extra margin here to account for
      // that. (The margin is +-10K)

      Numeric delta_min = tmin[i] - abs_t[idx0 + j] - 10;
      Numeric delta_max = tmax[i] - abs_t[idx0 + j] + 10;

      if (delta_min < mindev) mindev = delta_min;
      if (delta_max > maxdev) maxdev = delta_max;
    }
  }

  // We divide the interval between mindev and maxdev, so that the
  // steps are of size *step* or smaller. But we also need at least
  // *t_interp_order*+1 points.
  Index div = t_interp_order;
  Numeric effective_step;
  do {
    effective_step = (maxdev - mindev) / (Numeric)div;
    ++div;
  } while (effective_step > step);

  abs_t_pert =uniform_grid(mindev, div, effective_step);
}

//! Chose the H2O perturbations abs_nls_pert
/*!  
  This simple function creates a vector of fractional H2O VMR
  perturbations, relative to the reference H2O profile, that
  covers the minimum and maximum profile. 
  
  \author Stefan Buehler

  \param[out] abs_nls_pert H2O VMR perturbations
  \param[in] refprof       Reference profile
  \param[in] minprof       Minimum profile
  \param[in] maxprof       Maximum profile
  \param[in] step          Fractional perturbation step
*/
void choose_abs_nls_pert(Vector& abs_nls_pert,
                         ConstVectorView refprof,
                         ConstVectorView minprof,
                         ConstVectorView maxprof,
                         const Numeric& step,
                         const Index& p_interp_order,
                         const Index& nls_interp_order) {
  // The code to find out the range for perturbation is a bit
  // complicated. The problem is that, since we use higher order
  // interpolation for p, we may require humidities well outside the
  // local min/max range at the reference level. We solve this by
  // really looking at the min/max range for those levels required by
  // the p_interp_order.

  Numeric mindev = 0;
  Numeric maxdev = -1e9;

  // mindev is set to zero from the start, since we always want to
  // include 0.

  Vector the_grid=uniform_grid(0, refprof.nelem(), 1);
  for (Index i = 0; i < the_grid.nelem(); ++i) {
    //       cout << "!!!!!! i = " << i << "\n";
    //       cout << " min/ref/max = " << minprof[i] << " / "
    //            << refprof[i] << " / "
    //            << maxprof[i] << "\n";
    
    const Index idx0 = my_interp::pos_finder<true>(i, Numeric(i), the_grid, p_interp_order);

    for (Index j = 0; j < p_interp_order+1; ++j) {
      //           cout << "!!!!!! j = " << j << "\n";
      //           cout << "  ref[j] = " << refprof[gp.idx[j]] << "   ";
      //           cout << "  minfrac[j] = " << minprof[i] / refprof[gp.idx[j]] << "   ";
      //           cout << "  maxfrac[j] = " << maxprof[i] / refprof[gp.idx[j]] << "  \n";

      // Our pressure grid for the lookup table may be coarser than
      // the original one for the batch cases. This may lead to max/min
      // values in the original data that exceed those we assumed
      // above. We add some extra margin to the max value here to account for
      // that. (The margin is a factor of 2.)

      Numeric delta_min = minprof[i] / refprof[idx0 + j];
      Numeric delta_max = 2 * maxprof[i] / refprof[idx0 + j];

      if (delta_min < mindev) mindev = delta_min;
      // do not update maxdev, when delta_max is infinity (this results from
      // refprof being 0)
      if (!std::isinf(delta_max) && (delta_max > maxdev)) maxdev = delta_max;
    }
  }

  bool allownegative = false;
  if (mindev < 0) {
    allownegative = true;
  }

  if (!allownegative) {
    mindev = 0;
  }

  if (std::isinf(maxdev)) {
    std::ostringstream os;
    os << "Perturbation upper limit is infinity (likely due to the reference\n"
       << "profile being 0 at at least one pressure level). Can not work\n"
       << "with that.";
    throw std::runtime_error(os.str());
  }

  // We divide the interval between mindev and maxdev, so that the
  // steps are of size *step* or smaller. But we also need at least
  // *nls_interp_order*+1 points.
  Index div = nls_interp_order;
  Numeric effective_step;
  do {
    effective_step = (maxdev - mindev) / (Numeric)div;
    ++div;
  } while (effective_step > step);

  abs_nls_pert = uniform_grid(mindev, div, effective_step);

  // If there are negative values, we also add 0. The reason for this
  // is that 0 is a turning point.
  if (allownegative) {
    VectorInsertGridPoints(abs_nls_pert, abs_nls_pert, {0});
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const ArrayOfString& names) {
  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for (Index i = 0; i < names.size(); ++i) {
    abs_species.emplace_back(names[i]);
  }

  check_abs_species(abs_species);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd2(  // WS Output:
    const Workspace& ws,
    ArrayOfArrayOfSpeciesTag& abs_species,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    Index& propmat_clearsky_agenda_checked,
    // WS Generic Input:
    const Vector& rq_p_grid,
    const Vector& rq_lat_grid,
    const Vector& rq_lon_grid,
    // Control Parameters:
    const String& species,
    const String& mode) {
  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  // Add species to *abs_species*
  abs_species.emplace_back(species);

  check_abs_species(abs_species);

  // Do retrieval part
  jacobianAddAbsSpecies(ws,
                        jq,
                        jacobian_agenda,
                        rq_p_grid,
                        rq_lat_grid,
                        rq_lon_grid,
                        species,
                        mode,
                        1);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesInit(ArrayOfArrayOfSpeciesTag& abs_species) {
  abs_species.resize(0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesSet(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const ArrayOfString& names) {
  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  abs_species.resize(names.size());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for (Index i = 0; i < names.size(); ++i) {
    // This part has now been moved to array_species_tag_from_string.
    // Call this function.
    abs_species[i] = ArrayOfSpeciesTag(names[i]);
  }

  check_abs_species(abs_species);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupAdapt(GasAbsLookup& abs_lookup,
                     Index& abs_lookup_is_adapted,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const Vector& f_grid) {
  abs_lookup.Adapt(abs_species, f_grid);
  abs_lookup_is_adapted = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFromLookup(
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    const GasAbsLookup& abs_lookup,
    const Index& abs_lookup_is_adapted,
    const Index& abs_p_interp_order,
    const Index& abs_t_interp_order,
    const Index& abs_nls_interp_order,
    const Index& abs_f_interp_order,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const Numeric& extpolfac,
    const Index& no_negatives) {
  // Variables needed by abs_lookup.Extract:
  Matrix abs_scalar_gas, dabs_scalar_gas_df, dabs_scalar_gas_dt;

  // Check if the table has been adapted:
  if (1 != abs_lookup_is_adapted)
    throw std::runtime_error(
        "Gas absorption lookup table must be adapted,\n"
        "use method abs_lookupAdapt.");

  const bool do_jac = supports_lookup(jacobian_quantities);
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);

  const Vector a_vmr_list = [&]() {
    Vector vmr(abs_species.size());
    std::transform(abs_species.begin(), abs_species.end(), vmr.begin(),
                   [&](const ArrayOfSpeciesTag& spec) -> Numeric { return atm_point[spec]; });
    return vmr;
  }();

  // The combination of doing frequency jacobian together with an
  // absorption lookup table is quite dangerous. If the frequency
  // interpolation order for the table is zero, the Jacobian will be
  // zero, and the cause for this may be difficult for a user to
  // find. So we do not allow this combination.
  if (do_freq_jac and (1 > abs_f_interp_order))
    throw std::runtime_error("Wind/frequency Jacobian is not possible without at least first\n"
			     "order frequency interpolation in the lookup table.  Please use\n"
			     "abs_f_interp_order>0 or remove wind/frequency Jacobian.");
  
  // The function we are going to call here is one of the few helper
  // functions that adjust the size of their output argument
  // automatically.
  abs_lookup.Extract(abs_scalar_gas,
                     select_abs_species,
                     abs_p_interp_order,
                     abs_t_interp_order,
                     abs_nls_interp_order,
                     abs_f_interp_order,
                     atm_point.pressure,
                     atm_point.temperature,
                     a_vmr_list,
                     f_grid,
                     extpolfac);
  if (do_freq_jac) {
    Vector dfreq = f_grid;
    dfreq += df;
    abs_lookup.Extract(dabs_scalar_gas_df,
                       select_abs_species,
                       abs_p_interp_order,
                       abs_t_interp_order,
                       abs_nls_interp_order,
                       abs_f_interp_order,
                       atm_point.pressure,
                       atm_point.temperature,
                       a_vmr_list,
                       dfreq,
                       extpolfac);
  }
  if (do_temp_jac) {
    const Numeric dtemp = atm_point.temperature + dt;
    abs_lookup.Extract(dabs_scalar_gas_dt,
                       select_abs_species,
                       abs_p_interp_order,
                       abs_t_interp_order,
                       abs_nls_interp_order,
                       abs_f_interp_order,
                       atm_point.pressure,
                       dtemp,
                       a_vmr_list,
                       f_grid,
                       extpolfac);
  }

  if (no_negatives){
    //Check for negative values due to interpolation and set them to zero
    for (Index ir = 0; ir < abs_scalar_gas.nrows(); ir++){
      for (Index ic = 0; ic < abs_scalar_gas.ncols(); ic++){
        if (abs_scalar_gas(ir, ic)<0) abs_scalar_gas(ir, ic)=0;
      }
    }
  }

  // Now add to the right place in the absorption matrix.

  if (not do_jac) {
    for (Index ii = 0; ii < abs_scalar_gas.nrows(); ii++) {
      for (Index iv = 0; iv < abs_scalar_gas.ncols(); iv++)
        propmat_clearsky[iv].A() += abs_scalar_gas(ii, iv);
    }
  } else {
    for (Index isp = 0; isp < abs_scalar_gas.nrows(); isp++) {
      for (Index iv = 0; iv < abs_scalar_gas.ncols(); iv++) {
        propmat_clearsky[iv].A() += abs_scalar_gas(isp, iv);
        for (Index iq = 0; iq < jacobian_quantities.size(); iq++) {
          const auto& deriv = jacobian_quantities[iq];
          
          if (not deriv.propmattype()) continue;
          
          if (deriv == Jacobian::Atm::Temperature) {
            dpropmat_clearsky_dx(iq, iv).A() +=
                (dabs_scalar_gas_dt(isp, iv) - abs_scalar_gas(isp, iv)) / dt;
          } else if (is_frequency_parameter(deriv)) {
            dpropmat_clearsky_dx(iq, iv).A() +=
                (dabs_scalar_gas_df(isp, iv) - abs_scalar_gas(isp, iv)) / df;
          } else if (deriv == abs_species[isp]) {
            // WARNING:  If CIA in list, this scales wrong by factor 2
              dpropmat_clearsky_dx(iq, iv).A() += (std::isnormal(a_vmr_list[isp])) ? abs_scalar_gas(isp, iv) / a_vmr_list[isp] : 0;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromGasAbsLookup(Vector& f_grid,
                            const GasAbsLookup& abs_lookup) {
  const Vector& lookup_f_grid = abs_lookup.GetFgrid();
  f_grid.resize(lookup_f_grid.nelem());
  f_grid = lookup_f_grid;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromGasAbsLookup(Vector& p_grid,
                            const GasAbsLookup& abs_lookup) {
  const Vector& lookup_p_grid = abs_lookup.GetPgrid();
  p_grid.resize(lookup_p_grid.nelem());
  p_grid = lookup_p_grid;
}
