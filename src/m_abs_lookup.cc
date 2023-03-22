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
#include "agenda_class.h"
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "gas_abs_lookup.h"
#include "global_data.h"
#include "gridded_fields.h"
#include "interp.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "messages.h"
#include "physics_funcs.h"
#include "propagationmatrix.h"
#include "rng.h"
#include "species_tags.h"

using GriddedFieldGrids::GFIELD4_FIELD_NAMES;
using GriddedFieldGrids::GFIELD4_P_GRID;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupInit(GasAbsLookup& x, const Verbosity& verbosity) {
  ArtsOut2 out2(verbosity);
  // Nothing to do here.
  // That means, we rely on the default constructor.

  x = GasAbsLookup();
  out2 << "  Created an empty gas absorption lookup table.\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupCalc(  // Workspace reference:
    Workspace& ws,
    // WS Output:
    GasAbsLookup& abs_lookup,
    Index& abs_lookup_is_adapted,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfSpeciesTag& abs_nls,
    const Vector& f_grid,
    const Vector& abs_p,
    const Matrix& abs_vmrs,
    const Vector& abs_t,
    const Vector& abs_t_pert,
    const Vector& abs_nls_pert,
    const Agenda& propmat_clearsky_agenda,
    // GIN
    const Numeric& lowest_vmr,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  ARTS_USER_ERROR_IF(
      propmat_clearsky_agenda.name() not_eq "propmat_clearsky_agenda",
      R"--(
Not propmat_clearsky_agenda:

)--",
      propmat_clearsky_agenda)

  ARTS_USER_ERROR_IF(lowest_vmr <= 0, R"--(
You must give a lowest_vmr value above 0.  This is because
the computations of the absorption lookup table uses true
absorption calculations but the output is in cross-sections

Your current lowest_vmr value is: )--", lowest_vmr)

  // We need this temporary variable to make a local copy of all VMRs,
  // where we then perturb the H2O profile as needed
  Matrix these_all_vmrs = abs_vmrs;

  // We will be calling an absorption agenda one species at a
  // time. This is better than doing all simultaneously, because is
  // saves memory and allows for consistent treatment of nonlinear
  // species. But it means we need local copies of species, line list,
  // and line shapes for agenda communication.

  // 2. Determine various important sizes:
  const Index n_species = abs_species.nelem();  // Number of abs species
  const Index n_nls = abs_nls.nelem();          // Number of nonlinear species
  const Index n_f_grid = f_grid.nelem();      // Number of frequency grid points
  const Index n_p_grid = abs_p.nelem();       // Number of presure grid points
  const Index n_t_pert = abs_t_pert.nelem();  // Number of temp. perturbations
  const Index n_nls_pert = abs_nls_pert.nelem();  // Number of VMR pert. for NLS

  // 3. Input to absorption calculations:

  // Absorption vmrs and temperature:
  Matrix this_vmr(1, n_p_grid);
  Vector abs_h2o(n_p_grid);
  Vector this_t;  // Has same dimension, but is
                  // initialized by assignment later.
  const EnergyLevelMap this_nlte_dummy;

  // Species list, lines, and line shapes, all with only 1 element:
  ArrayOfArrayOfSpeciesTag this_species(1);

  // List of active species for agenda call. Will always be filled with only
  // one species.
  ArrayOfIndex abs_species_active(1);

  // Local copy of nls_pert and t_pert:
  Vector these_nls_pert;  // Is resized every time when used
  Vector these_t_pert;    // Is resized later on

  // 4. Checks of input parameter correctness:
  const Index h2o_index = find_first_species(abs_species, Species::fromShortName("H2O"));

  if (h2o_index < 0) {
    // If there are nonlinear species, then at least one species must be
    // H2O. We will use that to set h2o_abs, and to perturb in the case
    // of nonlinear species.
    ARTS_USER_ERROR_IF(n_nls > 0,
                       "If you have nonlinear species, at least one\n"
                       "species must be a H2O species.")

    out2 << "  You have no H2O species. Absorption continua will not work.\n"
         << "  You should get a runtime error if you try them anyway.\n";
  }

  // abs_species, f_grid, and p_grid should not be empty:
  if (0 == n_species || 0 == n_f_grid || 0 == n_p_grid) {
    ostringstream os;
    os << "One of the required input variables is empty:\n"
       << "abs_species.nelem() = " << n_species << ",\n"
       << "f_grid.nelem() = " << n_f_grid << ",\n"
       << "abs_p.nelem() = " << n_p_grid << ".";
    throw runtime_error(os.str());
  }

  // Set up the index array abs_nls from the tag array
  // abs_nls. Give an error message if these
  // tags are not included in abs_species.
  ArrayOfIndex abs_nls_idx;
  for (Index i = 0; i < n_nls; ++i) {
    Index s;
    for (s = 0; s < n_species; ++s) {
      if (abs_nls[i] == abs_species[s]) {
        abs_nls_idx.push_back(s);
        break;
      }
    }
    if (s == n_species) {
      ostringstream os;
      os << "Did not find *abs_nls* tag group \""
         << abs_nls[i] << "\" in *abs_species*.";
      throw runtime_error(os.str());
    }
  }

  // Furthermore abs_nls_idx must not contain duplicate values:
  if (!is_unique(abs_nls_idx)) {
    ostringstream os;
    os << "The variable *abs_nls* must not have duplicate species.\n"
       << "Value of *abs_nls*: " << abs_nls_idx;
    throw runtime_error(os.str());
  }

  // VMR matrix must match species list and pressure levels:
  chk_size("abs_vmrs", abs_vmrs, n_species, n_p_grid);

  // Temperature vector must match number of pressure levels:
  chk_size("abs_t", abs_t, n_p_grid);

  // abs_nls_pert should only be not-empty if we have nonlinear species:
  if (((0 == n_nls) && (0 != n_nls_pert)) ||
      ((0 != n_nls) && (0 == n_nls_pert))) {
    ostringstream os;
    os << "You have to set both abs_nls and abs_nls_pert, or none.";
    throw runtime_error(os.str());
  }

  // 4.a Set up a logical array for the nonlinear species.
  ArrayOfIndex non_linear(n_species, 0);
  for (Index s = 0; s < n_nls; ++s) {
    non_linear[abs_nls_idx[s]] = 1;
  }

  // 5. Set general lookup table properties:
  abs_lookup.species = abs_species;  // Species list
  abs_lookup.nonlinear_species =
      abs_nls_idx;             // Nonlinear species   (e.g., H2O, O2)
  abs_lookup.f_grid = f_grid;  // Frequency grid
  abs_lookup.p_grid = abs_p;   // Pressure grid
  abs_lookup.vmrs_ref = abs_vmrs;
  abs_lookup.t_ref = abs_t;
  abs_lookup.t_pert = abs_t_pert;
  abs_lookup.nls_pert = abs_nls_pert;

  // 5.a. Set log_p_grid:
  abs_lookup.log_p_grid.resize(n_p_grid);
  transform(abs_lookup.log_p_grid, log, abs_lookup.p_grid);

  // 6. Create abs_lookup.xsec with the right dimensions:
  {
    Index a, b, c, d;

    if (0 == n_t_pert)
      a = 1;
    else
      a = n_t_pert;

    b = n_species + n_nls * (n_nls_pert - 1);

    c = n_f_grid;

    d = n_p_grid;

    abs_lookup.xsec.resize(a, b, c, d);
    abs_lookup.xsec = NAN;
  }

  // 6.a. Set up these_t_pert. This is done so that we can use the
  // same loop over the perturbations, independent of
  // whether we have temperature perturbations or not.
  if (0 != n_t_pert) {
    out2 << "  With temperature perturbations.\n";
    these_t_pert.resize(n_t_pert);
    these_t_pert = abs_t_pert;
  } else {
    out2 << "  No temperature perturbations.\n";
    these_t_pert.resize(1);
    these_t_pert = 0;
  }

  const Index these_t_pert_nelem = these_t_pert.nelem();

  // 7. Now we have to fill abs_lookup.xsec with the right values!

  String fail_msg;
  bool failed = false;

  // Loop species:
  for (Index i = 0, spec = 0; i < n_species; ++i) {
    // Skipping Zeeman and free_electrons species.
    // (Mixed tag groups between those and other species are not allowed.)
    if (abs_species[i].Zeeman() || abs_species[i].FreeElectrons() || abs_species[i].Particles()) {
      spec++;
      continue;
    }

    // spec is the index for the second dimension of abs_lookup.xsec.

    // Prepare absorption agenda input for this species:
    out2 << "  Doing species " << i + 1 << " of " << n_species << ": "
         << abs_species[i] << ".\n";

    // Set active species:
    abs_species_active[0] = i;

    // Get a dummy list of tag groups with only the current element:
    this_species[0].resize(abs_species[i].nelem());
    this_species[0] = abs_species[i];

    // Set up these_nls_pert. This is done so that we can use the
    // same loop over the perturbations, independent of
    // whether we have nonlinear species or not.
    if (non_linear[i]) {
      out2 << "  This is a species with H2O VMR perturbations.\n";
      these_nls_pert.resize(n_nls_pert);
      these_nls_pert = abs_nls_pert;
    } else {
      these_nls_pert.resize(1);
      these_nls_pert = 1;
    }

    // Loop these_nls_pert:
    for (Index s = 0; s < these_nls_pert.nelem(); ++s, ++spec) {
      // Remember, spec is the index for the second dimension of
      // abs_lookup.xsec

      if (non_linear[i]) {
        out2 << "  Doing H2O VMR variant " << s + 1 << " of " << n_nls_pert
             << ": " << abs_nls_pert[s] << ".\n";
      }

      // Make a local copy of the VMRs, and manipulate the H2O VMR within it.
      // Note: We do not need a runtime error check that h2o_index is ok here,
      // because earlier on we throw an error if there is no H2O species although we
      // need it. So, if h2o_indes is -1, we here simply assume that there
      // should not be a perturbation
      if (h2o_index >= 0) {
        these_all_vmrs(h2o_index, joker) = abs_vmrs(h2o_index, joker);
        these_all_vmrs(h2o_index, joker) *=
            these_nls_pert[s];  // Add perturbation
      }

      // VMR for this species (still needed by interfact to continua):
      // FIXME: This variable may go away eventually, when the continuum
      // part no longer needs it.
      this_vmr(0, joker) = these_all_vmrs(i, joker);

      // For abs_h2o, we can always add the perturbations (it will
      // not make a difference if the species itself is also H2O).
      // Attention, we need to treat here also the case that there
      // is no H2O species. We will then set abs_h2o to
      // -1. Absorption routines that do not really need abs_h2o
      // will still run.
      //
      // FIXME: abs_h2o is currently still needed by the continuum part.
      // Should go away eventually.
      if (h2o_index == -1) {
        // The case without H2O species.
        abs_h2o.resize(1);
        abs_h2o = -1;
      } else {
        // The normal case.
        abs_h2o = these_all_vmrs(h2o_index, joker);
      }

      // Loop temperature perturbations
      // ------------------------------

      // We use a parallel for loop for this.

      // There is something strange here: abs_lookup seems to be
      // "shared" by default, although I have set default(none). I
      // suspect that the reason for this behavior is that
      // abs_lookup is a return by reference parameter of this
      // function. Anyway, shared is the correct setting for
      // abs_lookup, so there is no problem.

  WorkspaceOmpParallelCopyGuard wss{ws};

#pragma omp parallel for if (                                         \
    !arts_omp_in_parallel() &&                                        \
    these_t_pert_nelem >=                                             \
        arts_omp_get_max_threads()) private(this_t)                   \
    firstprivate(wss)
      for (Index j = 0; j < these_t_pert_nelem; ++j) {
        // Skip remaining iterations if an error occurred
        if (failed) continue;

        // The try block here is necessary to correctly handle
        // exceptions inside the parallel region.
        try {
          if (0 != n_t_pert) {
            // We first prepare the output in a string here,
            // so that we can write it to out3 with a single
            // operation. This avoids messy output from
            // multiple threads.
            ostringstream os;

            os << "  Doing temperature variant " << j + 1 << " of " << n_t_pert
               << ": " << these_t_pert[j] << ".\n";

            out3 << os.str();
          }

          // Create perturbed temperature profile:
          this_t = abs_lookup.t_ref;
          this_t += these_t_pert[j];

          // Store in the right place:
          // Loop through all altitudes
          for (Index p = 0; p < n_p_grid; ++p) {
            PropagationMatrix K;
            StokesVector S;
            ArrayOfPropagationMatrix dK;
            ArrayOfStokesVector dS;
            Vector rtp_vmr{these_all_vmrs(joker, p)};
            for (auto& x : rtp_vmr) x = std::max(lowest_vmr, x);

            // Perform the propagation matrix computations
            propmat_clearsky_agendaExecute(wss,
                                           K,
                                           S,
                                           dK,
                                           dS,
                                           {},
                                           abs_species[i],
                                           f_grid,
                                           {},
                                           AtmPoint{},  // FIXME: DUMMY VALUE
                                           propmat_clearsky_agenda);
            K.Kjj() /= rtp_vmr[i] * number_density(abs_p[p], this_t[p]);
            abs_lookup.xsec(j, spec, Range(joker), p) = K.Kjj();

            // There used to be a division by the number density
            // n here. This is no longer necessary, since
            // abs_xsec_per_species now contains true absorption
            // cross sections.

            // IMPORTANT: There was a bug in my old Matlab
            // function "create_lookup.m" to generate the lookup
            // table. (The number density was always the
            // reference one, and did not change for different
            // temperatures.) Patricks Atmlab function
            // "arts_abstable_from_arts1.m" did *not* have this bug.

            // Calculate the number density for the given pressure and
            // temperature:
            // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
            //                  const Numeric n = number_density( abs_lookup.p_grid[p],
            //                                                    this_t[p]   );
            //                  abs_lookup.xsec( j, spec, Range(joker), p ) /= n;
          }
        }  // end of try block
        catch (const std::runtime_error& e) {
#pragma omp critical(abs_lookupCalc_fail)
          {
            fail_msg = e.what();
            failed = true;
          }
        }
      }  // end of parallel for loop

      if (failed) throw runtime_error(fail_msg);
    }
  }

  // 6. Initialize fgp_default.
  abs_lookup.flag_default = my_interp::lagrange_interpolation_list<LagrangeInterpolation>(abs_lookup.f_grid, abs_lookup.f_grid, 0);

  // Set the abs_lookup_is_adapted flag. After all, the table fits the
  // current frequency grid and species selection.
  abs_lookup_is_adapted = 1;
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
                             const ArrayOfArrayOfSpeciesTag& abs_species,
                             const Verbosity& verbosity) {
  CREATE_OUT3;

  cont.resize(0);

  // This is quite complicated, unfortunately. The approach here
  // is copied from abs_xsec_per_speciesAddConts. For explanation,
  // see there.

  // Loop tag groups:
  for (Index i = 0; i < abs_species.nelem(); ++i) {
    // Loop tags in tag group
    for (Index s = 0; s < abs_species[i].nelem(); ++s) {
      // Check for continuum tags
      if (abs_species[i][s].type == Species::TagType::Predefined ||
          abs_species[i][s].type == Species::TagType::Cia) {
        const String thisname = abs_species[i][s].Name();
        // Ok, now we know this is a continuum tag.
        out3 << "  Continuum tag: " << thisname;

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
          out3 << " --> not added.\n";
          break;
        }

        // 2. Continua known to use h2o_abs
        if ("O2-" == thisname.substr(0, 3)) {
          cont.push_back(i);
          out3 << " --> added to abs_nls.\n";
          break;
        }

        // 3. abs_species tags that are NOT allowed in LUT
        // calculations
        if ("icecloud-" == thisname.substr(0, 9) ||
            "rain-" == thisname.substr(0, 5)) {
          ostringstream os;
          os << "Tag " << thisname << " not allowed in absorption "
             << "lookup tables.";
          throw runtime_error(os.str());
        }

        // If we get here, then the tag was neither in the
        // posivitive nor in the negative list. We through a
        // runtime error.
        out3 << " --> unknown.\n";
        ostringstream os;
        os << "Unknown whether tag " << thisname
           << " is a nonlinear species (i.e. uses h2o_abs) or not.\n"
           << "Cannot set abs_nls automatically.";
        throw runtime_error(os.str());
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
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const Verbosity& verbosity) {
  CREATE_OUT2;

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
  find_nonlinear_continua(cont, abs_species, verbosity);

  // Add these to abs_nls:
  for (Index i = 0; i < cont.nelem(); ++i) {
    abs_nls.push_back(abs_species[cont[i]]);
  }

  out2 << "  Species marked for nonlinear treatment:\n";
  for (Index i = 0; i < abs_nls.nelem(); ++i) {
    out2 << "  ";
    for (Index j = 0; j < abs_nls[i].nelem(); ++j) {
      if (j != 0) out2 << ", ";
      out2 << abs_nls[i][j].Name();
    }
    out2 << "\n";
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
                       const Index& t_interp_order,
                       const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

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

  out3 << "  abs_t_pert: mindev/maxdev : " << mindev << " / " << maxdev << "\n";

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

  out2 << "  abs_t_pert: " << abs_t_pert[0] << " K to "
       << abs_t_pert[abs_t_pert.nelem() - 1] << " K in steps of "
       << effective_step << " K (" << abs_t_pert.nelem() << " grid points)\n";
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
                         const Index& nls_interp_order,
                         const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

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

  out3 << "  abs_nls_pert: mindev/maxdev : " << mindev << " / " << maxdev
       << "\n";

  bool allownegative = false;
  if (mindev < 0) {
    out2
        << "  Warning: I am getting a negative fractional distance to the H2O\n"
        << "  reference profile. Some of your H2O fields may contain negative values.\n"
        << "  Will allow negative values also for abs_nls_pert.\n";
    allownegative = true;
  }

  if (!allownegative) {
    mindev = 0;
    out3 << "  Adjusted mindev : " << mindev << "\n";
  }

  if (std::isinf(maxdev)) {
    ostringstream os;
    os << "Perturbation upper limit is infinity (likely due to the reference\n"
       << "profile being 0 at at least one pressure level). Can not work\n"
       << "with that.";
    throw runtime_error(os.str());
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
    VectorInsertGridPoints(abs_nls_pert, abs_nls_pert, {0}, verbosity);
    out2
        << "  I am including also 0 in the abs_nls_pert, because it is a turning \n"
        << "  point. Consider to use a higher abs_nls_interp_order, for example 4.\n";
  }

  out2 << "  abs_nls_pert: " << abs_nls_pert[0] << " to "
       << abs_nls_pert[abs_nls_pert.nelem() - 1]
       << " (fractional units) in steps of "
       << abs_nls_pert[1] - abs_nls_pert[0] << " (" << abs_nls_pert.nelem()
       << " grid points)\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupSetup(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Vector& abs_t_pert,
    Matrix& abs_vmrs,
    ArrayOfArrayOfSpeciesTag& abs_nls,
    Vector& abs_nls_pert,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    //                     const Vector& lat_grid,
    //                     const Vector& lon_grid,
    const Tensor3& t_field,
    const Tensor4& vmr_field,
    const Index& atmfields_checked,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Index& abs_p_interp_order,
    const Index& abs_t_interp_order,
    const Index& abs_nls_interp_order,
    // Control Parameters:
    const Numeric& p_step10,
    const Numeric& t_step,
    const Numeric& h2o_step,
    const Verbosity& verbosity) {
  // Checks on input parameters:

  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");

  // We don't actually need lat_grid and lon_grid, but we have them as
  // input variables, so that we can use the standard functions to
  // check atmospheric fields and grids. A bit cheesy, but I don't
  // want to program all the checks explicitly.

  // Check grids (outcommented the ones that have been done by
  // atmfields_checkedCalc already):
  //chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  if (p_grid.nelem() < 2) {
    ostringstream os;
    os << "You need at least two pressure levels.";
    throw runtime_error(os.str());
  }

  // Check T field:
  //chk_atm_field("t_field", t_field, atmosphere_dim,
  //              p_grid, lat_grid, lon_grid);

  // Check VMR field (and abs_species):
  //chk_atm_field("vmr_field", vmr_field, atmosphere_dim,
  //              abs_species.nelem(), p_grid, lat_grid, lon_grid);

  // Check the keyword arguments:
  if (p_step10 <= 0 || t_step <= 0 || h2o_step <= 0) {
    ostringstream os;
    os << "The keyword arguments p_step, t_step, and h2o_step must be >0.";
    throw runtime_error(os.str());
  }

  // Ok, all input parameters seem to be reasonable.

  // For consistency with other code around arts (e.g., correlation
  // lengths in atmlab), p_step is given as log10(p[Pa]). However, we
  // convert it here to the natural log:
  const Numeric p_step = log(pow(10.0, p_step10));

  // We will need the log of the pressure grid:
  Vector log_p_grid(p_grid.nelem());
  transform(log_p_grid, log, p_grid);

  //  const Numeric epsilon = 0.01 * p_step; // This is the epsilon that
  //                                         // we use for comparing p grid spacings.

  // Construct abs_p
  // ---------------

  ArrayOfNumeric log_abs_p_a;  // We take log_abs_p_a as an array of
                               // Numeric, so that we can easily
      // build it up by appending new elements to the end.

  // Check whether there are pressure levels that are further apart
  // (in log(p)) than p_step, and insert additional levels if
  // necessary:

  log_abs_p_a.push_back(log_p_grid[0]);

  for (Index i = 1; i < log_p_grid.nelem(); ++i) {
    const Numeric dp =
        log_p_grid[i - 1] - log_p_grid[i];  // The grid is descending.

    const Numeric dp_by_p_step = dp / p_step;
    //          cout << "dp_by_p_step: " << dp_by_p_step << "\n";

    // How many times does p_step fit into dp?
    const Index n = (Index)ceil(dp_by_p_step);
    // n is the number of intervals that we want to have in the
    // new grid. The number of additional points to insert is
    // n-1. But we have to insert the original point as well.
    //          cout << n << "\n";

    const Numeric ddp = dp / (Numeric)n;
    //          cout << "ddp: " << ddp << "\n";

    for (Index j = 1; j <= n; ++j)
      log_abs_p_a.push_back(log_p_grid[i - 1] - (Numeric)j * ddp);
  }

  // Copy to a proper vector, we need this also later for
  // interpolation:
  Vector log_abs_p(log_abs_p_a.nelem());
  for (Index i = 0; i < log_abs_p_a.nelem(); ++i) log_abs_p[i] = log_abs_p_a[i];

  // Copy the new grid to abs_p, removing the log:
  abs_p.resize(log_abs_p.nelem());
  transform(abs_p, exp, log_abs_p);

  // Check that abs_p has enough points for the interpolation order
  if (abs_p.nelem() < abs_p_interp_order + 1) {
    ostringstream os;
    os << "Your pressure grid does not have enough levels for the desired interpolation order.";
    throw runtime_error(os.str());
  }

  // We will also have to interpolate T and VMR profiles to the new
  // pressure grid. We interpolate in log(p), as usual in ARTS.

  // Grid positions:
  ArrayOfGridPos gp(log_abs_p.nelem());
  gridpos(gp, log_p_grid, log_abs_p);

  // Interpolation weights:
  Matrix itw(gp.nelem(), 2);
  interpweights(itw, gp);

  // In the 1D case the lookup table is just a lookup table in
  // pressure. We treat this simple case first.
  if (1 == atmosphere_dim) {
    // Reference temperature,
    // interpolate abs_t from t_field:
    abs_t.resize(log_abs_p.nelem());
    interp(abs_t, itw, t_field(joker, 0, 0), gp);

    // Temperature perturbations:
    abs_t_pert.resize(0);

    // Reference VMR profiles,
    // interpolate abs_vmrs from vmr_field:
    abs_vmrs.resize(abs_species.nelem(), log_abs_p.nelem());
    for (Index i = 0; i < abs_species.nelem(); ++i)
      interp(abs_vmrs(i, joker), itw, vmr_field(i, joker, 0, 0), gp);

    // Species for which H2O VMR is perturbed:
    abs_nls.resize(0);

    // H2O VMR perturbations:
    abs_nls_pert.resize(0);
  } else {
    // 2D or 3D case. We have to set up T and nonlinear species variations.

    // Make an intelligent choice for the nonlinear species.
    choose_abs_nls(abs_nls, abs_species, verbosity);

    // Now comes a part where we analyse the atmospheric fields.
    // We need to find the max, min, and mean profile for
    // temperature and VMRs.
    // We do this on the original p grid, not on the refined p
    // grid, to be more efficient.

    // Temperature:
    Vector tmin(p_grid.nelem());
    Vector tmax(p_grid.nelem());
    Vector tmean(p_grid.nelem());

    for (Index i = 0; i < p_grid.nelem(); ++i) {
      tmin[i] = min(t_field(i, joker, joker));
      tmax[i] = max(t_field(i, joker, joker));
      tmean[i] = mean(t_field(i, joker, joker));
    }

    //       cout << "Tmin: " << tmin << "\n";
    //       cout << "Tmax: " << tmax << "\n";
    //       cout << "Tmean: " << tmean << "\n";

    // Calculate mean profiles of all species. (We need all for abs_vmrs
    // not only H2O.)
    Matrix vmrmean(abs_species.nelem(), p_grid.nelem());
    for (Index s = 0; s < abs_species.nelem(); ++s)
      for (Index i = 0; i < p_grid.nelem(); ++i) {
        vmrmean(s, i) = mean(vmr_field(s, i, joker, joker));
      }

    // If there are NLS, determine H2O statistics:

    // We have to define these here, outside the if block, because
    // we need the values later.
    Vector h2omin(p_grid.nelem());
    Vector h2omax(p_grid.nelem());
    const Index h2o_index = find_first_species(abs_species, Species::fromShortName("H2O"));
    // We need this inside the if clauses for nonlinear species
    // treatment. The function returns "-1" if there is no H2O
    // species. There is a check for that in the next if block, with
    // an appropriate runtime error.

    if (0 < abs_nls.nelem()) {
      // Check if there really is a H2O species.
      if (h2o_index < 0) {
        ostringstream os;
        os << "Some of your species require nonlinear treatment,\n"
           << "but you have no H2O species.";
        throw runtime_error(os.str());
      }

      for (Index i = 0; i < p_grid.nelem(); ++i) {
        h2omin[i] = min(vmr_field(h2o_index, i, joker, joker));
        h2omax[i] = max(vmr_field(h2o_index, i, joker, joker));
      }

      //           cout << "H2Omin: " << h2omin << "\n";
      //           cout << "H2Omax: " << h2omax << "\n";
      //           cout << "H2Omean: " << vmrmean(h2o_index,joker) << "\n";
    }

    // Interpolate in pressure, set abs_t, abs_vmr...

    // Reference temperature,
    // interpolate abs_t from tmean:
    abs_t.resize(log_abs_p.nelem());
    interp(abs_t, itw, tmean, gp);

    // Temperature perturbations:
    choose_abs_t_pert(abs_t_pert,
                      tmean,
                      tmin,
                      tmax,
                      t_step,
                      abs_p_interp_order,
                      abs_t_interp_order,
                      verbosity);
    //       cout << "abs_t_pert: " << abs_t_pert << "\n";

    // Reference VMR profiles,
    // interpolate abs_vmrs from vmrmean:
    abs_vmrs.resize(abs_species.nelem(), log_abs_p.nelem());
    for (Index i = 0; i < abs_species.nelem(); ++i)
      interp(abs_vmrs(i, joker), itw, vmrmean(i, joker), gp);

    if (0 < abs_nls.nelem()) {
      // Construct abs_nls_pert:
      choose_abs_nls_pert(abs_nls_pert,
                          vmrmean(h2o_index, joker),
                          h2omin,
                          h2omax,
                          h2o_step,
                          abs_p_interp_order,
                          abs_nls_interp_order,
                          verbosity);
    } else {
      // Empty abs_nls_pert:
      abs_nls_pert.resize(0);
    }
    //       cout << "abs_nls_pert: " << abs_nls_pert << "\n";
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupSetupBatch(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Vector& abs_t_pert,
    Matrix& abs_vmrs,
    ArrayOfArrayOfSpeciesTag& abs_nls,
    Vector& abs_nls_pert,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfGriddedField4& batch_fields,
    const Index& abs_p_interp_order,
    const Index& abs_t_interp_order,
    const Index& abs_nls_interp_order,
    const Index& atmosphere_dim,
    // Control Parameters:
    const Numeric& p_step10,
    const Numeric& t_step,
    const Numeric& h2o_step,
    const Vector& extremes,
    const Index& robust,
    const Index& check_gridnames,
    const Verbosity& verbosity) {
  CREATE_OUT1;
  CREATE_OUT2;
  CREATE_OUT3;

  // For consistency with other code around arts (e.g., correlation
  // lengths in atmlab), p_step is given as log10(p[Pa]). However, we
  // convert it here to the natural log:
  const Numeric p_step = log(pow(10.0, p_step10));

  // Derive which abs_species is H2O (required for nonlinear species handling)
  // returns -1 if no H2O present
  const Index h2o_index = find_first_species(abs_species, Species::fromShortName("H2O"));
  //  cout << "The " << h2o_index+1 << ". species in abs_species is H2O\n";
  //  cout << "That is, H2O is expected to be the " << indoff+h2o_index
  //       << ". column of the atmospheric fields\n";

  ArrayOfIndex batch_index(abs_species.nelem());
  Index T_index = -1;
  Index z_index = -1;

  ArrayOfString species_names(abs_species.nelem());
  for (Index i = 0; i < abs_species.nelem(); ++i)
    species_names[i] = Species::toShortName(abs_species[i].Species());

  const ArrayOfString field_names = batch_fields[0].get_string_grid(0);

  String species_type;
  String species_name;
  const String delim = "-";

  // Check that the field names in batch_fields are good. (We check
  // only the first field in the batch.)
  {
    ostringstream os;
    bool bail = false;

    // One of the abs_species has to be H2O
    // (ideally that should be checked later, when we know whether there are
    // nonlinear species present. however, currently batch mode REQUIRES H2O to
    // be present, so we check that immediately)
    if (h2o_index < 0) {
      os << "One of the atmospheric fields must be H2O.\n";
      bail = true;
    }

    // First simply check dimensions of abs_species and field_names
    if (field_names.nelem() < 3) {
      os << "Atmospheric states in batch_fields must have at\n"
         << "least three fields: T, z, and at least one absorption species.";
      throw runtime_error(os.str());
    }

    if (abs_species.nelem() < 1) {
      os << "At least one absorption species needs to be defined "
         << "in abs_species.";
      throw runtime_error(os.str());
    }

    // Check that all required fields are present.
    const Index nf = field_names.nelem();
    bool found;

    // Looking for temperature field:
    found = false;
    for (Index fi = 0; fi < nf; ++fi) {
      parse_atmcompact_speciestype(species_type, field_names[fi], delim);
      if (species_type == "T") {
        if (found) {
          os << "Only one temperature ('T') field allowed, "
             << "but found at least 2.\n";
          bail = true;
        } else {
          found = true;
          T_index = fi;
        }
      }
    }
    if (!found) {
      os << "One temperature ('T') field required, but none found.\n";
      bail = true;
    }

    // Looking for altitude field:
    found = false;
    for (Index fi = 0; fi < nf; ++fi) {
      parse_atmcompact_speciestype(species_type, field_names[fi], delim);
      if (species_type == "z") {
        if (found) {
          os << "Only one altitude ('z') field allowed, "
             << "but found at least 2.\n";
          bail = true;
        } else {
          found = true;
          z_index = fi;
        }
      }
    }
    if (!found) {
      os << "One altitude ('z') field required, but none found.\n";
      bail = true;
    }

    // Now going over all abs_species elements and match the corresponding
    // field_name (we always take the first occurence of a matching field). We
    // don't care that batch_fields contains further fields, too. And we can't
    // expect the fields being in the same order as the abs_species.
    // At the same time, we keep the indices of the fields corresponding to the
    // abs_species elements, because later we will need the field data.
    Index fi;
    for (Index j = 0; j < abs_species.nelem(); ++j) {
      found = false;
      fi = 0;
      while (!found && fi < nf) {
        parse_atmcompact_speciestype(species_type, field_names[fi], delim);
        // do we have an abs_species type field?
        if (species_type == "abs_species") {
          parse_atmcompact_speciesname(species_name, field_names[fi], delim);
          if (species_name == species_names[j]) {
            found = true;
            batch_index[j] = fi;
          }
        }
        fi++;
      }
      if (!found) {
        os << "No field for absorption species '" << species_names[j]
           << "'found.\n";
        bail = true;
      }
    }

    os << "Your field names are:\n" << field_names;

    if (bail) throw runtime_error(os.str());
  }

  // FIXME: Adjustment of min/max values for Jacobian perturbations is still missing.

  // Make an intelligent choice for the nonlinear species.
  choose_abs_nls(abs_nls, abs_species, verbosity);

  // Find out maximum and minimum pressure and check that pressure grid is decreasing.
  Numeric maxp = batch_fields[0].get_numeric_grid(GFIELD4_P_GRID)[0];
  Numeric minp = batch_fields[0].get_numeric_grid(
      GFIELD4_P_GRID)[batch_fields[0].get_numeric_grid(GFIELD4_P_GRID).nelem() -
                      1];

  ArrayOfIndex valid_field_indices;
  for (Index i = 0; i < batch_fields.nelem(); ++i) {
    // Local variables for atmfields_check.
    Index atmfields_checked;
    Tensor4 t4_dummy;
    Tensor3 t3_dummy;

    Vector p_grid;
    Vector lat_grid;
    Vector lon_grid;
    Tensor3 t_field;
    Tensor3 z_field;
    Tensor4 vmr_field;
    Tensor4 particle_bulkprop_field;
    ArrayOfString particle_bulkprop_names;
    GriddedField4 atm_fields_compact;
    Index abs_f_interp_order{0};

    // Extract fields from atmfield and check their validity.
    // This closes the loophole when only calculating lookup tables.
    atm_fields_compact = batch_fields[i];

    AtmFieldsAndParticleBulkPropFieldFromCompact(p_grid,
                                                 lat_grid,
                                                 lon_grid,
                                                 t_field,
                                                 z_field,
                                                 vmr_field,
                                                 particle_bulkprop_field,
                                                 particle_bulkprop_names,
                                                 abs_species,
                                                 atm_fields_compact,
                                                 atmosphere_dim,
                                                 "-",
                                                 0,
                                                 check_gridnames,
                                                 verbosity);

    try {
      atmfields_checkedCalc(atmfields_checked,
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid,
                            abs_species,
                            t_field,
                            vmr_field,
                            t3_dummy,
                            t3_dummy,
                            t3_dummy,
                            t3_dummy,
                            t3_dummy,
                            t3_dummy,
                            abs_f_interp_order,
                            0,
                            verbosity);
    } catch (const std::exception& e) {
      // If `robust`, skip field and continue, ...
      if (robust) {
        out1 << "  WARNING! Skipped invalid atmfield "
             << "at batch_atmfield index " << i << ".\n"
             << "The runtime error produced was:\n"
             << e.what() << "\n";
        continue;
      }
      // ... else throw an error.
      else {
        stringstream err;
        err << "Invalid atmfield at batch_atmfield index " << i << ".\n"
            << "The runtime error produced was:\n"
            << e.what() << "\n";
        throw std::runtime_error(err.str());
      }
    };
    valid_field_indices.push_back(i);  // Append index to list of valid fields.

    if (maxp < batch_fields[i].get_numeric_grid(GFIELD4_P_GRID)[0])
      maxp = batch_fields[i].get_numeric_grid(GFIELD4_P_GRID)[0];
    if (minp >
        batch_fields[i].get_numeric_grid(GFIELD4_P_GRID)
            [batch_fields[i].get_numeric_grid(GFIELD4_P_GRID).nelem() - 1])
      minp = batch_fields[i].get_numeric_grid(GFIELD4_P_GRID)
                 [batch_fields[i].get_numeric_grid(GFIELD4_P_GRID).nelem() - 1];
  }
  //  cout << "  minp/maxp: " << minp << " / " << maxp << "\n";

  // Information on the number of skipped atmospheres.
  if (batch_fields.nelem() > valid_field_indices.nelem()) {
    out1 << "  " << batch_fields.nelem() - valid_field_indices.nelem()
         << " out of " << batch_fields.nelem() << " atmospheres ignored.\n";
  }

  // Throw error if no atmfield passed the check.
  if (valid_field_indices.nelem() < 1) {
    stringstream err;
    err << "You need at least one valid profile.\n"
        << "It seems that no atmfield passed the checks!\n";
    throw std::runtime_error(err.str());
  }

  if (maxp == minp) {
    ostringstream os;
    os << "You need at least two pressure levels.";
    throw runtime_error(os.str());
  }

  // We construct the pressure grid as follows:
  // - Everything is done in log(p).
  // - Start with maxp and go down in steps of p_step until we are <= minp.
  // - Adjust the final pressure value to be exactly minp, otherwise
  //   we have problems in getting min, max, and mean values for this
  //   point later.
  Index np = (Index)ceil((log(maxp) - log(minp)) / p_step) + 1;
  // The +1 above has to be there because we must include also both
  // grid end points.

  // If np is too small for the interpolation order, we increase it:
  if (np < abs_p_interp_order + 1) np = abs_p_interp_order + 1;

  Vector log_abs_p=uniform_grid(log(maxp), np, -p_step);
  log_abs_p[np - 1] = log(minp);

  abs_p.resize(np);
  transform(abs_p, exp, log_abs_p);
  out2 << "  abs_p: " << abs_p[0] << " Pa to " << abs_p[np - 1]
       << " Pa in log10 steps of " << p_step10 << " (" << np
       << " grid points)\n";

  // Now we have to determine the statistics of T and VMRs, we need
  // profiles of min, max, and mean of these, on the abs_p grid.

  //  Index n_variables = batch_fields[0].data.nbooks();
  // we will do statistics for data fields excluding the scat_species fields
  Index n_variables = 2 + abs_species.nelem();

  // The first dimension of datamin, datamax, and datamean is the
  // variable (T,Z,H2O,O3,...). The second dimension is pressure. We
  // assume all elements of the batch have the same variables.

  Matrix datamin(n_variables, np, numeric_limits<Numeric>::max());
  Matrix datamax(n_variables, np, numeric_limits<Numeric>::min());
  // The limits here are from the header file <limits>
  Matrix datamean(n_variables, np, 0);
  Vector mean_norm(np, 0);  // Divide by this to get mean.

  // We will loop over all batch cases to analyze the statistics and
  // calculate reference profiles. As a little side project, we will
  // also calculate the absolute min and max of temperature and
  // humidity. This is handy as input to abs_lookupSetupWide.
  Numeric mint = +1e99;
  Numeric maxt = -1e99;
  Numeric minh2o = +1e99;
  Numeric maxh2o = -1e99;

  // Loop over valid atmfields.
  for (Index vi = 0; vi < valid_field_indices.nelem(); ++vi) {
    // Get internal field index.
    Index i = valid_field_indices[vi];

    // Check that really each case has the same variables (see
    // comment above.)
    if (batch_fields[i].get_string_grid(GFIELD4_FIELD_NAMES) !=
        batch_fields[0].get_string_grid(GFIELD4_FIELD_NAMES))
      throw runtime_error(
          "All batch atmospheres must contain the same field names.");

    for (Index j = 0;
         j < batch_fields[i].get_string_grid(GFIELD4_FIELD_NAMES).nelem();
         ++j)
      if (batch_fields[i].get_string_grid(GFIELD4_FIELD_NAMES)[j] !=
          batch_fields[0].get_string_grid(GFIELD4_FIELD_NAMES)[j])
        throw runtime_error(
            "All batch atmospheres must contain the same individual field names.");

    // Create convenient handles:
    const Vector& p_grid = batch_fields[i].get_numeric_grid(GFIELD4_P_GRID);
    const Tensor4& data = batch_fields[i].data;

    // Update our global max/min values for T and H2O:

    // We have to loop over pressure, latitudes, and longitudes
    // here. The dimensions of data are:
    // data[N_fields][N_p][N_lat][N_lon]

    for (Index ip = 0; ip < data.npages(); ++ip)
      for (Index ilat = 0; ilat < data.nrows(); ++ilat)
        for (Index ilon = 0; ilon < data.ncols(); ++ilon) {
          // Field T_index is temperature:
          if (data(T_index, ip, ilat, ilon) < mint)
            mint = data(T_index, ip, ilat, ilon);
          if (data(T_index, ip, ilat, ilon) > maxt)
            maxt = data(T_index, ip, ilat, ilon);
          // Field batch_index[h2o_index] is H2O:
          if (data(batch_index[h2o_index], ip, ilat, ilon) < minh2o) {
            minh2o = data(batch_index[h2o_index], ip, ilat, ilon);
          }
          if (data(batch_index[h2o_index], ip, ilat, ilon) > maxh2o) {
            maxh2o = data(batch_index[h2o_index], ip, ilat, ilon);
          }
        }

    // Interpolate the current batch fields to the abs_p grid. We
    // have to do this for each batch case, since the grids could
    // all be different.

    Vector log_p_grid(p_grid.nelem());
    transform(log_p_grid, log, p_grid);

    // There is a catch here: We can only use the part of abs_p that
    // is inside the current pressure grid p_grid, otherwise we
    // would have to extrapolate.
    // The eps_bottom and eps_top account for the little bit of
    // extrapolation that is allowed.

    const Numeric eps_bottom = (log_p_grid[0] - log_p_grid[1]) / 2.1;
    Index first_p = 0;
    while (log_abs_p[first_p] > log_p_grid[0] + eps_bottom) ++first_p;

    const Numeric eps_top = (log_p_grid[log_p_grid.nelem() - 2] -
                             log_p_grid[log_p_grid.nelem() - 1]) /
                            2.1;
    Index last_p = log_abs_p.nelem() - 1;
    while (log_abs_p[last_p] < log_p_grid[log_p_grid.nelem() - 1] - eps_top)
      --last_p;

    const Index extent_p = last_p - first_p + 1;

    // This was too complicated to get right:
    //      const Index first_p   = (Index) round ( (log_abs_p[0]       - log_p_grid[0])                    / p_step);
    //      const Index extent_p  = (Index) round ( (log_abs_p[first_p] - log_p_grid[log_p_grid.nelem()-1]) / p_step) + 1;

    ConstVectorView active_log_abs_p = log_abs_p[Range(first_p, extent_p)];

    //       cout << "first_p / last_p / extent_p : " << first_p << " / " << last_p << " / " << extent_p << "\n";
    //       cout << "log_p_grid: "        << log_p_grid << "\n";
    //       cout << "log_abs_p:  "        << log_abs_p << "\n";
    //       cout << "active_log_abs_p:  " << active_log_abs_p << "\n";
    //       cout << "=============================================================\n";
    //       arts_exit();

    // Grid positions:
    ArrayOfGridPos gp(active_log_abs_p.nelem());
    gridpos(gp, log_p_grid, active_log_abs_p);
    //      gridpos(gp, log_p_grid, active_log_abs_p, 100);
    // We allow much more extrapolation here than normal (0.5 is
    // normal). If we do not do this, then we get problems for
    // p_grids that are much finer than abs_p.

    // Interpolation weights:
    Matrix itw(gp.nelem(), 2);
    interpweights(itw, gp);

    // We have to loop over fields, latitudes, and longitudes here, doing the
    // pressure interpolation for all. The dimensions of data are:
    // data[N_fields][N_p][N_lat][N_lon]
    // For data_interp we reduce data by the particle related fields, i.e.,
    // the ones in between the first two and the last N=abs_species.nelem().
    //      Tensor4 data_interp(data.nbooks(),
    Tensor4 data_interp(
        n_variables, active_log_abs_p.nelem(), data.nrows(), data.ncols());

    for (Index lo = 0; lo < data.ncols(); ++lo)
      for (Index la = 0; la < data.nrows(); ++la) {
        //          for (Index fi=0; fi<data.nbooks(); ++fi)
        // we have to handle T/z and abs_species parts separately since they
        // can have scat_species in between, which we want to ignore
        interp(data_interp(0, joker, la, lo),
               itw,
               data(T_index, joker, la, lo),
               gp);
        interp(data_interp(1, joker, la, lo),
               itw,
               data(z_index, joker, la, lo),
               gp);
        for (Index fi = 0; fi < abs_species.nelem(); ++fi)
          interp(data_interp(fi + 2, joker, la, lo),
                 itw,
                 data(batch_index[fi], joker, la, lo),
                 gp);
      }

    // Now update our datamin, datamax, and datamean variables.
    // We need the min and max only for the T and H2O profile,
    // not for others. But we need the mean for all. We are just
    // hopping over the case that we do not need below. This is not
    // very clean, but efficient. And it avoids handling all the
    // different cases explicitly.
    for (Index lo = 0; lo < data_interp.ncols(); ++lo)
      for (Index la = 0; la < data_interp.nrows(); ++la) {
        for (Index fi = 0; fi < data_interp.nbooks(); ++fi) {
          if (1 != fi)  // We skip the z field, which we do not need
            for (Index pr = 0; pr < data_interp.npages(); ++pr) {
              // Min and max only needed for T and H2o
              if (0 == fi || (h2o_index + 2) == fi) {
                if (data_interp(fi, pr, la, lo) < datamin(fi, first_p + pr))
                  datamin(fi, first_p + pr) = data_interp(fi, pr, la, lo);
                if (data_interp(fi, pr, la, lo) > datamax(fi, first_p + pr))
                  datamax(fi, first_p + pr) = data_interp(fi, pr, la, lo);
              }

              datamean(fi, first_p + pr) += data_interp(fi, pr, la, lo);
            }
        }

        // The mean_norm is actually a bit tricky. It depends on
        // pressure, since different numbers of cases contribute
        // to the mean for different pressures. At the very eges
        // of the grid, typically only a single case contributes.

        mean_norm[Range(first_p, extent_p)] += 1;
      }
  }

  out2 << "  Global statistics:\n"
       << "  min(p)   / max(p)   [Pa]:  " << minp << " / " << maxp << "\n"
       << "  min(T)   / max(T)   [K]:   " << mint << " / " << maxt << "\n"
       << "  min(H2O) / max(H2O) [VMR]: " << minh2o << " / " << maxh2o << "\n";

  // Divide mean by mean_norm to get the mean:
  ARTS_ASSERT(np == mean_norm.nelem());
  for (Index fi = 0; fi < datamean.nrows(); ++fi)
    if (1 != fi)  // We skip the z field, which we do not need
      for (Index pi = 0; pi < np; ++pi) {
        // We do this in an explicit loop here, since we have to
        // check whether there really were data points to calculate
        // the mean at each level.
        if (0 < mean_norm[pi])
          datamean(fi, pi) /= mean_norm[pi];
        else {
          ostringstream os;
          os << "No data at pressure level " << pi + 1 << " of " << np << " ("
             << abs_p[pi] << " hPa).";
          throw runtime_error(os.str());
        }
      }

  // If we do higher order pressure interpolation, then we should
  // smooth the reference profiles with a boxcar filter of width
  // p_interp_order+1. Otherwise we get numerical problems if there
  // are any sharp spikes in the reference profiles.
  ARTS_ASSERT(log_abs_p.nelem() == np);
  Matrix smooth_datamean(datamean.nrows(), datamean.ncols(), 0);
  for (Index i = 0; i < np; ++i) {
    const Index idx0 = my_interp::pos_finder<false>(i, log_abs_p[i], log_abs_p, abs_p_interp_order);

    for (Index fi = 0; fi < datamean.nrows(); ++fi)
      if (1 != fi)  // We skip the z field, which we do not need
      {
        for (Index j = 0; j < abs_p_interp_order+1; ++j) {
          smooth_datamean(fi, i) += datamean(fi, idx0 + j);
        }
        smooth_datamean(fi, i) /= Numeric(abs_p_interp_order + 1);
      }
    //       cout << "H2O-raw / H2O-smooth: "
    //            << datamean(h2o_index+2,i) << " / "
    //            << smooth_datamean(h2o_index+2,i) << "\n";
  }

  // There is another complication: If the (smoothed) mean for the H2O
  // reference profile is 0, then we have to adjust both mean and max
  // value to a non-zero number, otherwise the reference profile will
  // be zero, and we will get numerical problems later on when we
  // divide by the reference value. So, we set it here to 1e-9.
  for (Index i = 0; i < np; ++i) {
    // Assert that really H2O has index batch_index[h2o_index] VMR field list
    // and h2o_index in abs_species
    parse_atmcompact_speciestype(
        species_type, field_names[batch_index[h2o_index]], delim);
    ARTS_ASSERT(species_type == "abs_species");
    parse_atmcompact_speciesname(
        species_name, field_names[batch_index[h2o_index]], delim);
    ARTS_ASSERT("H2O" == species_name);
    ARTS_ASSERT("H2O" == species_names[h2o_index]);

    // Find mean and max H2O for this level:
    Numeric& mean_h2o = smooth_datamean(h2o_index + 2, i);
    Numeric& max_h2o = datamax(h2o_index + 2, i);
    if (mean_h2o <= 0) {
      mean_h2o = 1e-9;
      max_h2o = 1e-9;
      out3 << "  H2O profile contained zero values, adjusted to 1e-9.\n";
    }
  }

  // Set abs_t:
  abs_t.resize(np);
  abs_t = smooth_datamean(0, joker);
  //   cout << "abs_t: " << abs_t << "\n\n";
  //   cout << "tmin:  " << datamin(0,joker) << "\n\n";
  //   cout << "tmax:  " << datamax(0,joker) << "\n";

  // Set abs_vmrs:
  ARTS_ASSERT(abs_species.nelem() == smooth_datamean.nrows() - 2);
  abs_vmrs.resize(abs_species.nelem(), np);
  abs_vmrs = smooth_datamean(Range(2, abs_species.nelem()), joker);
  //  cout << "\n\nabs_vmrs: " << abs_vmrs << "\n\n";

  // Construct abs_t_pert:
  ConstVectorView tmin = datamin(0, joker);
  ConstVectorView tmax = datamax(0, joker);
  choose_abs_t_pert(abs_t_pert,
                    abs_t,
                    tmin,
                    tmax,
                    t_step,
                    abs_p_interp_order,
                    abs_t_interp_order,
                    verbosity);
  //  cout << "abs_t_pert: " << abs_t_pert << "\n";

  // Construct abs_nls_pert:
  ConstVectorView h2omin = datamin(h2o_index + 2, joker);
  ConstVectorView h2omax = datamax(h2o_index + 2, joker);
  choose_abs_nls_pert(abs_nls_pert,
                      abs_vmrs(h2o_index, joker),
                      h2omin,
                      h2omax,
                      h2o_step,
                      abs_p_interp_order,
                      abs_nls_interp_order,
                      verbosity);
  //  cout << "abs_nls_pert: " << abs_nls_pert << "\n";

  // Append the explicitly given user extreme values, if necessary:
  if (0 != extremes.nelem()) {
    // There must be 4 values in this case: t_min, t_max, h2o_min, h2o_max
    if (4 != extremes.nelem()) {
      ostringstream os;
      os << "There must be exactly 4 elements in extremes:\n"
         << "min(abs_t_pert), max(abs_t_pert), min(abs_nls_pert), max(abs_nls_pert)";
      throw runtime_error(os.str());
    }

    // t_min:
    if (extremes[0] < abs_t_pert[0]) {
      Vector dummy = abs_t_pert;
      abs_t_pert.resize(abs_t_pert.nelem() + 1);
      abs_t_pert[0] = extremes[0];
      abs_t_pert[Range(1, dummy.nelem())] = dummy;
      out2 << "  Added min extreme value for abs_t_pert: " << abs_t_pert[0]
           << "\n";
    }

    // t_max:
    if (extremes[1] > abs_t_pert[abs_t_pert.nelem() - 1]) {
      Vector dummy = abs_t_pert;
      abs_t_pert.resize(abs_t_pert.nelem() + 1);
      abs_t_pert[Range(0, dummy.nelem())] = dummy;
      abs_t_pert[abs_t_pert.nelem() - 1] = extremes[1];
      out2 << "  Added max extreme value for abs_t_pert: "
           << abs_t_pert[abs_t_pert.nelem() - 1] << "\n";
    }

    // nls_min:
    if (extremes[2] < abs_nls_pert[0]) {
      Vector dummy = abs_nls_pert;
      abs_nls_pert.resize(abs_nls_pert.nelem() + 1);
      abs_nls_pert[0] = extremes[2];
      abs_nls_pert[Range(1, dummy.nelem())] = dummy;
      out2 << "  Added min extreme value for abs_nls_pert: " << abs_nls_pert[0]
           << "\n";
    }

    // nls_max:
    if (extremes[3] > abs_nls_pert[abs_nls_pert.nelem() - 1]) {
      Vector dummy = abs_nls_pert;
      abs_nls_pert.resize(abs_nls_pert.nelem() + 1);
      abs_nls_pert[Range(0, dummy.nelem())] = dummy;
      abs_nls_pert[abs_nls_pert.nelem() - 1] = extremes[3];
      out2 << "  Added max extreme value for abs_nls_pert: "
           << abs_nls_pert[abs_nls_pert.nelem() - 1] << "\n";
    }
  }
}

void abs_lookupSetupWide(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Vector& abs_t_pert,
    Matrix& abs_vmrs,
    ArrayOfArrayOfSpeciesTag& abs_nls,
    Vector& abs_nls_pert,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Index& abs_p_interp_order,
    const Index& abs_t_interp_order,
    const Index& abs_nls_interp_order,
    // Control Parameters:
    const Numeric& p_min,
    const Numeric& p_max,
    const Numeric& p_step10,
    const Numeric& t_min,
    const Numeric& t_max,
    const Numeric& h2o_min,
    const Numeric& h2o_max,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  // For consistency with other code around arts (e.g., correlation
  // lengths in atmlab), p_step is given as log10(p[Pa]). However, we
  // convert it here to the natural log:
  const Numeric p_step = log(pow(10.0, p_step10));

  // Make an intelligent choice for the nonlinear species.
  choose_abs_nls(abs_nls, abs_species, verbosity);

  // 1. Fix pressure grid abs_p
  // --------------------------

  // We construct the pressure grid as follows:
  // - Everything is done in log(p).
  // - Start with p_max and go down in steps of p_step until we are <= p_min.

  Index np = (Index)ceil((log(p_max) - log(p_min)) / p_step) + 1;
  // The +1 above has to be there because we must include also both
  // grid end points.

  // If np is too small for the interpolation order, we increase it:
  if (np < abs_p_interp_order + 1) np = abs_p_interp_order + 1;

  Vector log_abs_p=uniform_grid(log(p_max), np, -p_step);

  abs_p.resize(np);
  transform(abs_p, exp, log_abs_p);
  out2 << "  abs_p: " << abs_p[0] << " Pa to " << abs_p[np - 1]
       << " Pa in log10 steps of " << p_step10 << " (" << np
       << " grid points)\n";

  // 2. Fix reference temperature profile abs_t and temperature perturbations
  // ------------------------------------------------------------------------

  // We simply take a constant reference profile.

  Numeric t_ref = (t_min + t_max) / 2;

  abs_t.resize(np);
  abs_t = t_ref;  // Assign same value to all elements.

  // We have to make vectors out of t_min and t_max, so we can use
  // them in the choose_abs_t_pert function call:
  Vector min_prof(np), max_prof(np);
  min_prof = t_min;
  max_prof = t_max;

  // Chose temperature perturbations:
  choose_abs_t_pert(abs_t_pert,
                    abs_t,
                    min_prof,
                    max_prof,
                    20,
                    abs_p_interp_order,
                    abs_t_interp_order,
                    verbosity);

  // 3. Fix reference H2O profile and abs_nls_pert
  // ---------------------------------------------

  // We take a constant reference profile of 1000ppm (=1e-3) for H2O
  Numeric const h2o_ref = 1e-3;

  // And 1 ppt (1e-9) as default for all VMRs
  Numeric const other_ref = 1e-9;

  // We have to assign this value to all pressures of the H2O profile,
  // and 0 to all other profiles.

  // abs_vmrs has dimension [n_species, np].
  abs_vmrs.resize(abs_species.nelem(), np);
  abs_vmrs = other_ref;

  // We look for O2 and N2, and assign constant values to them.
  // The values are from Wallace&Hobbs, 2nd edition.

  const Index o2_index =
      find_first_species(abs_species, Species::fromShortName("O2"));
  if (o2_index >= 0) {
    abs_vmrs(o2_index, joker) = 0.2095;
  }

  const Index n2_index =
      find_first_species(abs_species, Species::fromShortName("N2"));
  if (n2_index >= 0) {
    abs_vmrs(n2_index, joker) = 0.7808;
  }

  // Which species is H2O?
  const Index h2o_index = find_first_species(
      abs_species, Species::fromShortName("H2O"));

  // The function returns "-1" if there is no H2O
  // species.
  if (0 < abs_nls.nelem()) {
    if (h2o_index < 0) {
      ostringstream os;
      os << "Some of your species require nonlinear treatment,\n"
         << "but you have no H2O species.";
      throw runtime_error(os.str());
    }

    // Assign constant reference value to all H2O levels:
    abs_vmrs(h2o_index, joker) = h2o_ref;

    // We have to make vectors out of h2o_min and h2o_max, so we can use
    // them in the choose_abs_nls_pert function call.
    // We re-use the vectors min_prof and max_prof that we have
    // defined above.
    min_prof = h2o_min;
    max_prof = h2o_max;

    // Construct abs_nls_pert:
    choose_abs_nls_pert(abs_nls_pert,
                        abs_vmrs(h2o_index, joker),
                        min_prof,
                        max_prof,
                        1e99,
                        abs_p_interp_order,
                        abs_nls_interp_order,
                        verbosity);
  } else {
    CREATE_OUT1;
    out1 << "  WARNING:\n"
         << "  You have no species that require H2O variations.\n"
         << "  This case might work, but it has never been tested.\n"
         << "  Please test it, then remove this warning.\n";
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const ArrayOfString& names,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  // Size of initial array
  Index n_gs = abs_species.nelem();

  // Temporary ArrayOfSpeciesTag
  ArrayOfSpeciesTag temp;

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for (Index i = 0; i < names.nelem(); ++i) {
    abs_species.emplace_back(names[i]);
  }

  check_abs_species(abs_species);

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Added tag groups:";
  for (Index i = n_gs; i < abs_species.nelem(); ++i) {
    out3 << "\n  " << i << ":";
    for (Index s = 0; s < abs_species[i].nelem(); ++s) {
      out3 << " " << abs_species[i][s].Name();
    }
  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd2(  // WS Output:
    Workspace& ws,
    ArrayOfArrayOfSpeciesTag& abs_species,
    ArrayOfRetrievalQuantity& jq,
    Agenda& jacobian_agenda,
    Index& propmat_clearsky_agenda_checked,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // WS Generic Input:
    const Vector& rq_p_grid,
    const Vector& rq_lat_grid,
    const Vector& rq_lon_grid,
    // Control Parameters:
    const String& species,
    const String& mode,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  // Add species to *abs_species*
  abs_species.emplace_back(species);

  check_abs_species(abs_species);

  // Print list of added tag group to the most verbose output stream:
  out3 << "  Appended tag group:";
  out3 << "\n  " << abs_species.nelem() - 1 << ":";
  for (Index s = 0; s < abs_species.back().nelem(); ++s) {
    out3 << " " << abs_species.back()[s].Name();
  }
  out3 << '\n';

  // Do retrieval part
  jacobianAddAbsSpecies(ws,
                        jq,
                        jacobian_agenda,
                        atmosphere_dim,
                        p_grid,
                        lat_grid,
                        lon_grid,
                        rq_p_grid,
                        rq_lat_grid,
                        rq_lon_grid,
                        species,
                        mode,
                        1,
                        verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesInit(ArrayOfArrayOfSpeciesTag& abs_species, const Verbosity&) {
  abs_species.resize(0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesSet(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const ArrayOfString& names,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  abs_species.resize(names.nelem());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for (Index i = 0; i < names.nelem(); ++i) {
    // This part has now been moved to array_species_tag_from_string.
    // Call this function.
    abs_species[i] = ArrayOfSpeciesTag(names[i]);
  }

  check_abs_species(abs_species);

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups: ";
  for (Index i = 0; i < abs_species.nelem(); ++i) {
    out3 << "\n  " << i << ":";
    for (Index s = 0; s < abs_species[i].nelem(); ++s)
      out3 << " " << abs_species[i][s].Name();
  }
  out3 << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupAdapt(GasAbsLookup& abs_lookup,
                     Index& abs_lookup_is_adapted,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const Vector& f_grid,
                     const Verbosity& verbosity) {
  abs_lookup.Adapt(abs_species, f_grid, verbosity);
  abs_lookup_is_adapted = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFromLookup(
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
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
    const Index& no_negatives,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Variables needed by abs_lookup.Extract:
  Matrix abs_scalar_gas, dabs_scalar_gas_df, dabs_scalar_gas_dt;

  // Check if the table has been adapted:
  if (1 != abs_lookup_is_adapted)
    throw runtime_error(
        "Gas absorption lookup table must be adapted,\n"
        "use method abs_lookupAdapt.");

  const bool do_jac = supports_lookup(jacobian_quantities);
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);

  const Vector a_vmr_list = [&]() {
    Vector vmr(abs_species.nelem());
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
      propmat_clearsky.Kjj() += abs_scalar_gas(ii, joker);
    }
  } else {
    for (Index isp = 0; isp < abs_scalar_gas.nrows(); isp++) {
      propmat_clearsky.Kjj() += abs_scalar_gas(isp, joker);

      for (Index iv = 0; iv < propmat_clearsky.NumberOfFrequencies();
           iv++) {
        for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
          const auto& deriv = jacobian_quantities[iq];
          
          if (not deriv.propmattype()) continue;
          
          if (deriv == Jacobian::Atm::Temperature) {
            dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                (dabs_scalar_gas_dt(isp, iv) - abs_scalar_gas(isp, iv)) / dt;
          } else if (is_frequency_parameter(deriv)) {
            dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                (dabs_scalar_gas_df(isp, iv) - abs_scalar_gas(isp, iv)) / df;
          } else if (deriv == abs_species[isp]) {
            // WARNING:  If CIA in list, this scales wrong by factor 2
              dpropmat_clearsky_dx[iq].Kjj()[iv] += (std::isnormal(a_vmr_list[isp])) ? abs_scalar_gas(isp, iv) / a_vmr_list[isp] : 0;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearsky_fieldCalc(Workspace& ws,
                                // WS Output:
                                Tensor7& propmat_clearsky_field,
                                // WS Input:
                                const Index& atmfields_checked,
                                const Vector& f_grid,
                                const Index& stokes_dim,
                                const Vector& p_grid,
                                const Vector& lat_grid,
                                const Vector& lon_grid,
                                const Tensor3& t_field,
                                const Tensor4& vmr_field,
                                const EnergyLevelMap& nlte_field,
                                const Tensor3& mag_u_field,
                                const Tensor3& mag_v_field,
                                const Tensor3& mag_w_field,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Agenda& abs_agenda,
                                // WS Generic Input:
                                const Vector& doppler,
                                const Vector& los,
                                const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");

  ArrayOfPropagationMatrix partial_abs;
  ArrayOfStokesVector partial_nlte;  // FIXME: This is not stored!
  PropagationMatrix abs;
  StokesVector nlte;
  Vector a_vmr_list;
  EnergyLevelMap a_nlte_list;

  // Get the number of species from the leading dimension of vmr_field:
  const Index n_species = vmr_field.nbooks();

  // Number of frequencies:
  const Index n_frequencies = f_grid.nelem();

  // Number of pressure levels:
  const Index n_pressures = p_grid.nelem();

  // Number of latitude grid points (must be at least one):
  const Index n_latitudes = max(Index(1), lat_grid.nelem());

  // Number of longitude grid points (must be at least one):
  const Index n_longitudes = max(Index(1), lon_grid.nelem());

  // Check that doppler is empty or matches p_grid
  if (0 != doppler.nelem() && p_grid.nelem() != doppler.nelem()) {
    ostringstream os;
    os << "Variable doppler must either be empty, or match the dimension of "
       << "p_grid.";
    throw runtime_error(os.str());
  }

  // Resize output field.
  // The dimension in lat and lon must be at least one, even if these
  // grids are empty.
  out2 << "  Creating propmat field with dimensions:\n"
       << "    " << n_species << "   gas species,\n"
       << "    " << n_frequencies << "   frequencies,\n"
       << "    " << stokes_dim << "   stokes dimension,\n"
       << "    " << stokes_dim << "   stokes dimension,\n"
       << "    " << n_pressures << "   pressures,\n"
       << "    " << n_latitudes << "   latitudes,\n"
       << "    " << n_longitudes << "   longitudes.\n";

  propmat_clearsky_field.resize(n_species,
                                n_frequencies,
                                stokes_dim,
                                stokes_dim,
                                n_pressures,
                                n_latitudes,
                                n_longitudes);
  
  // Fake jacobian_quantities to deal with partial absorption
  ArrayOfRetrievalQuantity jacobian_quantities;
  for (auto& species_list: abs_species) {
    jacobian_quantities.emplace_back();
    auto& rq = jacobian_quantities.back();
    rq.Subtag(species_list.Name());
    rq.Target(Jacobian::Target(Jacobian::Special::ArrayOfSpeciesTagVMR, species_list));
    rq.Target().perturbation = 0.001;
  }

  String fail_msg;
  bool failed = false;

  // Make local copy of f_grid, so that we can apply Dopler if we want.
  Vector this_f_grid = f_grid;

  WorkspaceOmpParallelCopyGuard wss{ws};

  // Now we have to loop all points in the atmosphere:
  if (n_pressures)
#pragma omp parallel for if (!arts_omp_in_parallel() &&                 \
                             n_pressures >= arts_omp_get_max_threads()) \
    firstprivate(wss, this_f_grid) private(                              \
        abs, nlte, partial_abs, partial_nlte, a_vmr_list)
    for (Index ipr = 0; ipr < n_pressures; ++ipr)  // Pressure:  ipr
    {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      // The try block here is necessary to correctly handle
      // exceptions inside the parallel region.
      try {
        Numeric a_pressure = p_grid[ipr];

        if (0 != doppler.nelem()) {
          this_f_grid = f_grid;
          this_f_grid += doppler[ipr];
        }

        {
          ostringstream os;
          os << "  p_grid[" << ipr << "] = " << a_pressure << "\n";
          out3 << os.str();
        }

        for (Index ila = 0; ila < n_latitudes; ++ila)     // Latitude:  ila
          for (Index ilo = 0; ilo < n_longitudes; ++ilo)  // Longitude: ilo
          {
            Numeric a_temperature = t_field(ipr, ila, ilo);
            a_vmr_list = vmr_field(Range(joker), ipr, ila, ilo);
            if (!nlte_field.value.empty())
              a_nlte_list = nlte_field(ipr, ila, ilo);

            Vector this_rtp_mag(3, 0.);

            if (mag_u_field.npages() != 0) {
              this_rtp_mag[0] = mag_u_field(ipr, ila, ilo);
            }
            if (mag_v_field.npages() != 0) {
              this_rtp_mag[1] = mag_v_field(ipr, ila, ilo);
            }
            if (mag_w_field.npages() != 0) {
              this_rtp_mag[2] = mag_w_field(ipr, ila, ilo);
            }

            // Execute agenda to calculate local absorption.
            // Agenda input:  f_index, a_pressure, a_temperature, a_vmr_list
            // Agenda output: abs, nlte
            propmat_clearsky_agendaExecute(wss,
                                           abs,
                                           nlte,
                                           partial_abs,
                                           partial_nlte,
                                           jacobian_quantities,
                                           {},
                                           this_f_grid,
                                           los,
                                           AtmPoint{},  // FIXME: DUMMY VALUE
                                           abs_agenda);
            
            // Convert from derivative to absorption
            for (Index ispec=0; ispec<partial_abs.nelem(); ispec++) {
              partial_abs[ispec] *= a_vmr_list[ispec];
            }

            // Verify, that the number of elements in abs matrix is
            // constistent with stokes_dim:
            ARTS_USER_ERROR_IF (stokes_dim != abs.StokesDimensions(),
                "propmat_clearsky_fieldCalc was called with stokes_dim = ",
                stokes_dim, ",\n"
                "but the stokes_dim returned by the agenda is ",
                abs.StokesDimensions(), ".")

            // Verify, that the number of species in abs is
            // constistent with vmr_field:
            ARTS_USER_ERROR_IF (n_species != partial_abs.nelem(),
                "The number of gas species in vmr_field is ", n_species,
                ",\n"
                "but the number of species returned by the agenda is ",
                partial_abs.nelem(), ".")

            // Verify, that the number of frequencies in abs is
            // constistent with f_extent:
            ARTS_USER_ERROR_IF (n_frequencies != abs.NumberOfFrequencies(),
                "The number of frequencies desired is ", n_frequencies,
                ",\n"
                "but the number of frequencies returned by the agenda is ",
                abs.NumberOfFrequencies(), ".")

            // Store the result in output field.
            // We have to transpose abs, because the dimensions of the
            // two variables are:
            // abs_field: [abs_species, f_grid, stokes, stokes, p_grid, lat_grid, lon_grid]
            // abs:       [abs_species][f_grid, stokes, stokes]
            for (Index i = 0; i < partial_abs.nelem(); i++) {
              partial_abs[i].GetTensor3(propmat_clearsky_field(
                  i, joker, joker, joker, ipr, ila, ilo));
            }
          }
      } catch (const std::runtime_error& e) {
#pragma omp critical(propmat_clearsky_fieldCalc_fail)
        {
          fail_msg = e.what();
          failed = true;
        }
      }
    }

  if (failed) throw runtime_error(fail_msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromGasAbsLookup(Vector& f_grid,
                            const GasAbsLookup& abs_lookup,
                            const Verbosity&) {
  const Vector& lookup_f_grid = abs_lookup.GetFgrid();
  f_grid.resize(lookup_f_grid.nelem());
  f_grid = lookup_f_grid;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromGasAbsLookup(Vector& p_grid,
                            const GasAbsLookup& abs_lookup,
                            const Verbosity&) {
  const Vector& lookup_p_grid = abs_lookup.GetPgrid();
  p_grid.resize(lookup_p_grid.nelem());
  p_grid = lookup_p_grid;
}
