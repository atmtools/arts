/**
  @file   m_cloudradar.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2010-10-31

  @brief  Workspace functions related to simulation of cloud radars.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>


#include "arts_constants.h"
#include "arts_omp.h"
#include "atm.h"
#include <workspace.h>
#include "debug.h"
#include "logic.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "matpack_view.h"
#include "montecarlo.h"
#include <rtepack.h>
#include "rte.h"
#include "sensor.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric LOG10_EULER_NUMBER=Constant::log10_euler;

/* Workspace method: Doxygen documentation will be auto-generated */
void RadarOnionPeelingTableCalc(
    const Workspace& ws,
    ArrayOfGriddedField3& invtable,
    const Vector& f_grid,
    const ArrayOfString& scat_species,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const ArrayOfAgenda& pnd_agenda_array,
    const ArrayOfArrayOfString& pnd_agenda_array_input_names,
    const Index& i_species,
    const Vector& dbze_grid,
    const Vector& t_grid,
    const Numeric& wc_min,
    const Numeric& wc_max,
    const Numeric& ze_tref,
    const Numeric& k2)
{
  // Some index and sizes
  const Index nss = scat_data.size();
  const Index ndb = dbze_grid.nelem();
  const Index nt = t_grid.nelem();
  const Index iss = i_species;

  // Check input
  ARTS_USER_ERROR_IF (scat_data[0][0].T_grid[0] < 220,
                      "First element in T_grid of scattering species is "
                      "below 220 K.\nThat seems to be too low for liquid.\n"
                      "Please note that the onion peeling function expects "
                      "rain\nas the first scattering element.");
  ARTS_USER_ERROR_IF (f_grid.nelem() != 1,
                      "This method requires that *f_grid* has length 1.");
  ARTS_USER_ERROR_IF (i_species < 0 || i_species > 1,
                      "*i_species* must either be 0 or 1.");
  ARTS_USER_ERROR_IF (nss != 2,
                      "*scat_data* must contain data for exactly two "
                      "scattering species.");
  ARTS_USER_ERROR_IF (scat_species.size() != nss,
        "*scat_data* and *scat_species* are inconsistent in size.");
  ARTS_USER_ERROR_IF (scat_meta.size() != nss,
        "*scat_data* and *scat_meta* are inconsistent in size.");
  ARTS_USER_ERROR_IF (scat_data[iss].size() != scat_meta[iss].size(),
                      "*scat_data* and *scat_meta* have inconsistent sizes.");
  ARTS_USER_ERROR_IF (scat_data[iss][0].f_grid.nelem() != 1,
                      "*scat_data* should just contain one frequency.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names[iss].size() != 1,
                      "The PSD applied must be of 1-moment type.");
  
  // Allocate
  if (invtable.empty())
    invtable.resize(2);
  invtable[iss].set_name("Radar inversion table created by *RadarOnionPeelingTableCalc*");
  invtable[iss].resize(2, ndb, nt);
  invtable[iss].data = 0;
  invtable[iss].set_grid_name(0, "Radiative properties");
  invtable[iss].set_grid(0, ArrayOfString{"Log of water content","Extinction"});
  invtable[iss].set_grid_name(1, "Radar reflectivity");
  invtable[iss].set_grid(1, dbze_grid);
  invtable[iss].set_grid_name(2, "Temperature");
  invtable[iss].set_grid(2, t_grid);  

  // Determine back-scattering and extinction on t_grid, at each size
  const Index nse = scat_data[iss].size();
  Matrix b(nt,nse), e(nt,nse);
  {
    Matrix itw(nt,2);
    ArrayOfGridPos gp(nt);
    for (Index i=0; i<nse; i++) {
      ARTS_USER_ERROR_IF (scat_data[iss][i].ptype != PTYPE_TOTAL_RND,
                          "So far only TRO is handled by this method.");
      const Index ib = scat_data[iss][i].za_grid.nelem() - 1;
      ARTS_USER_ERROR_IF (scat_data[iss][i].za_grid[ib] < 179.999,
          "All za_grid in scat_data must end with 180.\n"
          "This is not the case for scat_data[", iss, "][",
          i, "] (0-based)")
      gridpos(gp, scat_data[iss][i].T_grid, t_grid, 1);
      interpweights(itw, gp);
      interp(e(joker,i), itw, scat_data[iss][i].ext_mat_data(0,joker,0,0,0), gp);
      interp(b(joker,i), itw, scat_data[iss][i].pha_mat_data(0,joker,ib,0,0,0,0), gp);
    }
  }

  // Create test grid for water content
  Vector wc_grid;
  VectorLogSpace(wc_grid, wc_min, wc_max, 0.04);
  const Index nwc = wc_grid.nelem();
  
  // Calculate dBZe and extinction for wc_grid
  //
  Tensor3 D(2, nwc, nt, 0);
  //
  const Vector& pnd_agenda_input_t = t_grid;
  Matrix pnd_agenda_input(nt, 1);
  ArrayOfString dpnd_data_dx_names(0);
  //
  Vector cfac(1);
  ze_cfac(cfac, f_grid, ze_tref, k2);
  //
  for (Index w=0; w<nwc; w++) {
    // Get pnd
    pnd_agenda_input = wc_grid[w];
    Matrix pnd_data;
    Tensor3 dpnd_data_dx;    
    pnd_agenda_arrayExecute(ws,
                            pnd_data,
                            dpnd_data_dx,
                            iss,
                            pnd_agenda_input_t,
                            pnd_agenda_input,
                            pnd_agenda_array_input_names[iss],
                            dpnd_data_dx_names,
                            pnd_agenda_array);
    
    // Sum up to get bulk back-scattering and extinction
    for (Index t=0; t<nt; t++) {
      for (Index i=0; i<nse; i++) {
        ARTS_USER_ERROR_IF (b(t,i) < 0,
          "A negative back-scattering found for scat_species ", iss,
          ",\ntemperature ", t_grid[t], "K and element ", i)
        ARTS_USER_ERROR_IF (e(t,i) < 0,
          "A negative extinction found for scat_species ", iss, 
          ",\ntemperature ", t_grid[t], "K and element ", i)
        ARTS_USER_ERROR_IF (pnd_data(t,i) < 0,
          "A negative PSD value found for scat_species ", iss, 
          ",\ntemperature ", t_grid[t], "K and ", wc_grid[w],
          " kg/m3")
        D(0,w,t) += pnd_data(t,i) * b(t,i);
        D(1,w,t) += pnd_data(t,i) * e(t,i);
      }
      // Convert to dBZe
      D(0,w,t) = 10 * log10(cfac[0] * D(0,w,t));      
    }
  }

  // Get water content and extinction as a function of dBZe by interpolation
  Matrix itw(ndb,2);
  ArrayOfGridPos gp(ndb);
  // Water content interpolated in log
  Vector wc_log(nwc);
  transform(wc_log, log10, wc_grid);
  for (Index t=0; t<nt; t++) {
    if (!is_increasing(D(0,joker,t))) {
      for (Index w=0; w<nwc; w++) {
        std::cout << wc_grid[w] << " " << D(0,w,t) << std::endl;
      }
      ARTS_USER_ERROR (
        "A case found of non-increasing dBZe.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")
    }
    if (D(0,0,t) > dbze_grid[0]) {
      for (Index w=0; w<nwc; w++) {
        std::cout << wc_grid[w] << " " << D(0,w,t) << std::endl;
      }
      ARTS_USER_ERROR (
        "A case found where start of dbze_grid not covered.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")      
    }
    if (D(0,nwc-1,t) < dbze_grid[ndb-1]) {
      for (Index w=0; w<nwc; w++) {
        std::cout << wc_grid[w] << " " << D(0,w,t) << std::endl;
      }
      ARTS_USER_ERROR (
        "A case found where end of dbze_grid not covered.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")
    }
    //
    gridpos(gp, D(0,joker,t), dbze_grid);
    interpweights(itw, gp);
    interp(invtable[iss].data(0,joker,t), itw, wc_log, gp);
    interp(invtable[iss].data(1,joker,t), itw, D(1,joker,t), gp);
  }
}
