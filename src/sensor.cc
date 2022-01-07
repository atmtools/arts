/* Copyright (C) 2003-2012 Mattias Ekstr�m <ekstrom@rss.chalmers.se>

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

/*!
  \file   sensor.cc
  \author Mattias Ekstr�m <ekstrom@rss.chalmers.se>
  \date   2003-02-27

  \brief  Functions related to sensor modelling.

  Functions to model sensor behaviour and integration calculated as vector
  multiplication.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "sensor.h"
#include <cmath>
#include <list>
#include <stdexcept>
#include "arts.h"
#include "logic.h"
#include "matpackI.h"
#include "matpackII.h"
#include "messages.h"
#include "rte.h"
#include "sensor.h"
#include "sorting.h"

extern const Numeric PI;
extern const Numeric NAT_LOG_2;
extern const Numeric DEG2RAD;
extern const Index GFIELD1_F_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_ZA_GRID;
extern const Index GFIELD4_AA_GRID;

/*===========================================================================
  === The functions, besides the core integration and sum functions that are
  === placed together at the end
  ===========================================================================*/

void antenna1d_matrix(Sparse& H,
#ifndef NDEBUG
                      const Index& antenna_dim,
#else
                      const Index& antenna_dim _U_,
#endif
                      ConstVectorView antenna_dza,
                      const GriddedField4& antenna_response,
                      ConstVectorView za_grid,
                      ConstVectorView f_grid,
                      const Index n_pol,
                      const Index do_norm) {
  // Number of input pencil beam directions and frequency
  const Index n_za = za_grid.nelem();
  const Index n_f = f_grid.nelem();

  // Calculate number of antenna beams
  const Index n_ant = antenna_dza.nelem();

  // Asserts for variables beside antenna_response
  ARTS_ASSERT(antenna_dim == 1);
  ARTS_ASSERT(n_za >= 2);
  ARTS_ASSERT(n_pol >= 1);
  ARTS_ASSERT(do_norm >= 0 && do_norm <= 1);

  // Extract antenna_response grids
  const Index n_ar_pol =
      antenna_response.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
  ConstVectorView aresponse_f_grid =
      antenna_response.get_numeric_grid(GFIELD4_F_GRID);
  ConstVectorView aresponse_za_grid =
      antenna_response.get_numeric_grid(GFIELD4_ZA_GRID);
  DEBUG_ONLY(const Index n_ar_aa =
                 antenna_response.get_numeric_grid(GFIELD4_AA_GRID).nelem();)

  //
  const Index n_ar_f = aresponse_f_grid.nelem();
  const Index n_ar_za = aresponse_za_grid.nelem();
  const Index pol_step = n_ar_pol > 1;

  // Asserts for antenna_response
  ARTS_ASSERT(n_ar_pol == 1 || n_ar_pol >= n_pol);
  ARTS_ASSERT(n_ar_f);
  ARTS_ASSERT(n_ar_za > 1);
  ARTS_ASSERT(n_ar_aa == 1);

  // If response data extend outside za_grid is checked in
  // integration_func_by_vecmult

  // Some size(s)
  const Index nfpol = n_f * n_pol;

  // Resize H
  H.resize(n_ant * nfpol, n_za * nfpol);

  // Storage vectors for response weights
  Vector hrow(H.ncols(), 0.0);
  Vector hza(n_za, 0.0);

  // Antenna response to apply (possibly obtained by frequency interpolation)
  Vector aresponse(n_ar_za, 0.0);

  // Antenna beam loop
  for (Index ia = 0; ia < n_ant; ia++) {
    Vector shifted_aresponse_za_grid = aresponse_za_grid;
    shifted_aresponse_za_grid += antenna_dza[ia];

    // Order of loops assumes that the antenna response more often
    // changes with frequency than for polarisation

    // Frequency loop
    for (Index f = 0; f < n_f; f++) {
      // Polarisation loop
      for (Index ip = 0; ip < n_pol; ip++) {
        // Determine antenna pattern to apply
        //
        // Interpolation needed only if response has a frequency grid
        //
        Index new_antenna = 1;
        //
        if (n_ar_f == 1)  // No frequency variation
        {
          if (pol_step)  // Polarisation variation, update always needed
          {
            aresponse = antenna_response.data(ip, 0, joker, 0);
          } else if (f == 0 && ip == 0)  // Set fully constant pattern
          {
            aresponse = antenna_response.data(0, 0, joker, 0);
          } else  // The one set just above can be reused
          {
            new_antenna = 0;
          }
        } else {
          if (ip == 0 || pol_step)  // Interpolation required
          {
            // Interpolation (do this in "green way")
            ArrayOfGridPos gp_f(1), gp_za(n_ar_za);
            gridpos(gp_f, aresponse_f_grid, Vector(1, f_grid[f]));
            gridpos(gp_za, aresponse_za_grid, aresponse_za_grid);
            Tensor3 itw(1, n_ar_za, 4);
            interpweights(itw, gp_f, gp_za);
            Matrix aresponse_matrix(1, n_ar_za);
            interp(aresponse_matrix,
                   itw,
                   antenna_response.data(ip, joker, joker, 0),
                   gp_f,
                   gp_za);
            aresponse = aresponse_matrix(0, joker);
          } else  // Reuse pattern for ip==0
          {
            new_antenna = 0;
          }
        }

        // Calculate response weights
        if (new_antenna) {
          integration_func_by_vecmult(
              hza, aresponse, shifted_aresponse_za_grid, za_grid);

          // Normalisation?
          if (do_norm) {
            hza /= hza.sum();
          }
        }

        // Put weights into H
        //
        const Index ii = f * n_pol + ip;
        //
        hrow[Range(ii, n_za, nfpol)] = hza;
        //
        H.insert_row(ia * nfpol + ii, hrow);
        //
        hrow = 0;
      }
    }
  }
}



void antenna2d_gridded_dlos(Sparse& H,
#ifndef NDEBUG
                            const Index& antenna_dim,
#else
                            const Index& antenna_dim _U_,
#endif
                            ConstMatrixView antenna_dlos,
                            const GriddedField4& antenna_response,
                            ConstMatrixView mblock_dlos,
                            ConstVectorView f_grid,
                            const Index n_pol)
{
  // Number of input pencil beam directions and frequency
  const Index n_dlos = mblock_dlos.nrows();
  const Index n_f = f_grid.nelem();

  // Decompose mblock_dlos into za and aa grids, including checking
  ARTS_USER_ERROR_IF (mblock_dlos.ncols() != 2 ,
                      "For the gridded_dlos option, *mblock_dlos_grid* "
                      "must have two columns.");


  Index nza = 1;
  for(Index i=0; i<n_dlos-1 && mblock_dlos(i+1,0) > mblock_dlos(i,0); i++ ) {
    nza++;
  }
  ARTS_USER_ERROR_IF(nza < 2,
                     "For the gridded_dlos option, the number of za angles "
                     "(among dlos directions) must be >= 2.");
  ARTS_USER_ERROR_IF (n_dlos % nza > 0,
                      "For the gridded_dlos option, the number of dlos angles "
                      "must be a product of two integers.");
  const Index naa = n_dlos / nza; 
  const Vector za_grid = mblock_dlos(Range(0,nza),0);
  const Vector aa_grid = mblock_dlos(Range(0,naa,nza),1);
  for(Index i=0; i<n_dlos; i++ ) {
    ARTS_USER_ERROR_IF(i>=nza && abs(mblock_dlos(i,0)-mblock_dlos(i-nza,0)) > 1e-6 ,
        "Zenith angle of dlos ", i, " (0-based) differs to zenith " 
        "angle of dlos ", i-nza, ", while they need to be equal "
        "to form rectangular grid.")
    ARTS_USER_ERROR_IF(abs(mblock_dlos(i,1)-aa_grid[i/nza]) > 1e-6,
        "Azimuth angle of dlos ", i, " (0-based) differs to azimuth " 
        "angle ", (i/nza)*nza , ", while they need to be equal "
        "to form rectangular grid.")
  }

  // Calculate number of antenna beams
  const Index n_ant = antenna_dlos.nrows();

  // Asserts for variables beside antenna_response
  ARTS_ASSERT(antenna_dim == 2);
  ARTS_ASSERT(n_dlos >= 1);
  ARTS_ASSERT(n_pol >= 1);

  // Extract antenna_response grids
  const Index n_ar_pol =
      antenna_response.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
  ConstVectorView aresponse_f_grid =
      antenna_response.get_numeric_grid(GFIELD4_F_GRID);
  ConstVectorView aresponse_za_grid =
      antenna_response.get_numeric_grid(GFIELD4_ZA_GRID);
  ConstVectorView aresponse_aa_grid =
      antenna_response.get_numeric_grid(GFIELD4_AA_GRID);
  //
  const Index n_ar_f = aresponse_f_grid.nelem();
  const Index n_ar_za = aresponse_za_grid.nelem();
  const Index n_ar_aa = aresponse_aa_grid.nelem();
  const Index pol_step = n_ar_pol > 1;

  // Asserts for antenna_response
  ARTS_ASSERT(n_ar_pol == 1 || n_ar_pol >= n_pol);
  ARTS_ASSERT(n_ar_f);
  ARTS_ASSERT(n_ar_za > 1);
  ARTS_ASSERT(n_ar_aa > 1);
  ARTS_ASSERT(antenna_response.data.ncols() == n_ar_aa );
  ARTS_ASSERT(antenna_response.data.nrows() == n_ar_za );
  ARTS_ASSERT(antenna_response.data.npages() == n_ar_f );
  ARTS_ASSERT(antenna_response.data.nbooks() == n_ar_pol );

  // Include cos(za)-term in response
  Tensor4 aresponse_with_cos(n_ar_pol, n_ar_f, n_ar_za, n_ar_aa);
  for(Index i3=0; i3<n_ar_za; i3++) {
    const Numeric t = cos(DEG2RAD * aresponse_za_grid[i3]);
    for(Index i4=0; i4<n_ar_aa; i4++) {
      for(Index i2=0; i2<n_ar_f; i2++) {
        for(Index i1=0; i1<n_ar_pol; i1++) {
          aresponse_with_cos(i1,i2,i3,i4) = t * antenna_response.data(i1,i2,i3,i4);
        }
      }
    }
  }

  // Some size(s)
  const Index nfpol = n_f * n_pol;

  // Resize H
  H.resize(n_ant * nfpol, n_dlos * nfpol);

  // Storage vectors for response weights
  Vector hrow(H.ncols(), 0.0);
  Vector hza(n_dlos, 0.0);

  // Antenna response to apply (possibly obtained by frequency interpolation)
  Matrix aresponse(n_ar_za, n_ar_aa, 0.0);

  // If you find a bug or change something, likely also change other 2D antenna
  // function(s) as they are similar
  
  // Antenna beam loop
  for (Index ia = 0; ia < n_ant; ia++) {
    // Order of loops assumes that the antenna response more often
    // changes with frequency than for polarisation

    // Frequency loop
    for (Index f = 0; f < n_f; f++) {
      // Polarisation loop
      for (Index ip = 0; ip < n_pol; ip++) {
        // Determine antenna pattern to apply
        //
        // Interpolation needed only if response has a frequency grid
        //
        Index new_antenna = 1;
        //
        if (n_ar_f == 1)  // No frequency variation
        {
          if (pol_step)  // Polarisation variation, update always needed
          {
            aresponse = aresponse_with_cos(ip, 0, joker, joker);
          } else if (f == 0 && ip == 0)  // Set fully constant pattern
          {
            aresponse = aresponse_with_cos(0, 0, joker, joker);
          } else  // The one set just above can be reused
          {
            new_antenna = 0;
          }
        } else {
          if (ip == 0 || pol_step) {
            // Interpolation (do this in "green way")
            ArrayOfGridPos gp_f(1), gp_za(n_ar_za), gp_aa(n_ar_aa);
            gridpos(gp_f, aresponse_f_grid, Vector(1, f_grid[f]));
            gridpos(gp_za, aresponse_za_grid, aresponse_za_grid);
            gridpos(gp_aa, aresponse_aa_grid, aresponse_aa_grid);
            Tensor4 itw(1, n_ar_za, n_ar_aa, 8);
            interpweights(itw, gp_f, gp_za, gp_aa);
            Tensor3 aresponse_matrix(1, n_ar_za, n_ar_aa);
            interp(aresponse_matrix,
                   itw,
                   aresponse_with_cos(ip, joker, joker, joker),
                   gp_f,
                   gp_za,
                   gp_aa);
            aresponse = aresponse_matrix(0, joker, joker);
          } else  // Reuse pattern for ip==0
          {
            new_antenna = 0;
          }
        }

        // Calculate response weights, by using grid positions and "itw"
        if (new_antenna) {
          
          // za grid positions
          Vector zas = aresponse_za_grid;
          zas += antenna_dlos(ia, 0);
          ARTS_USER_ERROR_IF( zas[0] < za_grid[0],
              "The zenith angle grid in *mblock_dlos_grid* is too narrow. " 
              "It must be extended downwards with at least ",
              za_grid[0]-zas[0], " deg.")
          ARTS_USER_ERROR_IF( zas[n_ar_za-1] > za_grid[nza-1],
              "The zenith angle grid in *mblock_dlos_grid* is too narrow. " 
              "It must be extended upwards with at least ",
              zas[n_ar_za-1] - za_grid[nza-1], " deg.")
          
          ArrayOfGridPos gp_za(n_ar_za);
          gridpos(gp_za, za_grid, zas);
          
          // aa grid positions
          Vector aas = aresponse_aa_grid;
          if (antenna_dlos.ncols() > 1) { aas += antenna_dlos(ia, 1); }              
          ARTS_USER_ERROR_IF( aas[0] < aa_grid[0],
              "The azimuth angle grid in *mblock_dlos_grid* is too narrow. " 
              "It must be extended downwards with at least ",
              aa_grid[0]-aas[0], " deg.")
          ARTS_USER_ERROR_IF( aas[n_ar_aa-1] > aa_grid[naa-1],
              "The azimuth angle grid in *mblock_dlos_grid* is too narrow. " 
              "It must be extended upwards with at least ",
              aas[n_ar_aa-1] - aa_grid[naa-1], " deg.")
          
          ArrayOfGridPos gp_aa(n_ar_aa);
          gridpos(gp_aa, aa_grid, aas);


          // Derive interpolation weights
          Tensor3 itw(n_ar_za, n_ar_za, 4);
          interpweights(itw, gp_za, gp_aa);

          // Convert iwt to weights for H
          //
          hza = 0;   // Note that values in hza must be accumulated 
          //
          for (Index iaa = 0; iaa < n_ar_aa; iaa++) {
            const Index a = gp_aa[iaa].idx;
            const Index b = a + 1;
            
            for (Index iza = 0; iza < n_ar_za; iza++) {          

              const Index z = gp_za[iza].idx;
              const Index x = z + 1;

              if( itw(iza,iaa,0) > 1e-9 ) {
                hza[a*nza+z] += aresponse(iza,iaa) * itw(iza,iaa,0);
              }
              if( itw(iza,iaa,1) > 1e-9 ) {
                hza[b*nza+z] += aresponse(iza,iaa) * itw(iza,iaa,1);
              }
              if( itw(iza,iaa,2) > 1e-9 ) {
                hza[a*nza+x] += aresponse(iza,iaa) * itw(iza,iaa,2);
              }
              if( itw(iza,iaa,3) > 1e-9 ) {
                hza[b*nza+x] += aresponse(iza,iaa) * itw(iza,iaa,3);
              }
            }
          }

          // For 2D antennas we always normalise
          hza /= hza.sum();
        }

        // Put weights into H
        //
        const Index ii = f * n_pol + ip;
        //
        hrow[Range(ii, n_dlos, nfpol)] = hza;
        //
        H.insert_row(ia * nfpol + ii, hrow);
        //
        hrow = 0;
      }
    }
  }  
}



void antenna2d_interp_response(Sparse& H,
#ifndef NDEBUG
                               const Index& antenna_dim,
#else
                               const Index& antenna_dim _U_,
#endif
                               ConstMatrixView antenna_dlos,
                               const GriddedField4& antenna_response,
                               ConstMatrixView mblock_dlos,
                               ConstVectorView f_grid,
                               const Index n_pol)
{  
  // Number of input pencil beam directions and frequency
  const Index n_dlos = mblock_dlos.nrows();
  const Index n_f = f_grid.nelem();

  // Calculate number of antenna beams
  const Index n_ant = antenna_dlos.nrows();

  // Asserts for variables beside antenna_response
  ARTS_ASSERT(antenna_dim == 2);
  ARTS_ASSERT(n_dlos >= 1);
  ARTS_ASSERT(n_pol >= 1);

  // Extract antenna_response grids
  const Index n_ar_pol =
      antenna_response.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
  ConstVectorView aresponse_f_grid =
      antenna_response.get_numeric_grid(GFIELD4_F_GRID);
  ConstVectorView aresponse_za_grid =
      antenna_response.get_numeric_grid(GFIELD4_ZA_GRID);
  ConstVectorView aresponse_aa_grid =
      antenna_response.get_numeric_grid(GFIELD4_AA_GRID);
  //
  const Index n_ar_f = aresponse_f_grid.nelem();
  const Index n_ar_za = aresponse_za_grid.nelem();
  const Index n_ar_aa = aresponse_aa_grid.nelem();
  const Index pol_step = n_ar_pol > 1;

  // Asserts for antenna_response
  ARTS_ASSERT(n_ar_pol == 1 || n_ar_pol >= n_pol);
  ARTS_ASSERT(n_ar_f);
  ARTS_ASSERT(n_ar_za > 1);
  ARTS_ASSERT(n_ar_aa > 1);
  ARTS_ASSERT(antenna_response.data.ncols() == n_ar_aa );
  ARTS_ASSERT(antenna_response.data.nrows() == n_ar_za );
  ARTS_ASSERT(antenna_response.data.npages() == n_ar_f );
  ARTS_ASSERT(antenna_response.data.nbooks() == n_ar_pol );

  // Include cos(za)-term in response
  Tensor4 aresponse_with_cos(n_ar_pol, n_ar_f, n_ar_za, n_ar_aa);
  for(Index i3=0; i3<n_ar_za; i3++) {
    const Numeric t = cos(DEG2RAD * aresponse_za_grid[i3]);
    for(Index i4=0; i4<n_ar_aa; i4++) {
      for(Index i2=0; i2<n_ar_f; i2++) {
        for(Index i1=0; i1<n_ar_pol; i1++) {
          aresponse_with_cos(i1,i2,i3,i4) = t * antenna_response.data(i1,i2,i3,i4);
        }
      }
    }
  }
  
  // Some size(s)
  const Index nfpol = n_f * n_pol;

  // Resize H
  H.resize(n_ant * nfpol, n_dlos * nfpol);

  // Storage vectors for response weights
  Vector hrow(H.ncols(), 0.0);
  Vector hza(n_dlos, 0.0);

  // Antenna response to apply (possibly obtained by frequency interpolation)
  Matrix aresponse(n_ar_za, n_ar_aa, 0.0);

  // If you find a bug or change something, likely also change other 2D antenna
  // function(s) as they are similar

  // Antenna beam loop
  for (Index ia = 0; ia < n_ant; ia++) {
    // Order of loops assumes that the antenna response more often
    // changes with frequency than for polarisation

    // Frequency loop
    for (Index f = 0; f < n_f; f++) {
      // Polarisation loop
      for (Index ip = 0; ip < n_pol; ip++) {
        // Determine antenna pattern to apply
        //
        // Interpolation needed only if response has a frequency grid
        //
        Index new_antenna = 1;
        //
        if (n_ar_f == 1)  // No frequency variation
        {
          if (pol_step)  // Polarisation variation, update always needed
          {
            aresponse = aresponse_with_cos(ip, 0, joker, joker);
          } else if (f == 0 && ip == 0)  // Set fully constant pattern
          {
            aresponse = aresponse_with_cos(0, 0, joker, joker);
          } else  // The one set just above can be reused
          {
            new_antenna = 0;
          }
        } else {
          if (ip == 0 || pol_step) {
            // Interpolation (do this in "green way")
            ArrayOfGridPos gp_f(1), gp_za(n_ar_za), gp_aa(n_ar_aa);
            gridpos(gp_f, aresponse_f_grid, Vector(1, f_grid[f]));
            gridpos(gp_za, aresponse_za_grid, aresponse_za_grid);
            gridpos(gp_aa, aresponse_aa_grid, aresponse_aa_grid);
            Tensor4 itw(1, n_ar_za, n_ar_aa, 8);
            interpweights(itw, gp_f, gp_za, gp_aa);
            Tensor3 aresponse_matrix(1, n_ar_za, n_ar_aa);
            interp(aresponse_matrix,
                   itw,
                   aresponse_with_cos(ip, joker, joker, joker),
                   gp_f,
                   gp_za,
                   gp_aa);
            aresponse = aresponse_matrix(0, joker, joker);
          } else  // Reuse pattern for ip==0
          {
            new_antenna = 0;
          }
        }

        // Calculate response weights
        if (new_antenna) {
          for (Index l = 0; l < n_dlos; l++) {
            const Numeric za = mblock_dlos(l, 0) - antenna_dlos(ia, 0);
            Numeric aa = 0.0;
            if (mblock_dlos.ncols() > 1) {
              aa += mblock_dlos(l, 1);
            }
            if (antenna_dlos.ncols() > 1) {
              aa -= antenna_dlos(ia, 1);
            }

            // The response is zero if mblock_dlos is outside of
            // antennna pattern
            if (za < aresponse_za_grid[0] ||
                za > aresponse_za_grid[n_ar_za - 1] ||
                aa < aresponse_aa_grid[0] ||
                aa > aresponse_aa_grid[n_ar_aa - 1]) {
              hza[l] = 0;
            }
            // Otherwise we make an (blue) interpolation
            else {
              ArrayOfGridPos gp_za(1), gp_aa(1);
              gridpos(gp_za, aresponse_za_grid, Vector(1, za));
              gridpos(gp_aa, aresponse_aa_grid, Vector(1, aa));
              Matrix itw(1, 4);
              interpweights(itw, gp_za, gp_aa);
              Vector value(1);
              interp(value, itw, aresponse, gp_za, gp_aa);
              hza[l] = value[0];
            }
          }

          // For 2D antennas we always normalise
          hza /= hza.sum();
        }

        // Put weights into H
        //
        const Index ii = f * n_pol + ip;
        //
        hrow[Range(ii, n_dlos, nfpol)] = hza;
        //
        H.insert_row(ia * nfpol + ii, hrow);
        //
        hrow = 0;
      }
    }
  }
}



void gaussian_response_autogrid(Vector& x,
                                Vector& y,
                                const Numeric& x0,
                                const Numeric& fwhm,
                                const Numeric& xwidth_si,
                                const Numeric& dx_si) {
  ARTS_ASSERT(dx_si <= xwidth_si);

  const Numeric si = fwhm / (2 * sqrt(2 * NAT_LOG_2));

  // Number of points needed to enure that spacing is max dx_si
  const Index n = (Index)floor(2 * xwidth_si / dx_si) + 1;

  // Grid for response
  const Numeric dd = si * xwidth_si;
  nlinspace(x, -dd, dd, n);

  // Calculate response
  gaussian_response(y, x, x0, fwhm);
}



void gaussian_response(Vector& y,
                       const Vector& x,
                       const Numeric& x0,
                       const Numeric& fwhm) {
  const Numeric si = fwhm / (2 * sqrt(2 * NAT_LOG_2));
  const Numeric a = 1 / (si * sqrt(2 * PI));
  const Index n = x.nelem();

  y.resize(n);
  //
  for (Index i = 0; i < n; i++)
    y[i] = a * exp(-0.5 * pow((x[i] - x0) / si, 2.0));
}



void mixer_matrix(Sparse& H,
                  Vector& f_mixer,
                  const Numeric& lo,
                  const GriddedField1& filter,
                  ConstVectorView f_grid,
                  const Index& n_pol,
                  const Index& n_sp,
                  const Index& do_norm) {
  // Frequency grid of for sideband response specification
  ConstVectorView filter_grid = filter.get_numeric_grid(GFIELD1_F_GRID);

  DEBUG_ONLY(const Index nrp = filter.data.nelem();)

  // Asserts
  ARTS_ASSERT(lo > f_grid[0]);
  ARTS_ASSERT(lo < last(f_grid));
  ARTS_ASSERT(filter_grid.nelem() == nrp);
  ARTS_ASSERT(fabs(last(filter_grid) + filter_grid[0]) < 1e3);
  // If response data extend outside f_grid is checked in summation_by_vecmult

  // Find indices in f_grid where f_grid is just below and above the
  // lo frequency.
  /*
  Index i_low = 0, i_high = f_grid.nelem()-1, i_mean;
  while( i_high-i_low > 1 )
    {
      i_mean = (Index) (i_high+i_low)/2;
      if (f_grid[i_mean]<lo)
        { 
          i_low = i_mean; 
        }
      else
        {
          i_high = i_mean;
        }
    }
  if (f_grid[i_high]==lo)
    {
      i_high++;
    }
  const Numeric lim_low  = max( lo-f_grid[i_low], f_grid[i_high]-lo );
  */

  // Determine IF limits for new frequency grid
  const Numeric lim_low = 0;
  const Numeric lim_high = -filter_grid[0];

  // Convert sidebands to IF and use list to make a unique sorted
  // vector, this sorted vector is f_mixer.
  list<Numeric> l_mixer;
  for (Index i = 0; i < f_grid.nelem(); i++) {
    if (fabs(f_grid[i] - lo) >= lim_low && fabs(f_grid[i] - lo) <= lim_high) {
      l_mixer.push_back(fabs(f_grid[i] - lo));
    }
  }
  l_mixer.push_back(lim_high);  // Not necessarily a point in f_grid
  l_mixer.sort();
  l_mixer.unique();
  f_mixer.resize((Index)l_mixer.size());
  Index e = 0;
  for (list<Numeric>::iterator li = l_mixer.begin(); li != l_mixer.end();
       li++) {
    f_mixer[e] = *li;
    e++;
  }

  // Resize H
  H.resize(f_mixer.nelem() * n_pol * n_sp, f_grid.nelem() * n_pol * n_sp);

  // Calculate the sensor summation vector and insert the values in the
  // final matrix taking number of polarisations and zenith angles into
  // account.
  Vector row_temp(f_grid.nelem());
  Vector row_final(f_grid.nelem() * n_pol * n_sp);
  //
  Vector if_grid = f_grid;
  if_grid -= lo;
  //
  for (Index i = 0; i < f_mixer.nelem(); i++) {
    summation_by_vecmult(
        row_temp, filter.data, filter_grid, if_grid, f_mixer[i], -f_mixer[i]);

    // Normalise if flag is set
    if (do_norm) row_temp /= row_temp.sum();

    // Loop over number of polarisations
    for (Index p = 0; p < n_pol; p++) {
      // Loop over number of zenith angles/antennas
      for (Index a = 0; a < n_sp; a++) {
        // Distribute elements of row_temp to row_final.
        row_final = 0.0;
        row_final[Range(
            a * f_grid.nelem() * n_pol + p, f_grid.nelem(), n_pol)] = row_temp;
        H.insert_row(a * f_mixer.nelem() * n_pol + p + i * n_pol, row_final);
      }
    }
  }
}



void met_mm_polarisation_hmatrix(Sparse& H,
                                 const ArrayOfString& mm_pol,
                                 const Numeric dza,
                                 const Index stokes_dim,
                                 const String& iy_unit) {
  ARTS_ASSERT(stokes_dim > 1);

  // Set "Stokes vector weights" according to iy_unit
  Numeric w = 0.5;
  if (iy_unit == "PlanckBT" || iy_unit == "RJBT") {
    w = 1.0;
  }

  // Identify (basic) polarisation response and possible sensor specific
  // "rotation"
  const Index nch = mm_pol.nelem();  // Number of channels
  ArrayOfString pol(nch);
  ArrayOfString rot(nch);
  for (Index i = 0; i < nch; i++) {
    if (mm_pol[i] == "AMSU-H") {
      rot[i] = "AMSU";
      pol[i] = "H";
    } else if (mm_pol[i] == "AMSU-V") {
      rot[i] = "AMSU";
      pol[i] = "V";
    } else if (mm_pol[i] == "ISMAR-H") {
      rot[i] = "ISMAR";
      pol[i] = "H";
    } else if (mm_pol[i] == "ISMAR-V") {
      rot[i] = "ISMAR";
      pol[i] = "V";
    } else if (mm_pol[i] == "MARSS-H") {
      rot[i] = "MARSS";
      pol[i] = "H";
    } else if (mm_pol[i] == "MARSS-V") {
      rot[i] = "MARSS";
      pol[i] = "V";
    } else if (mm_pol[i] == "H" || mm_pol[i] == "V" || mm_pol[i] == "LHC" ||
               mm_pol[i] == "RHC") {
      rot[i] = "none";
      pol[i] = mm_pol[i];
    } else {
      ARTS_USER_ERROR ( "Unknown polarisation ", mm_pol[i])
    }
  }

  // Vectors representing standard cases of sensor polarisation response
  /*
    ArrayOfVector pv;
    stokes2pol( pv, w );
    */

  // Complete H, for all channels
  H = Sparse(nch, nch * stokes_dim);

  for (Index i = 0; i < nch; i++) {
    /*
        // See stokes2pol for index order used in pv
        Index ipv = -1;
        if( pol[i] == "V" )
        { ipv = 4; }
        else if( pol[i] == "H" )
        { ipv = 5; }
        else if( pol[i] == "LHC" )  // Left hand circular
        { ipv = 8; }
        else if( pol[i] == "RHC" )  // Right hand circular
        { ipv = 9; }
        else
        { ARTS_ASSERT( 0 ); }
      */
    // See instrument_pol for index order
    Index ipol = -1;
    if (pol[i] == "V") {
      ipol = 5;
    } else if (pol[i] == "H") {
      ipol = 6;
    } else if (pol[i] == "LHC")  // Left hand circular
    {
      ipol = 9;
    } else if (pol[i] == "RHC")  // Right hand circular
    {
      ipol = 10;
    } else {
      ARTS_ASSERT(0);
    }

    /*
        // Maybe this error messages should be mofified
        if( pv[ipv].nelem() > stokes_dim )
        {
            ostringstream os;
            os << "You have selected a channel polarisation that is not "
            << "covered by present value of *stokes_dim* (the later has to be "
            << "increased).";
            throw runtime_error(os.str());
            }*/

    // No rotation, just plane polarisation response
    if (rot[i] == "none") {
      // Here we just need to fill the row H
      Vector hrow(nch * stokes_dim, 0.0);
      /* Old code, matching older version of stokes2pol:
          hrow[Range(i*stokes_dim,pv[ipv].nelem())] = pv[ipv];
          */
      stokes2pol(hrow[Range(i * stokes_dim, stokes_dim)], stokes_dim, ipol, w);
      H.insert_row(i, hrow);
    }

    // Rotation + pol-response
    else {
      // Rotation part
      Sparse Hrot(stokes_dim, stokes_dim);
      if (rot[i] == "AMSU") {
        // No idea about the sign. Not important if U=0,
        // but matter for U != 0.
        muellersparse_rotation(Hrot, stokes_dim, abs(dza));
      } else if (rot[i] == "ISMAR") {
        // No rotation at -53 (= forward direction). But no idea about the
        // sign, as for AMSU above
        muellersparse_rotation(Hrot, stokes_dim, dza + 50);
      } else if (rot[i] == "MARSS") {
        // MARSS special, as 48 deg between different polarisation (not 90,
        // as for ISMAR. This is best interpretation of information
        // from Stuart. Should be double-checked with him at some point.
        if (pol[i] == "H") {
          muellersparse_rotation(Hrot, stokes_dim, dza + 42);
        } else {
          muellersparse_rotation(Hrot, stokes_dim, dza);
        }
      } else {
        ARTS_ASSERT(0);
      }

      // H-matrix matching polarization
      Sparse Hpol(1, stokes_dim);
      { /* Old code, matching older version of stokes2pol
            Vector hrow( stokes_dim, 0.0 );
            hrow[Range(0,pv[ipv].nelem())] = pv[ipv];
            */
        Vector hrow(stokes_dim);
        stokes2pol(hrow, stokes_dim, ipol, w);
        Hpol.insert_row(0, hrow);
      }

      // H for the individual channel
      Sparse Hc(1, stokes_dim);
      mult(Hc, Hpol, Hrot);

      // Put Hc into H
      Vector hrow(nch * stokes_dim, 0.0);
      const Index i0 = i * stokes_dim;
      for (Index s = 0; s < stokes_dim; s++) {
        hrow[i0 + s] = Hc(0, s);
      }
      H.insert_row(i, hrow);
    }
  }
}



void sensor_aux_vectors(Vector& sensor_response_f,
                        ArrayOfIndex& sensor_response_pol,
                        Matrix& sensor_response_dlos,
                        ConstVectorView sensor_response_f_grid,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        ConstMatrixView sensor_response_dlos_grid) {
  // Sizes
  const Index nf = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nlos = sensor_response_dlos_grid.nrows();
  const Index n = nf * npol * nlos;

  // Allocate
  sensor_response_f.resize(n);
  sensor_response_pol.resize(n);
  sensor_response_dlos.resize(n, sensor_response_dlos_grid.ncols());

  // Fill
  for (Index ilos = 0; ilos < nlos; ilos++) {
    const Index i2 = ilos * nf * npol;
    //
    for (Index ifr = 0; ifr < nf; ifr++) {
      const Index i3 = i2 + ifr * npol;
      //
      for (Index ip = 0; ip < npol; ip++) {
        const Index i = i3 + ip;
        //
        sensor_response_f[i] = sensor_response_f_grid[ifr];
        sensor_response_pol[i] = sensor_response_pol_grid[ip];
        sensor_response_dlos(i, joker) = sensor_response_dlos_grid(ilos, joker);
      }
    }
  }
}



void spectrometer_matrix(Sparse& H,
                         ConstVectorView ch_f,
                         const ArrayOfGriddedField1& ch_response,
                         ConstVectorView sensor_f,
                         const Index& n_pol,
                         const Index& n_sp,
                         const Index& do_norm) {
  // Check if matrix has one frequency column or one for every channel
  // frequency
  //
  ARTS_ASSERT(ch_response.nelem() == 1 || ch_response.nelem() == ch_f.nelem());
  //
  Index freq_full = ch_response.nelem() > 1;

  // If response data extend outside sensor_f is checked in
  // integration_func_by_vecmult

  // Reisze H
  //
  const Index nin_f = sensor_f.nelem();
  const Index nout_f = ch_f.nelem();
  const Index nin = n_sp * nin_f * n_pol;
  const Index nout = n_sp * nout_f * n_pol;
  //
  H.resize(nout, nin);

  // Calculate the sensor integration vector and put values in the temporary
  // vector, then copy vector to the transfer matrix
  //
  Vector ch_response_f;
  Vector weights(nin_f);
  Vector weights_long(nin, 0.0);
  //
  for (Index ifr = 0; ifr < nout_f; ifr++) {
    const Index irp = ifr * freq_full;

    //The spectrometer response is shifted for each centre frequency step
    ch_response_f = ch_response[irp].get_numeric_grid(GFIELD1_F_GRID);
    ch_response_f += ch_f[ifr];

    // Call *integration_func_by_vecmult* and store it in the temp vector
    integration_func_by_vecmult(
        weights, ch_response[irp].data, ch_response_f, sensor_f);

    // Normalise if flag is set
    if (do_norm) weights /= weights.sum();

    // Loop over polarisation and spectra (viewing directions)
    // Weights change only with frequency
    for (Index sp = 0; sp < n_sp; sp++) {
      for (Index pol = 0; pol < n_pol; pol++) {
        // Distribute the compact weight vector into weight_long
        weights_long[Range(sp * nin_f * n_pol + pol, nin_f, n_pol)] = weights;

        // Insert temp_long into H at the correct row
        H.insert_row(sp * nout_f * n_pol + ifr * n_pol + pol, weights_long);

        // Reset weight_long to zero.
        weights_long = 0.0;
      }
    }
  }
}



void stokes2pol(VectorView w,
                const Index& stokes_dim,
                const Index& ipol_1based,
                const Numeric nv) {
  ARTS_ASSERT(w.nelem() == stokes_dim);

  ARTS_USER_ERROR_IF (ipol_1based < 1 || ipol_1based > 10,
                      "Valid polarization indices are 1 to 10 (1-based).");

  ArrayOfVector s2p(10);
  //
  s2p[0] = {1};              // I
  s2p[1] = {0, 1};           // Q
  s2p[2] = {0, 0, 1};        // U
  s2p[3] = {0, 0, 0, 1};     // V
  s2p[4] = {nv, nv};         // Iv
  s2p[5] = {nv, -nv};        // Ih
  s2p[6] = {nv, 0, nv};      // I+45
  s2p[7] = {nv, 0, -nv};     // I-45
  s2p[8] = {nv, 0, 0, nv};   // Ilhc
  s2p[9] = {nv, 0, 0, -nv};  // Irhc

  const Index l = s2p[ipol_1based - 1].nelem();
  ARTS_USER_ERROR_IF (l > stokes_dim,
    "You have selected polarization with 1-based index: ", ipol_1based,
    "\n",
    "but this polarization demands stokes_dim >= ", l, "\n",
    "while the actual values of stokes_dim is ", stokes_dim)

  w[Range(0, l)] = s2p[ipol_1based - 1];
  if (l < stokes_dim) {
    w[Range(l, stokes_dim - l)] = 0;
  }
}

// Functions by Stefan, needed for HIRS:

//! Test if two instrument channels overlap, and if so, merge them.
/*!
  The channels boundaries are specified in two separate vectors, fmin
  and fmax. These vectors are both input and output. If merging has
  happened, they will each be one element shorter. 

  The positions of the channels to compare is given by the input
  parameters i and j. It is assumed that the minimum frequency of i
  is lower than or equal to that of j.

  Furthermore, it is assumed that i itself is lower than j.

  The range of the first channel (i) will have been extended to
  accomodate the second channel (j). The second channel will have been
  removed.

  The function also handles the updating of index j: If the two
  channels do not overlap, j is increased by one.

  Function returns true if merging has happened.

  \author Stefan Buehler
  
  \return True if channels were merged, otherwise false.
  \retval fmin Lower channel boundaries.
  \retval fmax Upper channel boundaries.
  \param i Index of first channel.
  \param j Index of second channel.
*/
bool test_and_merge_two_channels(Vector& fmin, Vector& fmax, Index i, Index j) {
  const Index nf = fmin.nelem();
  ARTS_ASSERT(fmax.nelem() == nf);
  ARTS_ASSERT(i >= 0 && i < nf);
  ARTS_ASSERT(j >= 0 && j < nf);
  ARTS_ASSERT(fmin[i] <= fmin[j]);
  ARTS_ASSERT(i < j);

  // There are three cases to consider:
  // a) The two channels are separate: fmax[i] <  fmin[j]
  // b) They overlap:                  fmax[i] >= fmin[j]
  // c) j is inside i:                 fmax[i] >  fmax[j]

  // In the easiest case (a), we do not have to do anything.
  if (fmax[i] >= fmin[j]) {
    // We are in case (b) or (c), so we know that we have to combine
    // the channels. The new minimum frequency is fmin[i]. The new
    // maximum frequency is the larger one of the two channels we
    // are combining:
    if (fmax[j] > fmax[i]) fmax[i] = fmax[j];

    // We now have to kick out element j from both vectors.

    // Number of elements behind j:
    Index n_behind = nf - 1 - j;

    Vector dummy = fmin;
    fmin.resize(nf - 1);
    fmin[Range(0, j)] = dummy[Range(0, j)];
    if (n_behind > 0) fmin[Range(j, n_behind)] = dummy[Range(j + 1, n_behind)];

    dummy = fmax;
    fmax.resize(nf - 1);
    fmax[Range(0, j)] = dummy[Range(0, j)];
    if (n_behind > 0) fmax[Range(j, n_behind)] = dummy[Range(j + 1, n_behind)];

    return true;
  }

  return false;
}

void find_effective_channel_boundaries(  // Output:
    Vector& fmin,
    Vector& fmax,
    // Input:
    const Vector& f_backend,
    const ArrayOfGriddedField1& backend_channel_response,
    const Numeric& delta,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  // How many channels in total:
  const Index n_chan = f_backend.nelem();

  // Checks on input quantities:

  // There must be at least one channel.
  ARTS_USER_ERROR_IF (n_chan < 1,
      "There must be at least one channel.\n"
      "(The vector *f_backend* must have at least one element.)")

  // There must be a response function for each channel.
  ARTS_USER_ERROR_IF (n_chan != backend_channel_response.nelem(),
      "Variables *f_backend_multi* and *backend_channel_response_multi*\n"
      "must have same number of bands for each LO.")

  // Frequency grids for response functions must be strictly increasing.
  for (Index i = 0; i < n_chan; ++i) {
    // Frequency grid for this response function:
    const Vector& backend_f_grid =
        backend_channel_response[i].get_numeric_grid(0);

    ARTS_USER_ERROR_IF (!is_increasing(backend_f_grid),
        "The frequency grid for the backend channel response of\n"
        "channel ", i, " is not strictly increasing.\n"
        "It is: ", backend_f_grid, "\n")
  }

  // Start the actual work.

  out2 << "  Original channel characteristics:\n"
       << "  min         nominal      max (all in Hz):\n";

  // count the number of passbands as defined by segments of filtershapes, backend_channel_response.data = [0,>0,>0,0]
  // Borders between passbands are identified as [...0,0...]

  // Get a list of original channel boundaries:
  Index numPB = 0;
  for (Index idx = 0; idx < n_chan; ++idx) {
    const Vector& backend_filter = backend_channel_response[idx].data;
    if (backend_filter.nelem() >
        2) {  // only run this code when there is more then two elements in the backend
      for (Index idy = 1; idy < backend_filter.nelem(); ++idy) {
        if ((backend_filter[idy] > 0) && (backend_filter[idy - 1] == 0)) {
          numPB++;  // Two consecutive zeros gives the border between two passbands
        }
      }
    } else {
      numPB++;
    }
  }

  ARTS_USER_ERROR_IF (!numPB,
        "No passbands found.\n"
        "*backend_channel_response* must be zero around the passbands.\n"
        "backend_channel_response.data = [0, >0, >0, 0]\n"
        "Borders between passbands are identified as [...0,0...]");

  Vector fmin_pb(numPB);
  Vector fmax_pb(numPB);
  Index pbIdx = 0;

  for (Index idx = 0; idx < n_chan; ++idx) {
    // Some handy shortcuts:
    //
    // We have to find the first and last frequency where the
    // response is actually different from 0. (No point in making
    // calculations for frequencies where the response is 0.)
    //       Index j=0;
    //       while (backend_response[j] <= 0) ++j;
    //       Numeric bf_min = backend_f_grid[j];

    //       j=nf-1;
    //       while (backend_response[j] <= 0) --j;
    //       Numeric bf_max = backend_f_grid[j];
    //
    // No, aparently the sensor part want values also where the
    // response is zero. So we simply take the grid boundaries here.
    //
    // We need to add a bit of extra margin at both sides,
    // otherwise there is a numerical problem in the sensor WSMs.
    //
    // PE 081003: The accuracy for me (double on 64 bit machine) appears to
    // be about 3 Hz. Select 1 MHz to have a clear margin. Hopefully OK
    // for other machines.
    //
    // SAB 2010-04-14: The approach with a constant delta does not seem to work
    // well in practice. What I do now is that I add a delta corresponding to a
    // fraction of the grid spacing. But that is done outside of this function.
    // So we set delta = 0 here.
    //
    // SAB 2010-05-03: Now we pass delta as a parameter (with a default value of 0).
    //
    // Isoz 2013-05-21: Added methods to ignore areas between passbands
    //
    const Vector& backend_f_grid =
        backend_channel_response[idx].get_numeric_grid(0);
    const Vector& backend_filter = backend_channel_response[idx].data;
    if (backend_filter.nelem() >=
        4)  // Is the passband frequency response given explicitly ? e.g. [0,>0,>0,0]
    {
      for (Index idy = 1; idy < backend_filter.nelem(); ++idy) {
        if (idy == 1) {
          fmin_pb[pbIdx] = f_backend[idx] + backend_f_grid[0];
        } else if ((backend_filter[idy] > 0) &&
                   (backend_filter[idy - 1] == 0)) {
          fmin_pb[pbIdx] = f_backend[idx] + backend_f_grid[idy - 1] - delta;
        }
        if ((backend_filter[idy] == 0) && (backend_filter[idy - 1] > 0)) {
          fmax_pb[pbIdx] = f_backend[idx] + backend_f_grid[idy] + delta;
          out2 << "  "
               << "fmin_pb " << fmin_pb[pbIdx] << "  "
               << "f_backend" << f_backend[idx] << "  "
               << "fmax_pb " << fmax_pb[pbIdx] << "  "
               << "diff " << fmax_pb[pbIdx] - fmin_pb[pbIdx] << "\n";
          pbIdx++;
        }
      }
      fmax_pb[pbIdx - 1] = f_backend[idx] + last(backend_f_grid);
    } else  // Or are the passbands given implicitly - such as the default for AMSUA and MHS
    {
      fmin_pb[pbIdx] = f_backend[idx] + backend_f_grid[0] - delta;  //delta;
      fmax_pb[pbIdx] = f_backend[idx] +
                       backend_f_grid[backend_f_grid.nelem() - 1] +
                       delta;  ///delta;
      out2 << "  "
           << "fmin_pb " << fmin_pb[pbIdx] << "  "
           << "f_backend" << f_backend[pbIdx] << "  "
           << "fmax_pb " << fmax_pb[pbIdx] << "\n";
      pbIdx++;
    }
  }

  // The problem is that channels may be overlapping. In that case, we
  // want to create a frequency grid that covers their entire range,
  // but we do not want to duplicate any frequencies.

  // To make matters worse, one or even several channels may be
  // completely inside another very broad channel.

  // Sort channels by frequency:
  // Caveat: A channel may be higher in
  // characteristic frequency f_backend, but also wider, so that it
  // has a lower minimum frequency fmin_orig. (This is the case for
  // some HIRS channels.) We sort by the minimum frequency here, not
  // by f_backend. This is necessary for function
  // test_and_merge_two_channels to work correctly.
  ArrayOfIndex isorted;
  get_sorted_indexes(isorted, fmin_pb);

  fmin.resize(numPB);
  fmax.resize(numPB);
  out2 << " resize numPb " << numPB << "\n";
  for (Index idx = 0; idx < numPB; ++idx) {
    fmin[idx] = fmin_pb[isorted[idx]];
    fmax[idx] = fmax_pb[isorted[idx]];
  }
  // We will be testing pairs of channels, and combine them if
  // possible. We have to test always only against the direct
  // neighbour. If that has no overlap, higher channels can not have
  // any either, due to the sorting by fmin.
  //
  // Note that fmin.nelem() changes, as the loop is
  // iterated. Nevertheless this is the correct stop condition.
  for (Index i = 0; i < fmin.nelem() - 1; ++i) {
    bool continue_checking = true;
    // The "i<fmin.nelem()" condition below is necessary, since
    // fmin.nelem() can decrease while the loop is executed, due to
    // merging.
    while (continue_checking && i < fmin.nelem() - 1) {
      continue_checking = test_and_merge_two_channels(fmin, fmax, i, i + 1);

      // Function returns true if merging has taken place.
      // In this case, we have to check again.
    }
  }

  out2 << "  New channel characteristics:\n"
       << "  min                       max (all in Hz):\n";
  for (Index i = 0; i < fmin.nelem(); ++i)
    out2 << "  " << fmin[i] << "               " << fmax[i] << "\n";
}



/*===========================================================================
  === Core integration and sum functions:
  ===========================================================================*/

void integration_func_by_vecmult(VectorView h,
                                 ConstVectorView f,
                                 ConstVectorView x_f_in,
                                 ConstVectorView x_g_in) {
  // Basic sizes
  const Index nf = x_f_in.nelem();
  const Index ng = x_g_in.nelem();

  // Asserts
  ARTS_ASSERT(h.nelem() == ng);
  ARTS_ASSERT(f.nelem() == nf);
  ARTS_ASSERT(is_increasing(x_f_in));
  ARTS_ASSERT(is_increasing(x_g_in) || is_decreasing(x_g_in));
  // More ARTS_ASSERTs below

  // End points of x_f
  Numeric xfmin = x_f_in[0];
  Numeric xfmax = x_f_in[nf - 1];

  // Handle possibly reversed x_g.
  Vector x_g;
  Index xg_reversed = 0;
  //
  if (is_decreasing(x_g_in)) {
    x_g = x_g_in[Range(ng - 1, ng, -1)];
    xg_reversed = 1;
  } else {
    x_g = x_g_in;
  }
  //
  ARTS_ASSERT(x_g[0] <= xfmin);
  ARTS_ASSERT(x_g[ng - 1] >= xfmax);

  // Normalise grids so x_f covers [0,1]
  const Numeric df = xfmax - xfmin;
  Vector x_f(nf);
  //
  for (Index i = 0; i < nf; i++) {
    x_f[i] = (x_f_in[i] - xfmin) / df;
  }
  for (Index i = 0; i < ng; i++) {
    x_g[i] = (x_g[i] - xfmin) / df;
  }
  xfmin = 0;
  xfmax = 1;
  // To test without normalisation, comment out lines above and use:
  //const Numeric df  = 1;
  //const Vector  x_f = x_f_in;

  // Create a reference grid vector, x_ref that containing the values
  // of x_f and x_g strictly sorted. Only g points inside the f range
  // are of concern.
  list<Numeric> l_x;
  for (Index it = 0; it < nf; it++) {
    l_x.push_back(x_f[it]);
  }
  for (Index it = 0; it < ng; it++) {
    if (x_g[it] > xfmin && x_g[it] < xfmax) {
      l_x.push_back(x_g[it]);
    }
  }
  l_x.sort();
  l_x.unique();
  //
  Vector x_ref(l_x.size());
  Index e = 0;
  for (list<Numeric>::iterator li = l_x.begin(); li != l_x.end(); li++) {
    x_ref[e] = *li;
    e++;
  }

  // Initiate output vector, with equal size as x_g, with zeros.
  // Start calculations
  h = 0.0;
  Index i_f = 0;
  Index i_g = 0;
  Numeric dx, a0, b0, c0, a1, b1, c1, x3, x2, x1;
  //
  for (Index i = 0; i < x_ref.nelem() - 1; i++) {
    // Find for what index in x_g (which is the same as for h) and f
    // calculation corresponds to
    while (x_g[i_g + 1] <= x_ref[i]) {
      i_g++;
    }
    while (x_f[i_f + 1] <= x_ref[i]) {
      i_f++;
    }

    // If x_ref[i] is out of x_f's range then that part of the integral is 0,
    // and no calculations should be done
    if (x_ref[i] >= xfmin && x_ref[i] < xfmax) {
      // Product of steps in x_f and x_g
      dx = (x_f[i_f + 1] - x_f[i_f]) * (x_g[i_g + 1] - x_g[i_g]);

      // Calculate a, b and c coefficients; h[i]=ax^3+bx^2+cx
      a0 = (f[i_f] - f[i_f + 1]) / 3.0;
      b0 = (-f[i_f] * (x_g[i_g + 1] + x_f[i_f + 1]) +
            f[i_f + 1] * (x_g[i_g + 1] + x_f[i_f])) /
           2.0;
      c0 = x_g[i_g + 1] * (f[i_f] * x_f[i_f + 1] - f[i_f + 1] * x_f[i_f]);

      a1 = -a0;
      b1 = (f[i_f] * (x_g[i_g] + x_f[i_f + 1]) -
            f[i_f + 1] * (x_g[i_g] + x_f[i_f])) /
           2.0;
      c1 = x_g[i_g] * (-f[i_f] * x_f[i_f + 1] + f[i_f + 1] * x_f[i_f]);

      x1 = x_ref[i + 1] - x_ref[i];
      // Simple way, but sensitive to numerical problems:
      //x2 = pow(x_ref[i+1],2) - pow(x_ref[i],2);
      //x3 = pow(x_ref[i+1],3) - pow(x_ref[i],3);
      // The same but a numerically better way:
      x2 = x1 * (2 * x_ref[i] + x1);
      x3 = x1 * (3 * x_ref[i] * (x_ref[i] + x1) + x1 * x1);

      // Calculate h[i] and h[i+1] increment
      // (df-factor to compensate for normalisation)
      h[i_g] += df * (a0 * x3 + b0 * x2 + c0 * x1) / dx;
      h[i_g + 1] += df * (a1 * x3 + b1 * x2 + c1 * x1) / dx;
    }
  }

  // Flip back if x_g was decreasing
  if (xg_reversed) {
    Vector tmp = h[Range(ng - 1, ng, -1)];  // Flip order
    h = tmp;
  }

  // The expressions are sensitive to numerical issues if two points in x_ref
  // are very close compared to the values in x_ref. A test trying to warn when
  // this happens:
  ARTS_USER_ERROR_IF (min(f) >= 0 && min(h) < -1e-15,
        "Significant negative response value obtained, "
        "despite sensor reponse is strictly positive. This "
        "indicates numerical problems. Is there any very "
        "small spacing of the sensor response grid?"
        "Please, send a report to Patrick if you see this!");
}



void integration_bin_by_vecmult(VectorView h,
                                ConstVectorView x_g_in,
                                const Numeric& limit1,
                                const Numeric& limit2) {
  // Basic sizes
  const Index ng = x_g_in.nelem();

  // Asserts
  ARTS_ASSERT(ng > 1);
  ARTS_ASSERT(h.nelem() == ng);
  ARTS_ASSERT(limit1 <= limit2);

  // Handle possibly reversed x_g.
  Vector x_g;
  Index xg_reversed = 0;
  //
  if (is_decreasing(x_g_in)) {
    x_g = x_g_in[Range(ng - 1, ng, -1)];
    xg_reversed = 1;
  } else {
    x_g = x_g_in;
  }
  //
  ARTS_ASSERT(x_g[0] <= limit1);
  ARTS_ASSERT(x_g[ng - 1] >= limit2);

  // Handle extreme cases
  // Bin has zero width
  if (limit1 == limit2) {
    h = 0.0;
    return;
  }
  // Bin covers complete x_g:
  else if (limit1 == x_g[0] && limit2 == x_g[ng - 1]) {
    h[0] = (x_g[1] - x_g[0]) / 2.0;
    for (Index i = 1; i < ng - 1; i++) {
      h[i] = (x_g[i + 1] - x_g[i - 1]) / 2.0;
    }
    h[ng - 1] = (x_g[ng - 1] - x_g[ng - 2]) / 2.0;
    return;
  }

  // The general case
  Numeric x1 = 0,
          x2 =
              0;  // End points of bin, inside basis range of present grid point
  for (Index i = 0; i < ng; i++) {
    bool inside = false;

    if (i == 0) {
      if (limit1 < x_g[1]) {
        inside = true;
        x1 = limit1;
        x2 = min(limit2, x_g[1]);
      }
    } else if (i == ng - 1) {
      if (limit2 > x_g[ng - 2]) {
        inside = true;
        x1 = max(limit1, x_g[ng - 2]);
        x2 = limit2;
      }
    } else {
      if ((limit1 < x_g[i + 1] && limit2 > x_g[i - 1]) ||
          (limit2 > x_g[i - 1] && limit1 < x_g[i + 1])) {
        inside = true;
        x1 = max(limit1, x_g[i - 1]);
        x2 = min(limit2, x_g[i + 1]);
      }
    }

    // h is zero if no overlap between bin and basis range of point
    if (!inside) {
      h[i] = 0.0;
    } else {
      // Lower part
      if (x1 < x_g[i]) {
        const Numeric r = 1.0 / (x_g[i] - x_g[i - 1]);
        const Numeric y1 = r * (x1 - x_g[i - 1]);
        const Numeric dx = min(x2, x_g[i]) - x1;
        const Numeric y2 = y1 + r * dx;
        h[i] = 0.5 * dx * (y1 + y2);
      } else {
        h[i] = 0.0;
      }

      // Upper part
      if (x2 > x_g[i]) {
        const Numeric r = 1.0 / (x_g[i + 1] - x_g[i]);
        const Numeric y2 = r * (x_g[i + 1] - x2);
        const Numeric dx = x2 - max(x1, x_g[i]);
        const Numeric y1 = y2 + r * dx;
        h[i] += 0.5 * dx * (y1 + y2);
      }
    }
  }

  // Flip back if x_g was decreasing
  if (xg_reversed) {
    Vector tmp = h[Range(ng - 1, ng, -1)];  // Flip order
    h = tmp;
  }
}


void summation_by_vecmult(VectorView h,
                          ConstVectorView f,
                          ConstVectorView x_f,
                          ConstVectorView x_g,
                          const Numeric x1,
                          const Numeric x2) {
  // Asserts
  ARTS_ASSERT(h.nelem() == x_g.nelem());
  ARTS_ASSERT(f.nelem() == x_f.nelem());
  ARTS_ASSERT(x_g[0] <= x_f[0]);
  ARTS_ASSERT(last(x_g) >= last(x_f));
  ARTS_ASSERT(x1 >= x_f[0]);
  ARTS_ASSERT(x2 >= x_f[0]);
  ARTS_ASSERT(x1 <= last(x_f));
  ARTS_ASSERT(x2 <= last(x_f));

  // Determine grid positions for point 1 (both with respect to f and g grids)
  // and interpolate response function.
  ArrayOfGridPos gp1g(1), gp1f(1);
  gridpos(gp1g, x_g, x1);
  gridpos(gp1f, x_f, x1);
  Matrix itw1(1, 2);
  interpweights(itw1, gp1f);
  Numeric f1;
  interp(f1, itw1, f, gp1f);

  // Same for point 2
  ArrayOfGridPos gp2g(1), gp2f(1);
  gridpos(gp2g, x_g, x2);
  gridpos(gp2f, x_f, x2);
  Matrix itw2(1, 2);
  interpweights(itw2, gp2f);
  Numeric f2;
  interp(f2, itw2, f, gp2f);

  //Initialise h at zero and store calculated weighting components
  h = 0.0;
  h[gp1g[0].idx] += f1 * gp1g[0].fd[1];
  h[gp1g[0].idx + 1] += f1 * gp1g[0].fd[0];
  h[gp2g[0].idx] += f2 * gp2g[0].fd[1];
  h[gp2g[0].idx + 1] += f2 * gp2g[0].fd[0];
}
