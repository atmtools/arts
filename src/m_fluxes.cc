/* Copyright (C) 2018
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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



/*===========================================================================
  ===  File description
  ===========================================================================*/
#include <stdexcept>
#include <iostream>
#include "matpackVII.h"
#include "workspace_ng.h"
#include "messages.h"
#include "legendre.h"
#include "math_funcs.h"
#include "sorting.h"

/*!
  \file   M_fluxes.cc
  \author Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2018-06-22

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/


extern const Numeric PI;
extern const Numeric DEG2RAD;


/*===========================================================================
  === The functions
  ===========================================================================*/



/* Workspace method: Doxygen documentation will be auto-generated */
void AngularGridsSetFluxCalc(
        Vector &scat_za_grid,
        Vector &scat_aa_grid,
        Vector &za_grid_weights,
        // Keywords:
        const Index &N_za_grid,
        const Index &N_aa_grid,
        const String &za_grid_type,
        const Verbosity &)
{
    // Azimuth angle grid
    if (N_aa_grid > 1)
        nlinspace(scat_aa_grid, 0, 360, N_aa_grid);
    else if (N_aa_grid < 1)
    {
        ostringstream os;
        os << "N_aa_grid must be > 0 (even for 1D).";
        throw std::runtime_error(os.str());
    } else
    {
        scat_aa_grid.resize(1);
        scat_aa_grid[0] = 0.;
    }

    if (N_za_grid % 2 == 1)
    {
        ostringstream os;
        os << "N_za_grid must be even.";
        throw runtime_error(os.str());
    }

    Index nph = N_za_grid / 2;


    //calculate zenith angle grid
    scat_za_grid.resize(N_za_grid);
    scat_za_grid = 0.;
    za_grid_weights.resize(N_za_grid);
    za_grid_weights = 0;

    if (za_grid_type == "double_gauss")
    {

        Vector x;
        Vector w;
        Vector xtemp;
        Vector wtemp;
        //Numeric theta;

        //calculate legendre weights and evaluation points
        gsl_integration_glfixed_table_alloc(xtemp, wtemp, nph);

        x.resize(nph);
        w.resize(nph);



        // reorder and expand weights and abscissa vectors
        // transform abscissa vector from cos(theta)-space to theta-space
        // and adjust the domain and weights
        if (nph % 2 == 1)
        {

            x[xtemp.nelem() - 1] = acos((xtemp[0] + 1) / 2) / DEG2RAD;
            w[wtemp.nelem() - 1] = wtemp[0] / 2;

            for (Index i = 0; i < xtemp.nelem() - 1; i++)
            {
                x[i] = acos((xtemp[xtemp.nelem() - 1 - i] + 1) / 2.) / DEG2RAD;
                x[xtemp.nelem() + i] = acos(1 - (xtemp[i + 1] + 1) / 2.) / DEG2RAD;

                w[i] = wtemp[wtemp.nelem() - 1 - i] / 2;
                w[wtemp.nelem() + i] = wtemp[i + 1] / 2;
            }
        } else
        {
            for (Index i = 0; i < xtemp.nelem(); i++)
            {
                x[i] = acos((xtemp[xtemp.nelem() - 1 - i] + 1) / 2.) / DEG2RAD;
                x[xtemp.nelem() + i] = acos(1 - (xtemp[i] + 1) / 2.) / DEG2RAD;

                w[i] = wtemp[wtemp.nelem() - 1 - i] / 2;
                w[wtemp.nelem() + i] = wtemp[i] / 2;
            }
        }


        for (Index i = 0; i < nph; i++)
        {
            //set the angles
            //theta=x[i];//acos((x[i]+1)/2)/DEG2RAD;
            scat_za_grid[i] = x[i];
            scat_za_grid[scat_za_grid.nelem() - 1 - i] = 180 - x[i];

            // set the weights to the right component
            za_grid_weights[i] = w[i];
            za_grid_weights[za_grid_weights.nelem() - 1 - i] = w[i];
        }

    } else if (za_grid_type == "linear")
    {

        Vector x;
        Vector w;
        calculate_weights_linear(x, w, nph);

        for (Index i = 0; i < N_za_grid; i++)
        {
            scat_za_grid[i] = (x[i] + 1) * 90.;

            // set the weights to the right component
            // by adjusting the domain, we also have to adjust the weights
            za_grid_weights[i] = w[i] * sin(scat_za_grid[i] * DEG2RAD);

        }
    } else if (za_grid_type == "linear_mu")
    {

        Vector x;
        Vector w;

        //calculate weights in cos(theta) space
        calculate_weights_linear(x, w, nph);

        //allocate
        Vector scat_za_grid_temp;
        scat_za_grid_temp.resize(x.nelem());


        for (Index i = 0; i < N_za_grid; i++)
        {
            scat_za_grid_temp[i] = acos(x[i]) / DEG2RAD;
        }

        //#sort weights and theta in increasing direction of scat_za_grid
        scat_za_grid = scat_za_grid_temp[Range(x.nelem() - 1, x.nelem(), -1)];
        za_grid_weights = w[Range(x.nelem() - 1, x.nelem(), -1)];


    } else
    {
        ostringstream os;
        os << "The selected grid type is not implemented";
        throw std::runtime_error(os.str());
    }

    //be sure that the first and the last angle are within the closed interval
    //between 0 and 180 deg, because ARTS is picky if the angles are due to numerics
    // slightly below and above,respectively.
    if (scat_za_grid[0] < 0)
    {
        scat_za_grid[0] = 0.;
    }

    if (scat_za_grid[scat_za_grid.nelem() - 1] > 180)
    {
        scat_za_grid[scat_za_grid.nelem() - 1] = 180.;
    }


}


/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_irradiance_fieldFromiyField(
        Tensor5 &spectral_irradiance_field,
        const Tensor7 &doit_i_field,
        const Vector &scat_za_grid,
        const Vector &scat_aa_grid,
        const Vector &za_grid_weights,
        const Verbosity &)
{

    // Number of zenith angles.
    const Index N_scat_za = scat_za_grid.nelem();
    const Index N_scat_aa = scat_aa_grid.nelem();


    Tensor5 iy_field_aa_integrated;

    //azimuth integration
    if (N_scat_aa == 1)  //1D case no azimuth dependency
    {

        iy_field_aa_integrated = doit_i_field(joker, joker, joker, joker, joker, 0, 0);
        iy_field_aa_integrated *= 2 * PI;

    } else //general case with azimuth dependency
    {

        iy_field_aa_integrated.resize(doit_i_field.nlibraries(), doit_i_field.nvitrines(),
                                      doit_i_field.nshelves(), doit_i_field.nbooks(), doit_i_field.npages());
        iy_field_aa_integrated = 0.;

        for (Index s = 0; s < iy_field_aa_integrated.nshelves(); s++)
        {
            for (Index b = 0; b < iy_field_aa_integrated.nbooks(); b++)
            {
                for (Index p = 0; p < iy_field_aa_integrated.npages(); p++)
                {
                    for (Index r = 0; r < iy_field_aa_integrated.nrows(); r++)
                    {
                        for (Index c = 0; c < iy_field_aa_integrated.ncols(); c++)
                        {
                            for (Index i = 0; i < N_scat_aa - 1; i++)
                            {
                                iy_field_aa_integrated(s, b, p, r, c) +=
                                        (doit_i_field(s, b, p, r, c, i, 0) +
                                         doit_i_field(s, b, p, r, c, i + 1, 0)) / 2. *
                                        abs(scat_aa_grid[i + 1] - scat_aa_grid[i]) * DEG2RAD;
                            }
                        }
                    }
                }
            }
        }
    }

    //allocate
    spectral_irradiance_field.resize(doit_i_field.nlibraries(), doit_i_field.nvitrines(),
                                     doit_i_field.nshelves(), doit_i_field.nbooks(), 2);
    spectral_irradiance_field = 0;


    // zenith angle integration
    for (Index s = 0; s < spectral_irradiance_field.nshelves(); s++)
    {
        for (Index b = 0; b < spectral_irradiance_field.nbooks(); b++)
        {
            for (Index p = 0; p < spectral_irradiance_field.npages(); p++)
            {
                for (Index r = 0; r < spectral_irradiance_field.nrows(); r++)
                {
                    for (Index i = 0; i < N_scat_za; i++)
                    {

                        if (scat_za_grid[i] <= 90.)
                        {
                            spectral_irradiance_field(s, b, p, r, 0) += iy_field_aa_integrated(s, b, p, r, i) *
                                                                        cos(scat_za_grid[i] * DEG2RAD) * (-1.) *
                                                                        za_grid_weights[i];
                        } else
                        {
                            spectral_irradiance_field(s, b, p, r, 1) += iy_field_aa_integrated(s, b, p, r, i) *
                                                                        cos(scat_za_grid[i] * DEG2RAD) * (-1.) *
                                                                        za_grid_weights[i];
                        }
                    }
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void irradiance_fieldFromRadiance(
        Tensor4 &irradiance_field,
        const Tensor5 &radiance_field,
        const Vector &scat_za_grid,
        const Vector &scat_aa_grid,
        const Vector &za_grid_weights,
        const Verbosity &)
{

    // Number of zenith angles.
    const Index N_scat_za = scat_za_grid.nelem();
    const Index N_scat_aa = scat_aa_grid.nelem();


    Tensor4 radiance_field_aa_integrated;

    //azimuth integration
    if (N_scat_aa == 1)  //1D case no azimuth dependency
    {

        radiance_field_aa_integrated = radiance_field(joker, joker, joker, joker, 0);
        radiance_field_aa_integrated *= 2 * PI;

    } else //general case with azimuth dependency
    {

        radiance_field_aa_integrated.resize(radiance_field.nshelves(), radiance_field.nbooks(),
                                            radiance_field.npages(), radiance_field.nrows());
        radiance_field_aa_integrated = 0.;


        for (Index b = 0; b < radiance_field_aa_integrated.nbooks(); b++)
        {
            for (Index p = 0; p < radiance_field_aa_integrated.npages(); p++)
            {
                for (Index r = 0; r < radiance_field_aa_integrated.nrows(); r++)
                {
                    for (Index c = 0; c < radiance_field_aa_integrated.ncols(); c++)
                    {
                        for (Index i = 0; i < N_scat_aa - 1; i++)
                        {
                            radiance_field_aa_integrated(b, p, r, c) +=
                                    (radiance_field(b, p, r, c, i) +
                                     radiance_field(b, p, r, c, i + 1)) / 2. *
                                    abs(scat_aa_grid[i + 1] - scat_aa_grid[i]) * DEG2RAD;
                        }
                    }
                }
            }
        }

    }

    //allocate
    irradiance_field.resize(radiance_field.nshelves(), radiance_field.nbooks(), radiance_field.npages(), 2);
    irradiance_field = 0;


    // zenith angle integration

    for (Index b = 0; b < irradiance_field.nbooks(); b++)
    {
        for (Index p = 0; p < irradiance_field.npages(); p++)
        {
            for (Index r = 0; r < irradiance_field.nrows(); r++)
            {
                for (Index i = 0; i < N_scat_za; i++)
                {

                    if (scat_za_grid[i] <= 90.)
                    {
                        irradiance_field(b, p, r, 0) += radiance_field_aa_integrated(b, p, r, i) *
                                                        cos(scat_za_grid[i] * DEG2RAD) * (-1.) *
                                                        za_grid_weights[i];
                    } else
                    {
                        irradiance_field(b, p, r, 1) += radiance_field_aa_integrated(b, p, r, i) *
                                                        cos(scat_za_grid[i] * DEG2RAD) * (-1.) *
                                                        za_grid_weights[i];
                    }
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void RadiationFieldSpectralIntegrate(
        Tensor4 &radiation_field,
        const Vector &f_grid,
        const Tensor5 &spectral_radiation_field,
        const Verbosity &)
{

    if (f_grid.nelem() != spectral_radiation_field.nshelves())
    {
        throw runtime_error("The length of f_grid does not match with\n"
                            " the first dimension of the spectral_radiation_field");
    }


    //allocate
    radiation_field.resize(spectral_radiation_field.nbooks(), spectral_radiation_field.npages(),
                           spectral_radiation_field.nrows(), spectral_radiation_field.ncols());
    radiation_field = 0;


    // frequency integration
    for (Index i = 0; i < spectral_radiation_field.nshelves() - 1; i++)
    {
        const Numeric df = f_grid[i + 1] - f_grid[i];

        for (Index b = 0; b < radiation_field.nbooks(); b++)
        {
            for (Index p = 0; p < radiation_field.npages(); p++)
            {
                for (Index r = 0; r < radiation_field.nrows(); r++)
                {
                    for (Index c = 0; c < radiation_field.ncols(); c++)
                    {

                        radiation_field(b, p, r, c) += (spectral_radiation_field(i + 1, b, p, r, c) +
                                                        spectral_radiation_field(i, b, p, r, c)) / 2 * df;

                    }
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void RadiationFieldSpectralIntegrate(
        Tensor5 &radiation_field,
        const Vector &f_grid,
        const Tensor7 &spectral_radiation_field,
        const Verbosity &)
{

    if (f_grid.nelem() != spectral_radiation_field.nlibraries())
    {
        throw runtime_error("The length of f_grid does not match with\n"
                            " the first dimension of the spectral_radiation_field");
    }


    //allocate
    radiation_field.resize(spectral_radiation_field.nvitrines(), spectral_radiation_field.nshelves(),
                           spectral_radiation_field.nbooks(), spectral_radiation_field.npages(),
                           spectral_radiation_field.nrows());
    radiation_field = 0;


    // frequency integration
    for (Index i = 0; i < spectral_radiation_field.nlibraries() - 1; i++)
    {
        const Numeric df = f_grid[i + 1] - f_grid[i];

        for (Index s = 0; s < radiation_field.nshelves(); s++)
        {
            for (Index b = 0; b < radiation_field.nbooks(); b++)
            {
                for (Index p = 0; p < radiation_field.npages(); p++)
                {
                    for (Index r = 0; r < radiation_field.nrows(); r++)
                    {
                        for (Index c = 0; c < radiation_field.ncols(); c++)
                        {

                            radiation_field(s, b, p, r, c) += (spectral_radiation_field(i + 1, s, b, p, r, c, 0) +
                                                               spectral_radiation_field(i, s, b, p, r, c, 0)) / 2 * df;


                        }
                    }
                }
            }
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void heating_ratesFromIrradiance(
        Tensor3 &heating_rates,
        const Vector &p_grid,
        const Tensor4 &irradiance_field,
        const Tensor3 &specific_heat_capacity,
        const Numeric &g0,
        const Verbosity &)
{


    //allocate
    heating_rates.resize(irradiance_field.nbooks(), irradiance_field.npages(),
                         irradiance_field.nrows());
    heating_rates = 0;

    // allocate some auxiliary variables
    Numeric net_flux_b; //net_flux bottom
    Numeric net_flux_c; //net_flux center
    Numeric net_flux_t; //net_flux top
    Index idx;


    // calculate heating rates, we skip the upper and lower boundary here, because to achieve the same
    // second order accuracy for the derivation of the net flux at the edged, we use
    // a differentiation based on polynomial interpolation
    for (Index b = 1; b < irradiance_field.nbooks() - 1; b++)
    {
        for (Index p = 0; p < irradiance_field.npages(); p++)
        {
            for (Index r = 0; r < irradiance_field.nrows(); r++)
            {

                net_flux_b = (irradiance_field(b - 1, p, r, 0) + irradiance_field(b - 1, p, r, 1));
                net_flux_t = (irradiance_field(b + 1, p, r, 0) + irradiance_field(b + 1, p, r, 1));

                heating_rates(b, p, r) = (net_flux_t - net_flux_b) /
                                         (p_grid[b + 1] - p_grid[b - 1]) *
                                         g0 / specific_heat_capacity(b, p, r);


            }
        }
    }

    idx = irradiance_field.nbooks();

    // now calculate the heating rates for the upper and lower boundary
    for (Index p = 0; p < irradiance_field.npages(); p++)
    {
        for (Index r = 0; r < irradiance_field.nrows(); r++)
        {
            // lower boundary
            net_flux_b = (irradiance_field(0, p, r, 0) + irradiance_field(0, p, r, 1));
            net_flux_c = (irradiance_field(1, p, r, 0) + irradiance_field(1, p, r, 1));
            net_flux_t = (irradiance_field(2, p, r, 0) + irradiance_field(0, p, r, 1));

            heating_rates(0, p, r) = (-3 * net_flux_b + 4 * net_flux_c - net_flux_t) / (p_grid[2] - p_grid[0]) *
                                     g0 / specific_heat_capacity(0, p, r);


            // lower boundary
            net_flux_t = (irradiance_field(idx - 1, p, r, 0) + irradiance_field(idx - 1, p, r, 1));
            net_flux_c = (irradiance_field(idx - 2, p, r, 0) + irradiance_field(idx - 2, p, r, 1));
            net_flux_b = (irradiance_field(idx - 3, p, r, 0) + irradiance_field(idx - 3, p, r, 1));

            heating_rates(idx - 1, p, r) = -(-3 * net_flux_t + 4 * net_flux_c - net_flux_b) / (p_grid[2] - p_grid[0]) *
                                           g0 / specific_heat_capacity(0, p, r);


        }
    }
}