/* Copyright (C) 2018 Oliver Lemke <oliver.lemke@uni-hamburg.de>

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
  \file   hitran_xsec.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/


#include "arts.h"

#ifdef ENABLE_FFTW

#include <fftw3.h>
#include <complex.h>

#endif  /* ENABLE_FFTW */

#include "hitran_xsec.h"
#include "absorption.h"
#include "check_input.h"


extern const Numeric PI;

Numeric func_2straights(const Numeric x, const Vector& coeffs)
{
    assert(coeffs.nelem() == 3);
    return (x <= coeffs[0]) ? coeffs[1] * x
                            : coeffs[2] * (x - coeffs[0]) +
                              coeffs[1] * coeffs[0];
}

Numeric lorentz_pdf(const Numeric x, const Numeric x0, const Numeric gamma)
{
    const Numeric xx0 = x - x0;
    return gamma / PI / (xx0 * xx0 + gamma * gamma);

}

String XsecRecord::SpeciesName() const
{
    // The function species_name_from_species_index internally does an assertion
    // that the species with this index really exists.
    return species_name_from_species_index(mspecies);
}


void convolve(Vector& result, ConstVectorView& xsec, ConstVectorView& lorentz)
{
    Index n_xsec = xsec.nelem();
    Index n_lorentz = lorentz.nelem();
    //    assert(n_xsec == n_lorentz);
    Vector temp(n_xsec + n_lorentz - 1);

    for (Index i = 0; i < n_xsec + n_lorentz - 1; ++i)
    {
        Numeric sum = 0.0;
        for (Index j = 0; j <= i; ++j)
        {
            sum += ((j < n_xsec) && (i - j < n_lorentz)) ? xsec[j] *
                                                           lorentz[i - j] : 0.0;
        }
        temp[i] = sum;
    }
    result = temp[Range(n_lorentz / 2, n_xsec, 1)];
}

#ifdef ENABLE_FFTW

void fftconvolve(VectorView& result, const Vector& xsec, const Vector& lorentz)
{
    int n_p = (int) (xsec.nelem() + lorentz.nelem() - 1);
    int n_p_2 = n_p / 2 + 1;

    double *xsec_in = fftw_alloc_real((size_t)n_p);
    fftw_complex *xsec_out = fftw_alloc_complex((size_t)n_p_2);
    memcpy(xsec_in, xsec.get_c_array(), sizeof(double) * xsec.nelem());
    memset(&xsec_in[xsec.nelem()], 0, sizeof(double) * (n_p - xsec.nelem()));

    fftw_plan plan;
#pragma omp critical(fftw_call)
    plan = fftw_plan_dft_r2c_1d(n_p, xsec_in, xsec_out, FFTW_ESTIMATE);

    fftw_execute(plan);

#pragma omp critical(fftw_call)
    fftw_destroy_plan(plan);

    fftw_free(xsec_in);

    double *lorentz_in = fftw_alloc_real((size_t)n_p);
    fftw_complex *lorentz_out = fftw_alloc_complex((size_t)n_p_2);
    memcpy(lorentz_in, lorentz.get_c_array(), sizeof(double) * lorentz.nelem());
    memset(&lorentz_in[lorentz.nelem()], 0,
           sizeof(double) * (n_p - lorentz.nelem()));

#pragma omp critical(fftw_call)
    plan = fftw_plan_dft_r2c_1d(n_p, lorentz_in, lorentz_out, FFTW_ESTIMATE);

    fftw_execute(plan);

#pragma omp critical(fftw_call)
    fftw_destroy_plan(plan);

    fftw_free(lorentz_in);

    fftw_complex *fft_in = fftw_alloc_complex((size_t)n_p_2);
    double *fft_out = fftw_alloc_real((size_t)n_p);
    memcpy(fft_in, xsec_out, sizeof(fftw_complex) * n_p_2);

    for (Index i = 0; i < n_p_2; i++)
    {
        fft_in[i][0] = xsec_out[i][0] * lorentz_out[i][0] -
                       xsec_out[i][1] * lorentz_out[i][1];
        fft_in[i][1] = xsec_out[i][0] * lorentz_out[i][1] +
                       xsec_out[i][1] * lorentz_out[i][0];
    }

#pragma omp critical(fftw_call)
    plan = fftw_plan_dft_c2r_1d(n_p, fft_in, fft_out, FFTW_ESTIMATE);

    fftw_execute(plan);

#pragma omp critical(fftw_call)
    fftw_destroy_plan(plan);

    fftw_free(fft_in);

    for (Index i = 0; i < xsec.nelem(); i++)
    {
        result[i] = fft_out[i + (int) lorentz.nelem() / 2] / n_p;
    }

    fftw_free(xsec_out);
    fftw_free(lorentz_out);
    fftw_free(fft_out);
}

#endif  /* ENABLE_FFTW */


void XsecRecord::Extract(VectorView result,
                         ConstVectorView f_grid,
                         const Numeric& pressure,
                         const Numeric& temperature,
                         const Index& apply_tfit,
                         const Verbosity& verbosity) const
{
    CREATE_OUTS;

    const Index nf = f_grid.nelem();

    // Assert that result vector has right size:
    assert(result.nelem() == nf);

    // Initialize result to zero (important for those frequencies outside the data grid).
    result = 0.;

    const Index ndatasets = mxsecs.nelem();
    for (Index this_dataset_i = 0; this_dataset_i < ndatasets; this_dataset_i++)
    {
        const Vector& data_f_grid = mfgrids[this_dataset_i];
        const Numeric fmin = data_f_grid[0];
        const Numeric fmax = data_f_grid[mfgrids[this_dataset_i].nelem() - 1];
        const Index data_nf = mfgrids[this_dataset_i].nelem();

        if (out3.sufficient_priority())
        {
            // Some detailed information to the most verbose output stream:
            ostringstream os;
            os << "    f_grid:      " << f_grid[0] << " - "
               << f_grid[nf - 1] << " Hz\n"
               << "    data_f_grid: " << fmin << " - " << fmax << " Hz\n"
               << "    pressure: " << pressure << " K\n";
            out3 << os.str();
        }

        // We want to return result zero for all f_grid points that are outside the
        // data_f_grid, because xsec datasets are defined only where the absorption
        // was measured. So, we have to find out which part of f_grid is inside
        // data_f_grid.
        Index i_fstart, i_fstop;

        for (i_fstart = 0; i_fstart < nf; ++i_fstart)
            if (f_grid[i_fstart] >= fmin) break;

        // Return directly if all frequencies are below data_f_grid:
        if (i_fstart == nf) continue;

        for (i_fstop = nf - 1; i_fstop >= 0; --i_fstop)
            if (f_grid[i_fstop] <= fmax) break;

        // Return directly if all frequencies are above data_f_grid:
        if (i_fstop == -1) continue;

        // Extent for active frequency vector:
        const Index f_extent = i_fstop - i_fstart + 1;

        if (out3.sufficient_priority())
        {
            ostringstream os;
            os << "    " << f_extent
               << " frequency extraction points starting at "
               << "frequency index " << i_fstart << ".\n";
            out3 << os.str();
        }

        // If f_extent is less than one, then the entire data_f_grid is between two
        // grid points of f_grid. (So that we do not have any f_grid points inside
        // data_f_grid.) Return also in this case.
        if (f_extent < 3) continue;

        // This is the part of f_grid for which we have to do the interpolation.
        ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];


        // We also need to determine the range in the xsec dataset that's inside
        // f_grid. We can ignore the remaining data.
        Index i_data_fstart, i_data_fstop;

        for (i_data_fstart = 0; i_data_fstart < data_nf; ++i_data_fstart)
            if (data_f_grid[i_data_fstart] >= fmin) break;

        for (i_data_fstop = data_nf - 1; i_data_fstop >= 0; --i_data_fstop)
            if (data_f_grid[i_data_fstop] <= fmax) break;

        // Extent for active data frequency vector:
        const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

        // This is the part of f_grid for which we have to do the interpolation.
        ConstVectorView data_f_grid_active = data_f_grid[Range(i_data_fstart,
                                                               data_f_extent)];

        // This is the part of the xsec dataset for which we have to do the
        // interpolation.
        Range active_range(i_data_fstart, data_f_extent);
        ConstVectorView xsec_active = mxsecs[this_dataset_i][active_range];

        Vector xsec_active_tfit;

        if (apply_tfit != 0)
        {
            xsec_active_tfit = mtslope[this_dataset_i][active_range];
            xsec_active_tfit *= temperature - mreftemperature[this_dataset_i];
            xsec_active_tfit += mtintersect[this_dataset_i][active_range];
            xsec_active_tfit /= 10000;
            xsec_active_tfit += mxsecs[this_dataset_i][active_range];

            xsec_active = xsec_active_tfit;
        }
        // We have to create a matching view on the result vector:
        VectorView result_active = result[Range(i_fstart, f_extent)];
        Vector xsec_interp(f_extent);


        // Decide on interpolation orders:
        const Index f_order = 3;

        // The frequency grid has to have enough points for this interpolation
        // order, otherwise throw a runtime error.
        if (data_f_grid.nelem() < f_order + 1)
        {
            ostringstream os;
            os << "Not enough frequency grid points in Hitran Xsec data.\n"
               << "You have only " << data_f_grid.nelem() << " grid points.\n"
               << "But need at least " << f_order + 1 << ".";
            throw runtime_error(os.str());
        }

        if (pressure > mrefpressure[this_dataset_i])
        {
            // Apply pressure dependent broadening and set negative values to zero.
            // (These could happen due to overshooting of the higher order interpolation.)
            const Numeric pdiff = pressure - mrefpressure[this_dataset_i];
            const Numeric fwhm = func_2straights(pdiff, mcoeffs);
            //        std::cout << mcoeffs << " - ";
            //        std::cout << "pdiff: " << pdiff << " - fwhm: " << fwhm << " - fstep: "
            //                  << f_grid[i_fstart] + f_grid[i_fstart + 1] << std::endl;


            Vector f_lorentz(data_f_extent);
            Numeric lsum = 0.;
            for (Index i = 0; i < data_f_extent; i++)
            {
                f_lorentz[i] = lorentz_pdf(data_f_grid[i_data_fstart + i],
                                           data_f_grid[i_data_fstart +
                                                       data_f_extent / 2],
                                           fwhm / 2.);
                lsum += f_lorentz[i];
            }

            f_lorentz /= lsum;

            Vector data_result(xsec_active.nelem());
#ifdef ENABLE_FFTW
            fftconvolve(data_result, xsec_active,
                        f_lorentz[Range(f_lorentz.nelem() / 4,
                                        f_lorentz.nelem() / 2, 1)]);
#else
            convolve(data_result, xsec_active,
                        f_lorentz[Range(f_lorentz.nelem() / 4,
                                        f_lorentz.nelem() / 2, 1)]);
#endif  /* ENABLE_FFTW */

            // TODO: Add to result_active here
            // Check if frequency is inside the range covered by the data:
            chk_interpolation_grids(
                    "Frequency interpolation for cross sections",
                    data_f_grid,
                    f_grid_active,
                    f_order);

            {
                // Find frequency grid positions:
                ArrayOfGridPosPoly f_gp(f_grid_active.nelem()), T_gp(1);
                gridpos_poly(f_gp, data_f_grid_active, f_grid_active, f_order);

                Matrix itw(f_gp.nelem(), f_order + 1);
                interpweights(itw, f_gp);
                interp(xsec_interp, itw, data_result, f_gp);
            }
        }
        else
        {
            // Find frequency grid positions:
            ArrayOfGridPosPoly f_gp(f_grid_active.nelem()), T_gp(1);
            gridpos_poly(f_gp, data_f_grid_active, f_grid_active, f_order);

            Matrix itw(f_gp.nelem(), f_order + 1);
            interpweights(itw, f_gp);
            interp(xsec_interp, itw, xsec_active, f_gp);
        }

        result_active += xsec_interp;
    }
}


/** Get the index in hitran_xsec_data for the given species.

 \param[in] hitran_xsec_data Hitran Xsec data array
 \param[in] species Species name

 \returns Index of this species in hitran_xsec_data. -1 if not found.
 */
Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const Index species)
{
    for (Index i = 0; i < xsec_data.nelem(); i++)
        if (xsec_data[i].Species() == species)
            return i;

    return -1;
}

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd)
{
    os << "Species: " << xd.Species() << std::endl;
    return os;
}
