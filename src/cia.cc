/* Copyright (C) 2012 Stefan Buehler <sbuehler@ltu.se>
 
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

/** \file
 Implementation file for work with HITRAN collision induced absorption (CIA).
 
 The CIA data are part of the HITRAN distribution. They are described in
 Richard, C., I. E. Gordon, L. S. Rothman, M. Abel, L. Frommhold, M. Gustafsson,
 J.-M. Hartmann, C. Hermans, W. J. Lafferty, G. S. Orton, K.M. Smith, and H. Tran (2012),
 New section of the HITRAN database: Collision-induced absorption (CIA),
 J. Quant. Spectrosc. Radiat. Transfer, 113, 1276-1285, doi:10.1016/j.jqsrt.2011.11.004.
 
 \author Stefan Buehler
 \date   2012-11-30
 */

#include "cia.h"
#include "interpolation_poly.h"
#include "abs_species_tags.h"
#include "file.h"

extern const Numeric SPEED_OF_LIGHT;

/** Interpolate CIA data.
 
 Interpolate CIA data to given frequency vector and given scalar temperature.
 Uses third order interpolation in both coordinates, if grid length allows,
 otherwise lower order or no interpolation.
 
 \param[out] result CIA value for given frequency grid and temperature.
 \param[in] frequency Frequency grid
 \param[in] temperature Scalar temperature
 \param[in] cia_data The CIA dataset to interpolate */
void cia_interpolation(VectorView result,
                       ConstVectorView frequency,
                       const Numeric& temperature,
                       const GriddedField2& cia_data)
{
    // Assert that result vector has right size:
    assert(result.nelem()==frequency.nelem());
    
    // Get data grids:
    ConstVectorView f_grid = cia_data.get_numeric_grid(0);
    ConstVectorView T_grid = cia_data.get_numeric_grid(1);

    // Decide on interpolation orders:
    const Index f_order = 3;
    
    // The frequency grid has to have enough points for this interpolation
    // order, otherwise throw a runtime error.
    if ( f_grid.nelem() < f_order+1 )
      {
        ostringstream os;
        os << "Not enough frequency grid points in CIA data.\n"
           << "You have only " << f_grid.nelem() << " grid points.\n"
           << "But need at least " << f_order+1 << ".";
        throw runtime_error(os.str());
      }


    // For T we have to be adaptive, since sometimes there is only one T in
    // the data
    Index T_order;
    switch (T_grid.nelem()) {
        case 1:
            T_order = 0;
            break;
        case 2:
            T_order = 1;
            break;
        case 3:
            T_order = 2;
            break;
        default:
            T_order = 3;
            break;
    }
    
    // Check if frequency is inside the range covered by the data:
    chk_interpolation_grids("Frequency interpolation for CIA continuum",
                                f_grid,
                                frequency,
                                f_order);
    

    // Check if temperature is inside the range covered by the data:
    if (T_order > 0) {
        chk_interpolation_grids("Temperature interpolation for CIA continuum",
                                T_grid,
                                temperature,
                                T_order);
    }
    
    // Find frequency grid positions:
    ArrayOfGridPosPoly f_gp(frequency.nelem()), T_gp(1);
    gridpos_poly(f_gp, f_grid, frequency, f_order);
 
    // Do the rest of the interpolation.
    if (T_order == 0)
      {
        // No temperature interpolation in this case, just a frequency interpolation.
        
        Matrix itw(f_gp.nelem(),f_order+1);
        interpweights(itw,f_gp);
        interp(result, itw, cia_data.data(joker,0), f_gp);
      }
    else
      {
        // Temperature and frequency interpolation.
        
        // Find temperature grid position:
        gridpos_poly(T_gp, T_grid, temperature, T_order);
    
        // Calculate combined interpolation weights:
        Tensor3 itw(f_gp.nelem(),T_gp.nelem(),(f_order+1)*(T_order+1));
        interpweights(itw,f_gp, T_gp);
        
        // Make a matrix view of the results vector:
        MatrixView result_matrix(result);
        
        // Actual interpolation:
        interp(result_matrix, itw, cia_data.data, f_gp, T_gp);
      }
}


// Documentation in header file.
String CIARecord::MoleculeName(const Index i) const
{
    // Assert that i is 0 or 1:
    assert(i>=0);
    assert(i<=1);

    // The function species_name_from_species_index internally does an assertion
    // that the species with this index really exists.
    return species_name_from_species_index( mspecies[i] );
}

// Documentation in header file.
void CIARecord::SetMoleculeName(const Index i,
                                const String& name)
{
    // Assert that i is 0 or 1:
    assert(i>=0);
    assert(i<=1);
    
    // Find out the species index for name:
    Index spec_ind = species_index_from_species_name(name);

    // Function species_index_from_species_name returns -1 if the species does
    // not exist. Check this:
    if ( spec_ind < 0 )
      {
        ostringstream os;
        os << "Species does not exist in ARTS: " << name;
        throw runtime_error(os.str());
      }

    // Assign species:
    mspecies[i] = spec_ind;
}


//! Read CIA catalog file.
/*!
 Reads the given CIA catalog file into this CIARecord.

 \param[in]  filename  Path of catalog file to read.
 \param[in]  verbosity.
 \return os
 */
void CIARecord::ReadFromCIA(const String& filename, const Verbosity& verbosity)
{
    CREATE_OUT2;

    ifstream is;

    out2 << "  Reading file: " << filename << "\n";
    open_input_file(is, filename);

    // Number of points for spectral range in current dataset
    Index npoints = -1;

    // Min/max wave numbers for this dataset
    Numeric wave_min = -1.;
    Numeric wave_max = -1.;

    // Current dataset index
    Index ndataset = -1;

    // Current set in dataset
    Index nset = -1;

    // Frequency, temp and cia values for current dataset
    Vector freq;
    ArrayOfNumeric temp;
    ArrayOfVector cia;

    // Keep track of current line in file
    Index nline = 0;

    mdata.resize(0);
    while (is)
    {
        String line;

        // Extract needed information from dataset header line
        //////////////////////////////////////////////////////
        nline++;
        getline(is,line);
        if (is.eof()) continue;

        if (line.nelem() < 100)
        {
            ostringstream os;
            os << "Error in line " << nline
            << " reading CIA catalog file " << filename << endl
            << "Header line unexpectedly short: " << endl << line;

            throw runtime_error(os.str());
        }

        if (is.bad())
        {
            ostringstream os;
            os << "Error in line " << nline
            << " reading CIA catalog file " << filename << endl;

            throw runtime_error(os.str());
        }

        line.erase(0, 20);

        istringstream istr(line);

        // Data for current set
        Index set_npoints;
        Numeric set_temp;
        Numeric set_wave_min, set_wave_max;

        istr >> set_wave_min >> set_wave_max >> set_npoints >> set_temp;
        if (!istr)
        {
            ostringstream os;
            os << "Error in line " << nline
            << " reading CIA catalog file " << filename << endl;

            throw runtime_error(os.str());
        }

        // If the min/max wave numbers of this set are different from the
        // previous one, a new dataset starts
        if (npoints == -1 || wave_min != set_wave_min || wave_max != set_wave_max)
        {
            if (ndataset != -1)
                AppendDataset(freq, temp, cia);

            npoints = set_npoints;
            ndataset++;
            wave_min = set_wave_min;
            wave_max = set_wave_max;
            nset = 0;
            freq.resize(set_npoints);
            temp.resize(0);
            cia.resize(0);
        }
        if (npoints != set_npoints)
        {
            ostringstream os;
            os << "Error in line " << nline
            << " reading CIA catalog file " << filename << endl
            << "Inconsistent number of data points. Expected " << npoints
            << ", got " << set_npoints;

            throw runtime_error(os.str());
        }

        temp.push_back(set_temp);
        cia.push_back(Vector(npoints));

        // Read dataset
        ////////////////////
        for (Index i = 0; i < npoints; i++)
        {
            Numeric w, c;

            nline++;
            getline(is,line);
            istr.str(line);
            istr >> w >> c;

            if (is.bad() || istr.bad())
            {
                ostringstream os;
                os << "Error in line " << nline
                << " reading CIA catalog file " << filename << ":" << endl
                << line;

                throw runtime_error(os.str());
            }

            // Convert wavenumbers to Herz
            freq[i] = 100. * w * SPEED_OF_LIGHT;

            cia[nset][i] = c;
        }

        nset++;
    }

    if (is.bad())
    {
        ostringstream os;
        os << "Error in line " << nline
        << " reading CIA catalog file " << filename << endl;

        throw runtime_error(os.str());
    }

    AppendDataset(freq, temp, cia);

//    // For debugging
//    for(Index i = 0; i < mdata.nelem(); i++)
//    {
//        cout << i << " ";
//        cout << mdata[i].get_numeric_grid(0).nelem() << " ";
//        cout << mdata[i].get_numeric_grid(1).nelem() << endl;
//    }
}


/** Append data dataset to mdata. */
void CIARecord::AppendDataset(const Vector& freq,
                              const ArrayOfNumeric& temp,
                              const ArrayOfVector& cia)
{
    GriddedField2 dataset;
    dataset.resize(freq.nelem(), temp.nelem());
    dataset.set_grid(0, freq);
    dataset.set_grid_name(0, "Frequency");

    Vector temp_t;
    temp_t = temp;
    dataset.set_grid(1, temp_t);
    dataset.set_grid_name(0, "Temperature");

    for (Index t = 0; t < temp.nelem(); t++)
        dataset.data(joker, t) = cia[t];
    mdata.push_back(dataset);
}


//! Output operator for CIARecord
/*!
 Outputs the grids for the given CIARecord.

 \param[in,out]  os  Output stream.
 \param[in]      cr  CIARecord.
 \return os
 */
ostream& operator<<(ostream& os, const CIARecord& /* cr */)
{
    os << "CIARecord output operator not yet implemented." << endl;
    return os;
}

