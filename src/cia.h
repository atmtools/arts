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
    Header file for work with HITRAN collision induced absorption (CIA). 

    The CIA data are part of the HITRAN distribution. They are described in
    Richard, C., I. E. Gordon, L. S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, 
    J.-M. Hartmann, C. Hermans, W. J. Lafferty, G. S. Orton, K.M. Smith, and H. Tran (2012), 
    New section of the HITRAN database: Collision-induced absorption (CIA),
    J. Quant. Spectrosc. Radiat. Transfer, 113, 1276-1285, doi:10.1016/j.jqsrt.2011.11.004.

    \author Stefan Buehler
    \date   2012-11-30
*/

#ifndef cia_h
#define cia_h

#include "arts.h"
#include "matpackI.h"
#include "mystring.h"
#include "gridded_fields.h"
#include "check_input.h"


class CIARecord;

typedef Array<CIARecord> ArrayOfCIARecord;


/* Header with implementation. */
void cia_interpolation(VectorView result,
                       ConstVectorView frequency,
                       const Numeric& temperature,
                       const GriddedField2& cia_data);

Index cia_get_index(ArrayOfCIARecord cia_data, Index sp1, Index sp2);


/** CIA data for a single pair of molecules.
 
 A variable of this class can hold the complete information from one HITRAN CIA file.
 A HITRAN CIA data file can hold several datasets (data for different temperatures but
 fixed frequency range). But all datasets are for the same pair of molecules.
 
 \author Stefan Buehler
 \date   2000-08-21  */
class CIARecord {
    
public:
    /** Return each molecule name (as a string) that is associated with this CIARecord.
     
     The CIARecord is defined for a pair of molecules!
     
     \param[in] i Must be either 0 or 1. Then the first or second name of the pair
                  is returned.
     
     */
    String MoleculeName(const Index i) const;
    
    
    /** Set each molecule name (from a string) that is associated with this CIARecord.
     
     The CIARecord is defined for a pair of molecules. The molecule names are 
     internally stored as species indices.
     
     \param[in] i Must be either 0 or 1. Then the first or second name of the pair
     is returned.
     \param[in] name The molecule name as a string, e.g., "H2O".
     
     */
    void SetMoleculeName(const Index i,
                         const String& name);

     
    /** Return CIA species index.
     
     \param[in] i Must be either 0 or 1. Then the first or second species index
     is returned.
     */
    Index Species(const Index i) const
    {
      // Assert that i is 0 or 1:
      assert(i>=0);
      assert(i<=1);
      
      return mspecies[i];
    }

    
    /** Return number of datasets in this record.
     */
    Index DatasetCount() const { return mdata.nelem(); }

    
    /** Return frequency grid for given dataset.
     */
    ConstVectorView FrequencyGrid(Index dataset) const
    {
        if (dataset > mdata.nelem())
            throw runtime_error("Invalid dataset index.");

        return mdata[dataset].get_numeric_grid(0);
    }
    
    
    /** Return temperatur grid for given dataset.
     */
    ConstVectorView TemperatureGrid(Index dataset) const
    {
        if (dataset > mdata.nelem())
            throw runtime_error("Invalid dataset index.");

        return mdata[dataset].get_numeric_grid(1);
    }


    /** Return CIA data.
     */
    const GriddedField2& Dataset(Index dataset) const
    {
        if (dataset > mdata.nelem())
            throw runtime_error("Invalid dataset index.");

        return mdata[dataset];
    }


    /** Set CIA species.
     \param[in] first CIA Species
     \param[in] second CIA Species
     */
    void SetSpecies(const Index first, const Index second)
    {
        mspecies[0] = first;
        mspecies[1] = second;
    }


    /** Vector version of extract.

     Check whether there is a suitable dataset in the CIARecord and do the 
     interpolation.
     
     \param[out] result CIA value for given frequency grid and temperature.
     \param[in] f_grid Frequency grid
     \param[in] temperature Scalar temparature
     \param[in] dataset Index of dataset to use */
    void Extract(VectorView      result,
                 ConstVectorView f_grid,
                 const Numeric&  temperature,
                 const Index& dataset) const;


    /** Scalar version of extract.
     
     Use the vector version, if you can, it is more efficient. This is just a 
     convenience wrapper for it.
     
     \return Scalar CIA value at given frequency and temperature.
     \param[in] frequency Scalar frequency
     \param[in] temperature Scalar temparature
     \param[in] dataset Index of dataset to use */
    Numeric Extract(const Numeric& frequency,
                    const Numeric& temperature,
                    const Index& dataset) const
    {
      Vector result(1);
      const Vector freqvec(1, frequency);
      
      Extract(result, freqvec, temperature, dataset);
      
      return result[0];
    }


    /** Read CIA catalog file. */
    void ReadFromCIA(const String& filename,
                     const Verbosity& verbosity);
    
private:

    /** Append dataset to mdata. */
    void AppendDataset(const Vector& freq,
                       const ArrayOfNumeric& temp,
                       const ArrayOfVector& cia);


    /** The data itself, directly from the HITRAN file. 
     
     Dimensions:
     Array dimension: Dataset. One file (one molecule pair) can have
                            different datasets, typically for different temperature
                            or frequency ranges.
     Gridded field dimension 1: Frequency [Hz].
     Gridded field dimension 2: Temperature [K].
     Data: FIXME: Document units!!!
     
     */
    ArrayOfGriddedField2 mdata;

    /** The pair of molecules associated with these CIA data.
     
     Molecules are specified by their ARTS internal mspecies index! (This has
     to be determined upon reading from a file. Should it ever be written out, it
     has to be mapped to a string again.)
     
     We use a plain C array here, since the length of this is always 2.
     */
    Index mspecies[2];
};


ostream& operator<<(ostream& os, const CIARecord& cr);

#endif // cia_h
