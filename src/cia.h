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
#include "gridded_fields.h"
#include "matpackI.h"
#include "mystring.h"
#include "species.h"
#include "messages.h"

// Declare existance of some classes:
class bifstream;
class CIARecord;

typedef Array<CIARecord> ArrayOfCIARecord;

/* Header with implementation. */
void cia_interpolation(VectorView result,
                       ConstVectorView frequency,
                       const Numeric& temperature,
                       const GriddedField2& cia_data,
                       const Numeric& T_extrapolfac,
                       const Index& robust,
                       const Verbosity& verbosity);

Index cia_get_index(const ArrayOfCIARecord& cia_data,
                    const Species::Species sp1,
                    const Species::Species sp2);

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
  void SetMoleculeName(const Index i, const String& name);

  /** Return CIA species index.
     
     \param[in] i Must be either 0 or 1. Then the first or second species index
     is returned.
     */
  Species::Species Species(const Index i) const {
    // Assert that i is 0 or 1:
    ARTS_ASSERT(i >= 0);
    ARTS_ASSERT(i <= 1);

    return mspecies[i];
  }

  /** Return number of datasets in this record.
     */
  Index DatasetCount() const { return mdata.nelem(); }

  /** Return frequency grid for given dataset.
     */
  ConstVectorView FrequencyGrid(Index dataset) const {
    ARTS_ASSERT(dataset >= 0);
    ARTS_ASSERT(dataset < mdata.nelem());

    return mdata[dataset].get_numeric_grid(0);
  }

  /** Return temperatur grid for given dataset.
     */
  ConstVectorView TemperatureGrid(Index dataset) const {
    ARTS_ASSERT(dataset >= 0);
    ARTS_ASSERT(dataset < mdata.nelem());

    return mdata[dataset].get_numeric_grid(1);
  }

  /** Return CIA dataset.
     */
  const GriddedField2& Dataset(Index dataset) const {
    ARTS_ASSERT(dataset >= 0);
    ARTS_ASSERT(dataset < mdata.nelem());

    return mdata[dataset];
  }

  /** Return CIA data.
     */
  const ArrayOfGriddedField2& Data() const { return mdata; }
  
  /** Return CIA data.
   */
  ArrayOfGriddedField2& Data() { return mdata; }

  /** Set CIA species.
     \param[in] first CIA Species.
     \param[in] second CIA Species.
     */
  void SetSpecies(const Species::Species first, const Species::Species second) {
    mspecies[0] = first;
    mspecies[1] = second;
  }

  /** Vector version of extract.

     Check whether there is a suitable dataset in the CIARecord and do the 
     interpolation.
     
     \param[out] result CIA value for given frequency grid and temperature.
     \param[in] f_grid Frequency grid.
     \param[in] temperature Scalar temparature.
     \param[in] dataset Index of dataset to use.
     \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
     \param[in] verbosity   Standard verbosity object.
     */
  void Extract(VectorView result,
               ConstVectorView f_grid,
               const Numeric& temperature,
               const Index& dataset,
               const Numeric& T_extrapolfac,
               const Index& robust,
               const Verbosity& verbosity) const;

  /** Scalar version of extract.
     
     Use the vector version, if you can, it is more efficient. This is just a 
     convenience wrapper for it.
     
     \return Scalar CIA value at given frequency and temperature.
     \param[in] frequency Scalar frequency
     \param[in] temperature Scalar temparature
     \param[in] dataset Index of dataset to use 
     \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
     \param[in] verbosity   Standard verbosity object.
     */
  Numeric Extract(const Numeric& frequency,
                  const Numeric& temperature,
                  const Index& dataset,
                  const Numeric& T_extrapolfac,
                  const Index& robust,
                  const Verbosity& verbosity) const {
    Vector result(1);
    const Vector freqvec(1, frequency);

    Extract(result,
            freqvec,
            temperature,
            dataset,
            T_extrapolfac,
            robust,
            verbosity);

    return result[0];
  }

  /** Read CIA catalog file. */
  void ReadFromCIA(const String& filename, const Verbosity& verbosity);

  friend void xml_read_from_stream(istream& is_xml,
                                   CIARecord& cr,
                                   bifstream* pbifs,
                                   const Verbosity& verbosity);

  /** Append other CIARecord to this. */
  void AppendDataset(const CIARecord& c2);

std::array<Species::Species, 2> TwoSpecies() const {return mspecies;}
std::array<Species::Species, 2>& TwoSpecies() {return mspecies;}

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
     Data: Binary absorption cross-sections in m^5 molec^(-2) 
     
     */
  ArrayOfGriddedField2 mdata;

  /** The pair of molecules associated with these CIA data.
     
     Molecules are specified by their ARTS internal mspecies index! (This has
     to be determined upon reading from a file. Should it ever be written out, it
     has to be mapped to a string again.)
     
     We use a plain C array here, since the length of this is always 2.
     */
  std::array<Species::Species, 2> mspecies;
};

ostream& operator<<(ostream& os, const CIARecord& cr);

#endif  // cia_h
