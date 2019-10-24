/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>
   Axel von Engeln <engeln@uni-bremen.de>

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
    Declarations required for the calculation of absorption coefficients.

    This is the file from arts-1-0, back-ported to arts-1-1.

    \author Stefan Buehler, Axel von Engeln
*/

#ifndef absorption_h
#define absorption_h

#include <cmath>
#include <stdexcept>
#include "abs_species_tags.h"
#include "array.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "matpackI.h"
#include "messages.h"
#include "mystring.h"
#include "absorptionlines.h"

/** Contains the lookup data for one isotopologue.
    \author Stefan Buehler */
class IsotopologueRecord {
 public:
  /** Default constructor. Needed by make_array. */
  IsotopologueRecord() = default;
  IsotopologueRecord(const IsotopologueRecord&) = default;
  IsotopologueRecord(IsotopologueRecord&&) = default;
  IsotopologueRecord& operator=(const IsotopologueRecord&) = default;
  IsotopologueRecord& operator=(IsotopologueRecord&&) = default;

  /** Constructor that sets the values. */
  IsotopologueRecord(const String& name,
                     const Numeric& abundance,
                     const Numeric& mass,
                     const Index& mytrantag,
                     const Index& hitrantag,
                     const ArrayOfIndex& jpltags)
      : mname(name),
        mabundance(abundance),
        mmass(mass),
        mmytrantag(mytrantag),
        mhitrantag(hitrantag),
        mjpltags(jpltags),
        mqcoeff(),
        mqcoefftype(PF_NOTHING),
        mqcoeffgrid(),
        mqcoeffinterporder() {
    // With Matpack, initialization of mjpltags from jpltags should now work correctly.

    // Some consistency checks whether the given data makes sense.
#ifndef NDEBUG
    {
      /* 1. All the tags must be positive or -1 */
      assert((0 < mmytrantag) || (-1 == mmytrantag));
      assert((0 < mhitrantag) || (-1 == mhitrantag));
      for (Index i = 0; i < mjpltags.nelem(); ++i)
        assert((0 < mjpltags[i]) || (-1 == mjpltags[i]));
    }
#endif  // ifndef NDEBUG
  }

  /** Isotopologue name. */
  const String& Name() const { return mname; }
  /** Normal abundance ( = isotopologue ratio). (Absolute number.) */
  const Numeric& Abundance() const { return mabundance; }
  /** Mass of the isotopologue. (In unified atomic mass units u)
      If I understand this correctly this is the same as g/mol. */
  const Numeric& Mass() const { return mmass; }
  /** MYTRAN2 tag numbers for all isotopologues. -1 means not included. */
  const Index& MytranTag() const { return mmytrantag; }
  /** HITRAN-96 tag numbers for all isotopologues. -1 means not included. */
  const Index& HitranTag() const { return mhitrantag; }
  /** JPL tag numbers for all isotopologues. Empty array means not included. There
      can be more than one JPL tag for an isotopologue species, because in
      JPL different vibrational states have different tags. */
  const ArrayOfIndex& JplTags() const { return mjpltags; }

  //! Check if isotopologue is actually a continuum.
  /*!
   \return True if this is a continuum.
   */
  bool isContinuum() const { return mname.length() && !isdigit(mname[0]); }

  //! Return the partition function coefficients.
  const Vector& GetCoeff() const { return mqcoeff; }

  //! Return the partition function coefficients.
  const Vector& GetCoeffGrid() const { return mqcoeffgrid; }

  //! Return the partition function coefficient types.
  Index GetCoeffType() const { return mqcoefftype; }

  void SetPartitionFctCoeff(const ArrayOfNumeric& qcoeff,
                            const ArrayOfNumeric& temp_range,
                            const Index& qcoefftype) {
    mqcoeff = qcoeff;
    mqcoeffgrid = temp_range;
    mqcoefftype = qcoefftype;
  }

  //! Calculate partition function ratio.
  /*!
    This computes the partition function ratio Q(Tref)/Q(T). 

    Unfortunately, we have to recalculate also Q(Tref) for each
    spectral line, because the reference temperatures can be
    different!
    
    \param reference_temperature The reference temperature.
    \param actual_temperature The actual temperature.
  
    \return The ratio.
  */
  Numeric CalculatePartitionFctRatio(Numeric reference_temperature,
                                     Numeric actual_temperature) const {
    Numeric qcoeff_at_t_ref, qtemp;

    switch (mqcoefftype) {
      case PF_FROMCOEFF:
        qcoeff_at_t_ref = CalculatePartitionFctAtTempFromCoeff_OldUnused(
            reference_temperature);
        qtemp =
            CalculatePartitionFctAtTempFromCoeff_OldUnused(actual_temperature);
        break;
      case PF_FROMTEMP:
        qcoeff_at_t_ref = CalculatePartitionFctAtTempFromData_OldUnused(
            reference_temperature);
        qtemp =
            CalculatePartitionFctAtTempFromData_OldUnused(actual_temperature);
        break;
      default:
        throw runtime_error("The partition functions are incorrect.\n");
        break;
    }
    /*        cout << "ref_t: " << reference_temperature << ", act_t:" <<
          actual_temperature << "\n";
        cout << "ref_q: " << qcoeff_at_t_ref << ", act_q:" <<
          qtemp << "\n";
*/
    if (qtemp > 0.)
      return qcoeff_at_t_ref / qtemp;
    else {
      std::ostringstream os;
      os << "Partition function of "
         << "Isotopologue "
         << mname
         //               << " is unknown.";
         << " at T=" << actual_temperature << "K is zero or negative.";
      throw std::runtime_error(os.str());
    }
  }

  enum {
    PF_FROMCOEFF,  // Partition function will be from coefficients
    PF_FROMTEMP,   // Partition function will be from temperature field
    PF_NOTHING     // This will be the designated starter value
  };

 private:
  // calculate the partition fct at a certain temperature
  // this is only the prototyping
  Numeric CalculatePartitionFctAtTempFromCoeff_OldUnused(
      Numeric temperature) const;
  Numeric CalculatePartitionFctAtTempFromData_OldUnused(
      Numeric temperature) const;

  String mname;
  Numeric mabundance;
  Numeric mmass;
  Index mmytrantag;
  Index mhitrantag;
  ArrayOfIndex mjpltags;
  Vector mqcoeff;
  Index mqcoefftype;
  Vector mqcoeffgrid;
  Index mqcoeffinterporder;
};

/** Contains the lookup data for one species.

    \author Stefan Buehler  */
class SpeciesRecord {
 public:
  /** Default constructor. */
  SpeciesRecord()
      : mname(), mdegfr(-1), misotopologue() { /* Nothing to do here */
  }

  /** The constructor used in define_species_data. */
  SpeciesRecord(const char name[],
                const Index degfr,
                const Array<IsotopologueRecord>& isotopologue)
      : mname(name), mdegfr(degfr), misotopologue(isotopologue) {
    // Thanks to Matpack, initialization of misotopologue with isotopologue
    // should now work correctly.

#ifndef NDEBUG
    {
      /* Check that the isotopologues are correctly sorted. */
      for (Index i = 0; i < misotopologue.nelem() - 1; ++i) {
        assert(std::isnan(misotopologue[i].Abundance()) ||
               std::isnan(misotopologue[i + 1].Abundance()) ||
               misotopologue[i].Abundance() >=
                   misotopologue[i + 1].Abundance());
      }

      /* Check that the Mytran tags are correctly sorted. */
      for (Index i = 0; i < misotopologue.nelem() - 1; ++i) {
        if ((0 < misotopologue[i].MytranTag()) &&
            (0 < misotopologue[i + 1].MytranTag())) {
          assert(misotopologue[i].MytranTag() <
                 misotopologue[i + 1].MytranTag());

          // Also check that the tags have the same base number:
          assert(misotopologue[i].MytranTag() / 10 ==
                 misotopologue[i].MytranTag() / 10);
        }
      }

      /* Check that the Hitran tags are correctly sorted. */
      for (Index i = 0; i < misotopologue.nelem() - 1; ++i) {
        if ((0 < misotopologue[i].HitranTag()) &&
            (0 < misotopologue[i + 1].HitranTag())) {
          //                assert( misotopologue[i].HitranTag() < misotopologue[i+1].HitranTag() );

          // Also check that the tags have the same base number:
          assert(misotopologue[i].HitranTag() / 10 ==
                 misotopologue[i + 1].HitranTag() / 10);
        }
      }
    }
#endif  // #ifndef NDEBUG
  }

  const String& Name() const { return mname; }
  Index Degfr() const { return mdegfr; }
  const Array<IsotopologueRecord>& Isotopologue() const {
    return misotopologue;
  }
  Array<IsotopologueRecord>& Isotopologue() { return misotopologue; }

 private:
  /** Species name. */
  String mname;
  /** Degrees of freedom. */
  Index mdegfr;
  /** Isotopologue data. */
  Array<IsotopologueRecord> misotopologue;
};

/** Auxiliary data for isotopologues */
class SpeciesAuxData {
 public:
  typedef enum {
    AT_NONE,
    AT_ISOTOPOLOGUE_RATIO,
    AT_ISOTOPOLOGUE_QUANTUM,
    AT_PARTITIONFUNCTION_TFIELD,
    AT_PARTITIONFUNCTION_COEFF,
    AT_PARTITIONFUNCTION_COEFF_VIBROT,
    AT_FINAL_ENTRY
  } AuxType;

  typedef Array<Array<AuxType> > ArrayOfArrayOfAuxType;
  typedef Array<GriddedField1> AuxData;
  typedef Array<Array<AuxData> > ArrayOfArrayOfAuxData;

  /** Default constructor. */
  SpeciesAuxData(){};

  void InitFromSpeciesData();

  /** Get single parameter value. */
  AuxType getType(Index species, Index isotopologue) const;

  /** Returns number of species. */
  Index nspecies() const { return mparams.nelem(); };

  /** Returns number of isotopologues for a certain species. */
  Index nisotopologues(const Index species) const {
    return mparams[species].nelem();
  };

  /** Return a constant reference to the parameters. */
  const ArrayOfGriddedField1& getParam(const Index species,
                                       const Index isotopologue) const;

  /** Returns mparams[st.Species()][st.Isotopologue()][0].data[0] if st.Isotopologue() > 0, else 1 */
  Numeric getIsotopologueRatio(const SpeciesTag& st) const;
  
  /** Returns mparams[qid.Species()][qid.Isotopologue()][0].data[0] */
  Numeric getIsotopologueRatio(const QuantumIdentifier& qid) const;
  
  /** Return a constant reference to the parameters. */
  const ArrayOfGriddedField1& getParam(const QuantumIdentifier& qid) const {
    return getParam(qid.Species(), qid.Isotopologue());
  }

  /** Return a parameter type as string. */
  String getTypeString(const Index species, const Index isotopologue) const;

  /** Set parameter. */
  void setParam(const Index species,
                const Index isotopologue,
                const AuxType auxtype,
                const ArrayOfGriddedField1& auxdata);

  /** Set parameter by ARTS tag. */
  void setParam(const String& artstag,
                const String& auxtype,
                const ArrayOfGriddedField1& auxdata);

  /** Return a constant reference to the parameter types. */
  const AuxType& getParamType(const Index species,
                              const Index isotopologue) const {
    return mparam_type[species][isotopologue];
  }

  /** Return a constant reference to the parameter types. */
  const AuxType& getParamType(const QuantumIdentifier& qid) const {
    return getParamType(qid.Species(), qid.Isotopologue());
  }

  /** Read parameters from input stream (only for version 1 format). */
  bool ReadFromStream(String& artsid,
                      istream& is,
                      Index nparams,
                      const Verbosity& verbosity);

 private:
  ArrayOfArrayOfAuxData mparams;
  ArrayOfArrayOfAuxType mparam_type;
};

/** Check that isotopologue ratios for the given species are correctly defined. */
void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& sad);

/** Check that partition functions for the given species are correctly defined. */
void checkPartitionFunctions(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const SpeciesAuxData& partfun);

/** Fill SpeciesAuxData with default isotopologue ratios from species data. */
void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(
    SpeciesAuxData& sad);

/** Fill SpeciesAuxData with default partition functions from species data. */
void fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(
    SpeciesAuxData& sad);

// is needed to map jpl tags/arts identifier to the species/isotopologue data within arts
class SpecIsoMap {
 public:
  SpecIsoMap() : mspeciesindex(0), misotopologueindex(0) {}
  SpecIsoMap(const Index& speciesindex, const Index& isotopologueindex)
      : mspeciesindex(speciesindex), misotopologueindex(isotopologueindex) {}

  // Return the index to the species
  const Index& Speciesindex() const { return mspeciesindex; }
  // Return the index to the isotopologue
  const Index& Isotopologueindex() const { return misotopologueindex; }

 private:
  Index mspeciesindex;
  Index misotopologueindex;
};

/** Output operator for SpeciesRecord. Incomplete version: only writes
    SpeciesName.

    \author Jana Mendrok */
ostream& operator<<(ostream& os, const SpeciesRecord& sr);

/** Output operator for SpeciesAuxData.
    \author Oliver Lemke */
ostream& operator<<(ostream& os, const SpeciesAuxData& sad);

// A helper function for energy conversion:
Numeric wavenumber_to_joule(Numeric e);

//======================================================================
//             Functions related to species
//======================================================================

Index species_index_from_species_name(String name);

String species_name_from_species_index(const Index spec_ind);

//======================================================================
//             Functions to convert the accuracy index
//======================================================================

// ********* for HITRAN database *************
// convert index for the frequency accuracy.
void convHitranIERF(Numeric& mdf, const Index& df);

// convert to percents index for intensity and halfwidth accuracy.

void convHitranIERSH(Numeric& mdh, const Index& dh);

// ********* for MYTRAN database *************
// convert index for the halfwidth accuracy.
void convMytranIER(Numeric& mdh, const Index& dh);

// Functions to set abs_n2 and abs_h2o:

void abs_n2Set(Vector& abs_n2,
               const ArrayOfArrayOfSpeciesTag& abs_species,
               const Matrix& abs_vmrs,
               const Verbosity&);

void abs_h2oSet(Vector& abs_h2o,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Matrix& abs_vmrs,
                const Verbosity&);

void xsec_species(Matrix& xsec,
                  Matrix& source,
                  Matrix& phase,
                  ArrayOfMatrix& dxsec_dx,
                  ArrayOfMatrix& dsource_dx,
                  ArrayOfMatrix& dphase_dx,
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const ArrayOfIndex& jacobian_propmat_positions,
                  const Vector& f_grid,
                  const Vector& abs_p,
                  const Vector& abs_t,
                  const EnergyLevelMap& abs_nlte,
                  const Matrix& all_vmrs,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const AbsorptionLines& lines,
                  const Numeric& isot_ratio,
                  const SpeciesAuxData::AuxType& partfun_type,
                  const ArrayOfGriddedField1& partfun_data);

/** Returns the species data */
const SpeciesRecord& SpeciesDataOfLines(const AbsorptionLines& lines);

#endif  // absorption_h
