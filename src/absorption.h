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

#include <stdexcept>
#include <cmath>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "messages.h"
#include "abs_species_tags.h"
#include "linerecord.h"
#include "gridded_fields.h"
#include "jacobian.h"


/** The type that is used to store pointers to lineshape
    functions.  */
typedef void (*lsf_type)(Vector&,//attenuation
                         Vector&,//phase
                         Vector&,//partial attenuation over frequency term
                         Vector&,//partial phase over frequency term
                         Vector&,//partial attenuation over pressure term
                         Vector&,//partial phase over pressure term
                         Vector&,//xvector helper (?)
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         const Numeric,
                         ConstVectorView,
                         const bool,//do phase
                         const bool //do partial derivative
                        );
/**  The types that are used to store pointers to lineshape
     functions internal derivatives.  */
typedef void (*lsf_type_dT)(Vector&,      // dx/dT
                            Numeric&,     // dy/dT
                            Numeric&,     // dFu/dT
                            ConstVectorView,   // frequency
                            const Numeric&,    // line center + line shifts
                            const Numeric&,    // sigma
                            const Numeric&,    // derivative of pressure shift with temperature
                            const Numeric&,    // derivative of line mixing shift with temperature
                            const Numeric&,    // derivative of sigma with temperature
                            const Numeric&,    // gamma
                            const Numeric&     // derivative of gamma with temperature
);
typedef void (*lsf_type_dF)(Numeric&,      // dx/dF
                            const Numeric& // sigma
);
typedef void (*lsf_type_dF0)(Numeric&,      // dx/dF
                             const Numeric& // sigma
);
typedef void (*lsf_type_dgamma)(Numeric&,      // dy/dgamma
                                const Numeric& // sigma
);
typedef void (*lsf_type_dH)(Numeric&,       // dx/dH
                            const Numeric&, // sigma
                            const Numeric&  // derivative of magnetic shift to magnetic strength
);
typedef void (*lsf_type_dDF)(Numeric&,      // dx/dDF
                             const Numeric& // sigma
);

/** Lineshape related information. There is one LineshapeRecord for
    each available lineshape function.

    \author Stefan Buehler
    \date   2000-08-21  */
class LineshapeRecord{
public:

  /** Default constructor. */
  LineshapeRecord() : mname(),
                      mdescription(),
                      mphase(),
                      mpartials(),
                      mfunction(),
                      mfunction_dT(),
                      mfunction_dF(),
                      mfunction_dgamma(),
                      mfunction_dH(),
                      mfunction_dDF()
  { /* Nothing to do here. */ }

  /** Initializing constructor, used to build the lookup table. No partials. */
  LineshapeRecord(const String& name,
                  const String& description,
                  lsf_type      function,
                  const bool    phase,
                  const bool    partials)
    : mname(name),
      mdescription(description),
      mphase(phase),
      mpartials(partials),
      mfunction(function)
  { if(mpartials) throw std::runtime_error("Warning to developers:\n"
      "The line shape is poorly designed since partials" 
      "are expected but none given in the definition.\nThis will fail.\n"); }
  /** Initializing constructor, used to build the lookup table. With partials. */
  LineshapeRecord(const String&   name,
                  const String&   description,
                  lsf_type        function,
                  lsf_type_dT     function_dT,
                  lsf_type_dF     function_dF,
                  lsf_type_dF0    function_dF0,
                  lsf_type_dgamma function_dgamma,
                  lsf_type_dH     function_dH,
                  lsf_type_dDF    function_dDF,
                  const bool      phase,
                  const bool      partials)
  : mname(name),
  mdescription(description),
  mphase(phase),
  mpartials(partials),
  mfunction(function),
  mfunction_dT(function_dT),
  mfunction_dF(function_dF),
  mfunction_dF0(function_dF0),
  mfunction_dgamma(function_dgamma),
  mfunction_dH(function_dH),
  mfunction_dDF(function_dDF)
  { /* Nothing to do here. */ }
  
  /** Return the name of this lineshape. */
  const String&  Name()        const { return mname;        }   
  /** Return the description text. */
  const String&  Description() const { return mdescription; }
  /** Return pointer to lineshape function. */
  lsf_type Function() const { return mfunction; }
  lsf_type_dT dInput_dT() const { return mfunction_dT; }
  lsf_type_dF dInput_dF() const { return mfunction_dF; }
  lsf_type_dF dInput_dF0() const { return mfunction_dF0; }
  lsf_type_dgamma dInput_dgamma() const { return mfunction_dgamma; }
  lsf_type_dH dInput_dH() const { return mfunction_dH; }
  lsf_type_dDF dInput_dDF() const { return mfunction_dDF; }
  /** Returns true if lineshape function calculates phase information. */
  bool Phase() const { return mphase; }
  bool Partials() const { return mpartials; }
private:        
  String  mname;        ///< Name of the function (e.g., Lorentz).
  String  mdescription; ///< Short description.
  bool    mphase;       ///< Does this lineshape calculate phase information?
  bool    mpartials;    ///< Does this lineshape calculate partial derivatives?
  lsf_type mfunction;   ///< Pointer to lineshape function.
  lsf_type_dT mfunction_dT;   ///< Pointer to lineshape function derivative.
  lsf_type_dF mfunction_dF;   ///< Pointer to lineshape function derivative.
  lsf_type_dF0 mfunction_dF0;   ///< Pointer to lineshape function derivative.
  lsf_type_dgamma mfunction_dgamma;   ///< Pointer to lineshape function derivative.
  lsf_type_dH mfunction_dH;   ///< Pointer to lineshape function derivative.
  lsf_type_dDF mfunction_dDF;   ///< Pointer to lineshape function derivative.

};

/** The type that is used to store pointers to lineshape
    normalization functions.  */
typedef void (*lsnf_type)(Vector&,
                          const Numeric,
                          ConstVectorView,
                          const Numeric);

/** Lineshape related normalization function information. There is one
    LineshapeNormRecord for each available lineshape normalization
    function.

    \author Axel von Engeln
    \date   2000-11-30  */
class LineshapeNormRecord{
public:

  /** Default constructor. */
  LineshapeNormRecord() : mname(),
                          mdescription(),
                          mfunction(),
                          mfunction_dT(),
                          mfunction_dF(),
                          mfunction_dF0()
  { /* Nothing to do here. */ }

  /** Initializing constructor, used to build the lookup table. */
  LineshapeNormRecord(const String& name,
                      const String& description,
                      lsnf_type      function,
                      lsnf_type      dfunction_dT,
                      lsnf_type      dfunction_dF,
                      lsnf_type      dfunction_dF0
                     )
    : mname(name),
      mdescription(description),
      mfunction(function),
      mfunction_dT(dfunction_dT),
      mfunction_dF(dfunction_dF),
      mfunction_dF0(dfunction_dF0)
  { /* Nothing to do here. */ }
  /** Return the name of this lineshape. */
  const String&  Name()        const { return mname;        }   
  /** Return the description text. */
  const String&  Description() const { return mdescription; }
  /** Return pointer to lineshape normalization function. */
  lsnf_type Function() const { return mfunction; }
  
  lsnf_type dFunction_dT()  const { return mfunction_dT; }
  lsnf_type dFunction_dF()  const { return mfunction_dF; }
  lsnf_type dFunction_dF0() const { return mfunction_dF0; }
private:        
  String  mname;        ///< Name of the function (e.g., linear).
  String  mdescription; ///< Short description.
  lsnf_type mfunction;  ///< Pointer to lineshape normalization function.
  lsnf_type mfunction_dT;  ///< Pointer to lineshape normalization function partial derivative temperature.
  lsnf_type mfunction_dF;  ///< Pointer to lineshape normalization function partial derivative frequency.
  lsnf_type mfunction_dF0;  ///< Pointer to lineshape normalization function partial derivative frequency.
};

/** Lineshape related specification like which lineshape to use, the
normalizationfactor, and the cutoff.

    \author Axel von Engeln
    \date   2001-01-05  */
class LineshapeSpec{
public:

  /** Default constructor. */
  LineshapeSpec() : mind_ls(-1),
                    mind_lsn(-1),
                    mcutoff(0.)
  { /* Nothing to do here. */ }

  /** Initializing constructor. */
  LineshapeSpec(const Index&    ind_ls,
                const Index&    ind_lsn,
                const Numeric&   cutoff)
    : mind_ls(ind_ls),
      mind_lsn(ind_lsn),
      mcutoff(cutoff)
  { /* Nothing to do here. */ }

  /** Return the index of this lineshape. */
  const Index&  Ind_ls()        const { return mind_ls; }   
  /** Set it. */
  void SetInd_ls( Index ind_ls ) { mind_ls = ind_ls; }

  /** Return the index of the normalization factor. */
  const Index&  Ind_lsn()       const { return mind_lsn; }
  /** Set it. */
  void SetInd_lsn( Index ind_lsn ) { mind_lsn = ind_lsn; }

  /** Return the cutoff frequency (in Hz). This is the distance from
      the line center outside of which the lineshape is defined to be
      zero. Negative means no cutoff.*/
  const Numeric& Cutoff() const { return mcutoff; }
  /** Set it. */
  void SetCutoff( Numeric cutoff ) { mcutoff = cutoff; }
private:        
  Index  mind_ls;
  Index  mind_lsn;
  Numeric mcutoff;
};

ostream& operator<< (ostream& os, const LineshapeSpec& lsspec);

/** Holds a list of lineshape specifications: function, normalization, cutoff.
    \author Axel von Engeln */
typedef Array<LineshapeSpec> ArrayOfLineshapeSpec;



/** Contains the lookup data for one isotopologue.
    \author Stefan Buehler */
class IsotopologueRecord{
public:

  /** Default constructor. Needed by make_array. */
  IsotopologueRecord() = default;
  IsotopologueRecord(const IsotopologueRecord&) = default;
  IsotopologueRecord(IsotopologueRecord&&) = default;
  IsotopologueRecord& operator=(const IsotopologueRecord&) = default;
  IsotopologueRecord& operator=(IsotopologueRecord&&) = default;

  /** Constructor that sets the values. */
  IsotopologueRecord(const String&  name,
                     const Numeric&      abundance,
                     const Numeric&      mass,
                     const Index&        mytrantag,
                     const Index&        hitrantag,
                     const ArrayOfIndex& jpltags) :
    mname(name),
    mabundance(abundance),
    mmass(mass),
    mmytrantag(mytrantag),
    mhitrantag(hitrantag),
    mjpltags(jpltags),
    mqcoeff(),
    mqcoefftype(PF_NOTHING),
    mqcoeffgrid(),
    mqcoeffinterporder()
  {
    // With Matpack, initialization of mjpltags from jpltags should now work correctly.

    // Some consistency checks whether the given data makes sense.
#ifndef NDEBUG
      {
        /* 1. All the tags must be positive or -1 */
        assert( (0<mmytrantag) || (-1==mmytrantag) );
        assert( (0<mhitrantag) || (-1==mhitrantag) );
        for ( Index i=0; i<mjpltags.nelem(); ++i )
          assert( (0<mjpltags[i]) || (-1==mjpltags[i]) );
      }
#endif // ifndef NDEBUG
  }

  /** Isotopologue name. */
  const String&       Name()         const { return mname;  }
  /** Normal abundance ( = isotopologue ratio). (Absolute number.) */
  const Numeric&      Abundance()    const { return mabundance; }
  /** Mass of the isotopologue. (In unified atomic mass units u)
      If I understand this correctly this is the same as g/mol. */
  const Numeric&      Mass()         const { return mmass;    }
  /** MYTRAN2 tag numbers for all isotopologues. -1 means not included. */
  const Index&          MytranTag()    const { return mmytrantag;    }
  /** HITRAN-96 tag numbers for all isotopologues. -1 means not included. */
  const Index&          HitranTag()    const { return mhitrantag;    }
  /** JPL tag numbers for all isotopologues. Empty array means not included. There
      can be more than one JPL tag for an isotopologue species, because in
      JPL different vibrational states have different tags. */
  const ArrayOfIndex&   JplTags()      const { return mjpltags;      }

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
    
  void SetPartitionFctCoeff( const ArrayOfNumeric& qcoeff, const ArrayOfNumeric& temp_range,const Index& qcoefftype )
  {
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
  Numeric CalculatePartitionFctRatio( Numeric reference_temperature,
                                      Numeric actual_temperature ) const
  {
      
      Numeric qcoeff_at_t_ref, qtemp;
      
      switch(mqcoefftype)
        {
          case PF_FROMCOEFF:
            qcoeff_at_t_ref =
            CalculatePartitionFctAtTempFromCoeff_OldUnused( reference_temperature );
            qtemp =
            CalculatePartitionFctAtTempFromCoeff_OldUnused( actual_temperature    );
            break;
          case PF_FROMTEMP:
              qcoeff_at_t_ref =
              CalculatePartitionFctAtTempFromData_OldUnused( reference_temperature    );
              qtemp =
              CalculatePartitionFctAtTempFromData_OldUnused( actual_temperature    );
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
        if ( qtemp > 0. ) 
            return qcoeff_at_t_ref / qtemp;
        else
          {
            std::ostringstream os;
            os << "Partition function of "
               << "Isotopologue " << mname
//               << " is unknown.";
               << " at T=" << actual_temperature << "K is zero or negative.";
            throw std::runtime_error(os.str());
          }
  }
  
  
  enum {
      PF_FROMCOEFF,   // Partition function will be from coefficients
      PF_FROMTEMP,    // Partition function will be from temperature field
      PF_NOTHING      // This will be the designated starter value
  };

private:

  // calculate the partition fct at a certain temperature
  // this is only the prototyping
  Numeric CalculatePartitionFctAtTempFromCoeff_OldUnused( Numeric temperature ) const;
  Numeric CalculatePartitionFctAtTempFromData_OldUnused( Numeric temperature ) const;

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
class SpeciesRecord{
public:

  /** Default constructor. */
  SpeciesRecord() : mname(),
                    mdegfr(-1),
                    misotopologue() { /* Nothing to do here */ }
  
  /** The constructor used in define_species_data. */
  SpeciesRecord(const char name[],
                const Index degfr,
                const Array<IsotopologueRecord>& isotopologue)
    : mname(name),
      mdegfr(degfr),
      misotopologue(isotopologue)
  {

    // Thanks to Matpack, initialization of misotopologue with isotopologue
    // should now work correctly.  

#ifndef NDEBUG
      {
        /* Check that the isotopologues are correctly sorted. */
        for ( Index i=0; i<misotopologue.nelem()-1; ++i )
          {
              assert(std::isnan(misotopologue[i].Abundance()) || std::isnan(misotopologue[i+1].Abundance())
                     || misotopologue[i].Abundance() >= misotopologue[i+1].Abundance());
          }

        /* Check that the Mytran tags are correctly sorted. */
        for ( Index i=0; i<misotopologue.nelem()-1; ++i )
          {
            if ( (0<misotopologue[i].MytranTag()) && (0<misotopologue[i+1].MytranTag()) )
              {
                assert( misotopologue[i].MytranTag() < misotopologue[i+1].MytranTag() );
            
                // Also check that the tags have the same base number:
                assert( misotopologue[i].MytranTag()/10 == misotopologue[i].MytranTag()/10 );
              }
          }

        /* Check that the Hitran tags are correctly sorted. */
        for ( Index i=0; i<misotopologue.nelem()-1; ++i )
          {
            if ( (0<misotopologue[i].HitranTag()) && (0<misotopologue[i+1].HitranTag()) )
              {
//                assert( misotopologue[i].HitranTag() < misotopologue[i+1].HitranTag() );
            
                // Also check that the tags have the same base number:
                assert( misotopologue[i].HitranTag()/10 == misotopologue[i+1].HitranTag()/10 );
              }
          }
      }
#endif // #ifndef NDEBUG
  }

  const String&               Name()     const { return mname;     }   
  Index                         Degfr()    const { return mdegfr;    }
  const Array<IsotopologueRecord>& Isotopologue()  const { return misotopologue;  }
  Array<IsotopologueRecord>&       Isotopologue()        { return misotopologue;  }
  
private:
  /** Species name. */
  String mname;
  /** Degrees of freedom. */
  Index mdegfr;
  /** Isotopologue data. */
  Array<IsotopologueRecord> misotopologue;
};


/** Auxiliary data for isotopologues */
class SpeciesAuxData
{
public:
    typedef enum
    {
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
    SpeciesAuxData() { };

    void InitFromSpeciesData();

    /** Get single parameter value. */
    AuxType getType(Index species, Index isotopologue) const;

    /** Returns number of species. */
    Index nspecies() const { return mparams.nelem(); };

    /** Returns number of isotopologues for a certain species. */
    Index nisotopologues(const Index species) const { return mparams[species].nelem(); };

    /** Return a constant reference to the parameters. */
    const ArrayOfGriddedField1& getParam(const Index species, const Index isotopologue) const;
    
    /** Returns mparams[st.Species()][st.Isotopologue()][0].data[0] if st.Isotopologue() > 0, else 1 */
    Numeric getIsotopologueRatio(const SpeciesTag& st) const;
    
    /** Return a constant reference to the parameters. */
    const ArrayOfGriddedField1& getParam(const LineRecord& lr) const
    { return getParam(lr.Species(), lr.Isotopologue()); }

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
                                const Index isotopologue) const
    { return mparam_type[species][isotopologue]; }

    /** Return a constant reference to the parameter types. */
    const AuxType& getParamType(const LineRecord& lr) const
    { return getParamType(lr.Species(), lr.Isotopologue()); }

    /** Read parameters from input stream (only for version 1 format). */
    bool ReadFromStream(String& artsid, istream& is, Index nparams,
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
void fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(SpeciesAuxData& sad);

/** Fill SpeciesAuxData with default partition functions from species data. */
void fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(SpeciesAuxData& sad);


// is needed to map jpl tags/arts identifier to the species/isotopologue data within arts
class SpecIsoMap{
public:
  SpecIsoMap():mspeciesindex(0), misotopologueindex(0){}
  SpecIsoMap(const Index& speciesindex,
                const Index& isotopologueindex)
    : mspeciesindex(speciesindex),
      misotopologueindex(isotopologueindex) 
  {}

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
ostream& operator<< (ostream& os, const SpeciesRecord& sr);

/** Output operator for SpeciesAuxData.
    \author Oliver Lemke */
ostream& operator<< (ostream& os, const SpeciesAuxData& sad);



void xsec_species(    MatrixView               xsec_attenuation,
		      MatrixView               xsec_source,
		      MatrixView               xsec_phase,
		      ConstVectorView          f_grid,
		      ConstVectorView          abs_p,
		      ConstVectorView          abs_t,
                      ConstMatrixView          abs_t_nlte,
		      ConstMatrixView          all_vmrs,
		      const ArrayOfArrayOfSpeciesTag& abs_species,
		      const ArrayOfLineRecord& abs_lines,
		      const Index              ind_ls,
		      const Index              ind_lsn,
		      const Numeric            cutoff,
		      const SpeciesAuxData&    isotopologue_ratios,
                      const SpeciesAuxData&    partition_functions,
		      const Verbosity&         verbosity );


void xsec_single_line(// Output:
                      VectorView xsec_accum_attenuation, 
                      VectorView xsec_accum_source, 
                      VectorView xsec_accum_phase, 
                      // Helper variables
                      Vector& attenuation, 
                      Vector& phase,
                      Vector& da_dF,
                      Vector& dp_dF,
                      Vector& da_dP,
                      Vector& dp_dP,
                      Range&  this_f_range,
                      //Vector& dC_dT, for normalization factor
                      Vector& fac, 
                      Vector& aux, 
                      // Frequency grid:
                      Vector& f_local, 
                      const ConstVectorView f_grid, 
                      const Index nf, 
                      const Numeric cutoff,
                      Numeric F0, 
                      // Line strength:
                      Numeric intensity, 
                      const Numeric part_fct_ratio,  
                      const Numeric boltzmann_ratio,
                      const Numeric abs_nlte_ratio,
                      const Numeric src_nlte_ratio,
                      const Numeric Isotopologue_Ratio,
                      // Atmospheric state
                      const Numeric temperature, 
                      // Line shape:
                      const Index ind_ls, 
                      const Index ind_lsn,  
                      // Line broadening:
                      const Numeric gamma_0,
                      const Numeric gamma_2,
                      const Numeric eta,
                      const Numeric df_0,
                      const Numeric df_2,
                      const Numeric sigma,
                      const Numeric f_VC,
                      // Line mixing
                      const Numeric LM_DF,
                      const Numeric LM_Y, 
                      const Numeric LM_G,
                      // Feature flags
                      const bool calc_cut, 
                      const bool calc_phase,
                      const bool calc_partials,
                      const bool calc_src);


void xsec_species_line_mixing_wrapper(      MatrixView               xsec_attenuation,
					    MatrixView               xsec_source,
                                            MatrixView               xsec_phase,
                                            ArrayOfMatrix&           partial_xsec_attenuation,
                                            ArrayOfMatrix&           partial_xsec_source,
                                            ArrayOfMatrix&           partial_xsec_phase,
                                            const ArrayOfRetrievalQuantity& flag_partials,
                                            const ArrayOfIndex& flag_partials_position,
                                            ConstVectorView          f_grid,
                                            ConstVectorView          abs_p,
                                            ConstVectorView          abs_t,
                                            ConstMatrixView          abs_t_nlte,
                                            ConstMatrixView          all_vmrs,
                                            const ArrayOfArrayOfSpeciesTag& abs_species,
                                            const Index              this_species,
                                            const ArrayOfLineRecord& abs_lines,
                                            const Numeric            H_magnitude_Zeeman,
                                            const Index              ind_ls,
                                            const Index              ind_lsn,
                                            const Numeric            cutoff,
                                            const SpeciesAuxData&    isotopologue_ratios,
                                            const SpeciesAuxData&    partition_functions);


void calc_gamma_and_deltaf_artscat4(Numeric& gamma,
                                    Numeric& deltaf,
                                    const Numeric p,
                                    const Numeric t,
                                    ConstVectorView vmrs,
                                    const Index this_species,
                                    const ArrayOfIndex& broad_spec_locations,
                                    const Numeric& T0,
                                    const Numeric& Sgam,
                                    const Numeric& Nself,
                                    const Vector&  Gamma_foreign,
                                    const Vector&  N_foreign,
                                    const Vector&  Delta_foreign,
                                    const Verbosity& verbosity);

void calc_gamma_and_deltaf_artscat4_old_unused(Numeric& gamma,
                                               Numeric& deltaf,
                                               const Numeric p,
                                               const Numeric t,
                                               ConstVectorView vmrs,
                                               const Index this_species,
                                               const ArrayOfIndex& broad_spec_locations,
                                               const LineRecord& l_l,
                                               const Verbosity& verbosity);


// A helper function for energy conversion:
Numeric wavenumber_to_joule(Numeric e);


//======================================================================
//             Functions related to species
//======================================================================

Index species_index_from_species_name( String name );

String species_name_from_species_index( const Index spec_ind );

//======================================================================
//             Functions to convert the accuracy index
//======================================================================

// ********* for HITRAN database *************
// convert index for the frequency accuracy.
void convHitranIERF(     
                    Numeric&     mdf,
              const Index&       df 
                    );

// convert to percents index for intensity and halfwidth accuracy.

void convHitranIERSH(     
                    Numeric&     mdh,
              const Index&       dh 
                    );

// ********* for MYTRAN database *************
// convert index for the halfwidth accuracy.
void convMytranIER(     
                    Numeric&     mdh,
              const Index  &      dh 
                    );


// Functions to set abs_n2 and abs_h2o:

void abs_n2Set(Vector&            abs_n2,
               const ArrayOfArrayOfSpeciesTag& abs_species,
               const Matrix&    abs_vmrs,
               const Verbosity&);

void abs_h2oSet(Vector&          abs_h2o,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Matrix&    abs_vmrs,
                const Verbosity&);

void xsec_species2(Matrix& xsec, 
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
                   const Matrix& abs_t_nlte,
                   const Matrix& all_vmrs,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const ArrayOfLineRecord& abs_lines,
                   const SpeciesAuxData& isotopologue_ratios,
                   const SpeciesAuxData& partition_functions);

#endif // absorption_h
