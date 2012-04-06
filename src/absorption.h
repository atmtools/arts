/* Copyright (C) 2000-2008
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
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "make_array.h"
#include "messages.h"

/** The type that is used to store pointers to lineshape
    functions.  */
typedef void (*lsf_type)(Vector&,
                         Vector&,
                         Numeric,
                         Numeric,
                         Numeric,
                         VectorView,
                         const Index);

/** Lineshape related information. There is one LineshapeRecord for
    each available lineshape function.

    \author Stefan Buehler
    \date   2000-08-21  */
class LineshapeRecord{
public:

  /** Default constructor. */
  LineshapeRecord() : mname(),
                      mdescription(),
                      mfunction()
  { /* Nothing to do here. */ }

  /** Initializing constructor, used to build the lookup table. */
  LineshapeRecord(const String& name,
                  const String& description,
                  lsf_type      function)
    : mname(name),
      mdescription(description),
      mfunction(function)
  { /* Nothing to do here. */ }
  /** Return the name of this lineshape. */
  const String&  Name()        const { return mname;        }   
  /** Return the description text. */
  const String&  Description() const { return mdescription; }
  /** Return pointer to lineshape function. */
  lsf_type Function() const { return mfunction; }
private:        
  String  mname;        ///< Name of the function (e.g., Lorentz).
  String  mdescription; ///< Short description.
  lsf_type mfunction;   ///< Pointer to lineshape function.

};

/** The type that is used to store pointers to lineshape
    normalization functions.  */
typedef void (*lsnf_type)(Vector&,
                          Numeric,
                          VectorView,
                          const Numeric,
                          const Index);

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
                          mfunction()
  { /* Nothing to do here. */ }

  /** Initializing constructor, used to build the lookup table. */
  LineshapeNormRecord(const String& name,
                      const String& description,
                      lsnf_type      function)
    : mname(name),
      mdescription(description),
      mfunction(function)
  { /* Nothing to do here. */ }
  /** Return the name of this lineshape. */
  const String&  Name()        const { return mname;        }   
  /** Return the description text. */
  const String&  Description() const { return mdescription; }
  /** Return pointer to lineshape normalization function. */
  lsnf_type Function() const { return mfunction; }
private:        
  String  mname;        ///< Name of the function (e.g., linear).
  String  mdescription; ///< Short description.
  lsnf_type mfunction;  ///< Pointer to lineshape normalization function.
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



/** Contains the lookup data for one isotope.
    \author Stefan Buehler */
class IsotopeRecord{
public:

  /** Default constructor. Needed by make_array. */
  IsotopeRecord() : mname(),
                    mabundance(0.),
                    mmass(0.),
                    mmytrantag(-1),
                    mhitrantag(-1),
                    mjpltags(),
                    mqcoeff()
  { /* Nothing left to do here. */ }

  /** Copy constructor. We need this, since operator= does not work
      correctly for Arrays. (Target Array has to be resized first.) */
  IsotopeRecord(const IsotopeRecord& x) :
    mname(x.mname),
    mabundance(x.mabundance),
    mmass(x.mmass),
    mmytrantag(x.mmytrantag),
    mhitrantag(x.mhitrantag),
    mjpltags(x.mjpltags),
    mqcoeff()
  { /* Nothing left to do here. */ }

  /** Constructor that sets the values. */
  IsotopeRecord(const String&           name,
                const Numeric&          abundance,
                const Numeric&          mass,
                const Index&            mytrantag,
                const Index&            hitrantag,
                const MakeArray<Index>& jpltags) :
    mname(name),
    mabundance(abundance),
    mmass(mass),
    mmytrantag(mytrantag),
    mhitrantag(hitrantag),
    mjpltags(jpltags),
    mqcoeff()
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

  /** Isotope name. */
  const String&       Name()         const { return mname;  }
  /** Normal abundance ( = isotopic ratio). (Absolute number.) */
  const Numeric&      Abundance()    const { return mabundance; }
  /** Mass of the isotope. (In unified atomic mass units u)
      If I understand this correctly this is the same as g/mol. */
  const Numeric&      Mass()         const { return mmass;    }
  /** MYTRAN2 tag numbers for all isotopes. -1 means not included. */
  const Index&          MytranTag()    const { return mmytrantag;    }
  /** HITRAN-96 tag numbers for all isotopes. -1 means not included. */
  const Index&          HitranTag()    const { return mhitrantag;    }
  /** JPL tag numbers for all isotopes. Empty array means not included. There
      can be more than one JPL tag for an isotopic species, because in
      JPL different vibrational states have different tags. */
  const ArrayOfIndex&   JplTags()      const { return mjpltags;      }

  void SetPartitionFctCoeff( const ArrayOfNumeric& qcoeff )
  {
    mqcoeff = qcoeff;
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
    Numeric qcoeff_at_t_ref =
      CalculatePartitionFctAtTemp( reference_temperature );
    Numeric qtemp =
      CalculatePartitionFctAtTemp( actual_temperature    );

    if ( qtemp < 0. ) 
      {
        ostringstream os;
        os << "Partition function of "
           << "Isotope = " << mname
           << "is unknown.";
        throw runtime_error(os.str());
      }
    return qcoeff_at_t_ref / qtemp;
  }

private:

  // calculate the partition fct at a certain temperature
  // this is only the prototyping
  Numeric CalculatePartitionFctAtTemp( Numeric temperature ) const;

  String mname;
  Numeric mabundance;
  Numeric mmass;
  Index mmytrantag;
  Index mhitrantag;
  ArrayOfIndex mjpltags;
  ArrayOfNumeric mqcoeff;
};


/** Contains the lookup data for one species.

    \author Stefan Buehler  */
class SpeciesRecord{
public:

  /** Default constructor. */
  SpeciesRecord() : mname(),
                    mdegfr(-1),
                    misotope() { /* Nothing to do here */ }
  
  /** The constructor used in define_species_data. */
  SpeciesRecord(const char name[],
                const Index degfr,
                const MakeArray<IsotopeRecord>& isotope)
    : mname(name),
      mdegfr(degfr),
      misotope(isotope)
  {

    // Thanks to Matpack, initialization of misotope with isotope
    // should now work correctly.  

#ifndef NDEBUG
      {
        /* Check that the isotopes are correctly sorted. */
        for ( Index i=0; i<misotope.nelem()-1; ++i )
          {
            assert( misotope[i].Abundance() >= misotope[i+1].Abundance() );
          }

        /* Check that the Mytran tags are correctly sorted. */
        for ( Index i=0; i<misotope.nelem()-1; ++i )
          {
            if ( (0<misotope[i].MytranTag()) && (0<misotope[i+1].MytranTag()) )
              {
                assert( misotope[i].MytranTag() < misotope[i+1].MytranTag() );
            
                // Also check that the tags have the same base number:
                assert( misotope[i].MytranTag()/10 == misotope[i].MytranTag()/10 );
              }
          }

        /* Check that the Hitran tags are correctly sorted. */
        for ( Index i=0; i<misotope.nelem()-1; ++i )
          {
            if ( (0<misotope[i].HitranTag()) && (0<misotope[i+1].HitranTag()) )
              {
                assert( misotope[i].HitranTag() < misotope[i+1].HitranTag() );
            
                // Also check that the tags have the same base number:
                assert( misotope[i].HitranTag()/10 == misotope[i+1].HitranTag()/10 );
              }
          }
      }
#endif // #ifndef NDEBUG
  }

  const String&               Name()     const { return mname;     }   
  Index                         Degfr()    const { return mdegfr;    }
  const Array<IsotopeRecord>& Isotope()  const { return misotope;  }
  Array<IsotopeRecord>&       Isotope()        { return misotope;  }
  
private:
  /** Species name. */
  String mname;
  /** Degrees of freedom. */
  Index mdegfr;
  /** Isotope data. */
  Array<IsotopeRecord> misotope;
};

/** Spectral line catalog data. 

    Below is a description of the ARTS catalogue format.
    The file starts with the usual XML header:

    \verbatim
    <?xml version="1.0"?>
    <arts format="ascii" version="1">
    <ArrayOfLineRecord version="ARTSCAT-3" nelem="8073">
    \endverbatim

    The ARTSCAT version number is there to keep track of catalogue format changes.
    The "nelem" tag contains the total number of lines in the file.

    The file ends with the usual XML closing tags:

    \verbatim
    </ArrayOfLineRecord>
    </arts>
    \endverbatim

    In-between the header and the footer are the actual spectroscopic
    data. Each new entry, corresponding to one spectral line, starts
    with the `@' character.

    The line catalogue should not have any fixed column widths because
    the precision of the parameters should not be limited by the
    format.  The catalogue can then be stored as binary or ASCII. In
    the ASCII version the columns are separated by one or more
    blanks. The line format is then specified by only the order and
    the units of the columns. As the catalogue entry for each
    transition can be quite long, it can be broken across lines in the
    ASCII file. That is why each new transition is marked with a `@'
    character.

    The first column will contain the species and isotope, following
    the naming scheme described below.  Scientific notation is
    allowed, e.g. 501.12345e9.  

    Note that starting with ARTSCAT-2, the intensity is per molecule,
    i.e., it does not contain the isotope ratio. This is similar to
    JPL, but different to HITRAN.

    The line format is:

    \verbatim
    Col  Variable                Label    Unit     Comment
    ------------------------------------------------------------------      
     0   `@'                         ENTRY        -    marks start of entry
     1   name                         NAME        -    e.g. O3-666
     2   center frequency                F       Hz    e.g. 501.12345e9 
     3   pressure shift of F           PSF    Hz/Pa    
     4   line intensity                 I0    Hz*m^2   per isotope, not per species
     5   reference temp. for I0       T_I0        K
     6   lower state energy           ELOW        J    
     7   air broadened width          AGAM    Hz/Pa    values around 20000 Hz/Pa
     8   self broadened width         SGAM    Hz/Pa
     9   AGAM temp. exponent          NAIR        -    values around .5
    10   SGAM temp. exponent         NSELF        - 
    11   ref. temp. for AGAM, SGAM   T_GAM        K
    12   number of aux. parameters   N_AUX        -
    13   auxiliary parameter          AUX1        -
    14   ... 
    15   error for F                    DF       Hz
    16   error for I0                  DI0        %
    17   error for AGAM              DAGAM        %
    18   error for SGAM              DSGAM        %
    19   error for NAIR              DNAIR        %
    20   error for NSELF            DNSELF        %
    21   error for PSF                DPSF        %
    \endverbatim

    The parameters 0-12 must be present, the others can be missing,
    since they are not needed for the calculation.

    For the error fields (15-21), a -1 means that no value exist.

    A valid ARTS line file would be:
    \verbatim
    <?xml version="1.0"?>
    <arts format="ascii" version="1">
    <ArrayOfLineRecord version="ARTSCAT-3" nelem="2">
    @ O3-676 80015326542.0992 0 3.70209114155527e-19 296 7.73661776567701e-21 21480.3182341969
    28906.7092490501 0.76 0.76 296 0 300000 0.1 0.1 0.1 0.2 -1 -1.24976056865038e-11
    @ O3-676 80015476438.3282 0 3.83245786810611e-19 296 7.73661577919822e-21 21480.3182341969
    28906.7092490501 0.76 0.76 296 0 300000 0.1 0.1 0.1 0.2 -1 -1.24975822742335e-11
    </ArrayOfLineRecord>
    </arts>
    \endverbatim

    Some species need special parameters that are not needed by other
    species (for example overlap coefficients for O2). In the case of
    oxygen two parameters are sufficient to describe the overlap, but
    other species, e.g., methane, may need more coefficients. The
    default for \c N_AUX is zero. In that case, no further
    \c AUX fields are present. [FIXME: Check Oxygen.]

    The names of the private members and public access functions of
    this data structure follow the above table. The only difference is
    that underscores are omited and only the first letter of each name
    is capitalized. This is for consistency with the notation
    elsewhere in the program.

    \author Stefan Buehler 
*/
class LineRecord {
public:

  /** Make the helper function linesElowToJoule a friend, so that it
      can change the lower state energy. FIXME: Remove this when we
      don't need the helper function anymore. */
  friend void linesElowToJoule(Array<LineRecord>& abs_lines);

  /** Default constructor. Initialize to default values. The indices
      are initialized to large numbers, so that we at least get range
      errors when we try to used un-initialized data. */
  LineRecord()
    : mversion (3),
      mspecies (1000000),
      misotope (1000000),
      mf       (0.     ),
      mpsf     (0.     ),
      mi0      (0.     ),
      mti0     (0.     ),
      melow    (0.     ),
      magam    (0.     ),
      msgam    (0.     ),
      mnair    (0.     ),
      mnself   (0.     ),
      mtgam    (0.     ),
      maux     (       ),
      mdf      (-1.    ),
      mdi0     (-1.    ),
      mdagam   (-1.    ),
      mdsgam   (-1.    ),
      mdnair   (-1.    ),
      mdnself  (-1.    ),
      mdpsf    (-1.    )
 { /* Nothing to do here. */ }

  /** Constructor that sets all data elements explicitly. If
      assertions are not disabled (i.e., if NDEBUG is not \#defined),
      assert statements check that the species and isotope data
      exists. */
  LineRecord( Index                 species,
              Index                 isotope,
              Numeric               f,
              Numeric               psf,
              Numeric               i0,
              Numeric               ti0,
              Numeric               elow,
              Numeric               agam,
              Numeric               sgam,
              Numeric               nair,
              Numeric               nself,
              Numeric               tgam,
              const ArrayOfNumeric& aux,
              Numeric               /* df */,
              Numeric               /* di0 */,
              Numeric               /* dagam */,
              Numeric               /* dsgam */,
              Numeric               /* dnair */,
              Numeric               /* dnself */,
              Numeric               /* dpsf */)
    : mversion (3),
      mspecies (species    ),
      misotope (isotope    ),
      mf       (f          ),
      mpsf     (psf        ),
      mi0      (i0         ),
      mti0     (ti0        ),
      melow    (elow       ),
      magam    (agam       ),
      msgam    (sgam       ),
      mnair    (nair       ),
      mnself   (nself      ),
      mtgam    (tgam       ), 
      maux     (aux        ),
      mdf      (-1.    ),
      mdi0     (-1.    ),
      mdagam   (-1.    ),
      mdsgam   (-1.    ),
      mdnair   (-1.    ),
      mdnself  (-1.    ),
      mdpsf    (-1.    )
  {
    // Thanks to Matpack, initialization of misotope with isotope
    // should now work correctly.  

    // Check if this species is legal, i.e., if species and isotope
    // data exists.
    ////    extern Array<SpeciesRecord> species_data;
    //assert( mspecies < species_data.nelem() );
    //assert( misotope < species_data[mspecies].Isotope().nelem() );
  }

  /** Return the version String. */
  String VersionString() const
  {
    ostringstream os;
    os << "ARTSCAT-" << mversion;
    return os.str();
  }
  
  /** Return the version number. */
  Index Version() const { return mversion; }
  
  /** The index of the molecular species that this line belongs
      to. The species data can be accessed by species_data[Species()]. */
  Index Species() const { return mspecies; }

  /** The index of the isotopic species that this line belongs
      to. The isotopic species data can be accessed by
      species_data[Species()].Isotope()[Isotope()].  */
  Index Isotope() const { return misotope; }

  String Name() const;

  const SpeciesRecord& SpeciesData() const;

  const IsotopeRecord& IsotopeData() const;

  /** The line center frequency in <b> Hz</b>. */
  Numeric F() const     { return mf; }

  /** Set the line center frequency in <b> Hz</b>. */
  void setF( Numeric new_mf ) { mf = new_mf; }

  /** The pressure shift parameter in <b> Hz/Pa</b>. */
  Numeric Psf() const   { return mpsf; }

  /** Set the pressure shift parameter in <b> Hz/Pa</b>. */
  void setPsf( Numeric new_mpsf ) { mpsf = new_mpsf; }

  /** The line intensity in <b> m^2*Hz</b> at the reference temperature \c Ti0. 

    The line intensity \f$I_0\f$ is defined by:
    
    \f[
    \alpha(\nu) = n \, x \, I_0(T) \, F(\nu)
    \f]

    where \f$\alpha\f$ is the absorption coefficient (in <b>
    m^-1</b>), \f$\nu\f$ is frequency, \f$n\f$ is the
    total number density, \f$x\f$ is the volume mixing ratio, and
    \f$F(\nu)\f$ is the lineshape function. */
  Numeric I0() const    { return mi0; }

 /** Set Intensity */
  void setI0( Numeric new_mi0 ) { mi0 = new_mi0; }

  /** Reference temperature for I0 in <b> K</b>: */
  Numeric Ti0() const   { return mti0; }

  /** Lower state energy in <b> cm^-1</b>: */
  Numeric Elow() const  { return melow; }

  /** Air broadened width in <b> Hz/Pa</b>: */
  Numeric Agam() const  { return magam; }

   /** Set Air broadened width in <b> Hz/Pa</b>: */
  void setAgam( Numeric new_agam ) { magam = new_agam; }

  /** Self broadened width in <b> Hz/Pa</b>: */
  Numeric Sgam() const  { return msgam; }

  /** Set Self  broadened width in <b> Hz/Pa</b>: */
  void setSgam( Numeric new_sgam ) { msgam = new_sgam; }

  /** AGAM temperature exponent (dimensionless): */
  Numeric Nair() const  { return mnair; }

  /** Set AGAM temperature exponent (dimensionless): */
  void setNair( Numeric new_mnair ) { mnair = new_mnair; }

  /** SGAM temperature exponent (dimensionless): */
  Numeric Nself() const { return mnself; }

 /** Set SGAM temperature exponent (dimensionless): */
  void setNself( Numeric new_mnself ) { mnself = new_mnself; }

  /** Reference temperature for AGAM and SGAM in <b> K</b>: */
  Numeric Tgam() const  { return mtgam; }

  /** Number of auxiliary parameters. This function is actually
      redundant, since the number of auxiliary parameters can also be
      obtained directly with Aux.nelem(). I just added the function in
      order to have consistency of the interface with the catalgue
      format. */
  Index Naux() const   { return maux.nelem(); }

  /** Auxiliary parameters. */
  const ArrayOfNumeric& Aux() const { return maux; }
  //

  /** Accuracy for line position in <b> Hz </b>: */
  Numeric dF() const  { return mdf; }

  /** Accuracy for line intensity in <b> relative value </b>: */
  Numeric dI0() const  { return mdi0; }

 /** Accuracy for air broadened width in <b> relative value </b>: */
  Numeric dAgam() const  { return mdagam; }

  /** Accuracy for self broadened width in <b> relative value </b>: */
  Numeric dSgam() const  { return mdsgam; }

  /** Accuracy for AGAM temperature exponent in <b> relative value </b>: */
  Numeric dNair() const  { return mdnair; }

  /** Accuracy for SGAM temperature exponent in <b> relative value</b>: */
  Numeric dNself() const { return mdnself; }

  /** Accuracy for pressure shift in <b> relative value </b>: */
  Numeric dPsf() const { return mdpsf; }

  /** Einstein A-coefficient in <b> 1/s </b>: */
  Numeric A() const { return ma; }
  
  /** Upper state stat. weight: */
  Numeric G_upper() const { return mgupper; }
  
  /** Lower state stat. weight: */
  Numeric G_lower() const { return mglower; }
  
  /** Broadening parameter self in <b> Hz/Pa </b>: */
  Numeric Gamma_self() const { return mgamma_self; }
  
  /** Broadening parameter N2 in <b> Hz/Pa </b>: */
  Numeric Gamma_N2() const { return mgamma_n2; }
  
  /** Broadening parameter O2 in <b> Hz/Pa </b>: */
  Numeric Gamma_O2() const { return mgamma_o2; }
  
  /** Broadening parameter H2O in <b> Hz/Pa </b>: */
  Numeric Gamma_H2O() const { return mgamma_h2o; }
  
  /** Broadening parameter CO2 in <b> Hz/Pa </b>: */
  Numeric Gamma_CO2() const { return mgamma_co2; }
  
  /** Broadening parameter H2 in <b> Hz/Pa </b>: */
  Numeric Gamma_H2() const { return mgamma_h2; }
  
  /** Broadening parameter He in <b> Hz/Pa </b>: */
  Numeric Gamma_He() const { return mgamma_he; }
  
  /** GAM temp. exponent N self: */
  Numeric Gam_N_self() const { return mn_self; }
  
  /** GAM temp. exponent N N2: */
  Numeric Gam_N_N2() const { return mn_n2; }
  
  /** GAM temp. exponent N O2: */
  Numeric Gam_N_O2() const { return mn_o2; }
  
  /** GAM temp. exponent N H2O: */
  Numeric Gam_N_H2O() const { return mn_h2o; }
  
  /** GAM temp. exponent N CO2: */
  Numeric Gam_N_CO2() const { return mn_co2; }
  
  /** GAM temp. exponent N H2: */
  Numeric Gam_N_H2() const { return mn_h2; }
  
  /** GAM temp. exponent N He: */
  Numeric Gam_N_He() const { return mn_he; }
  
  /** F Pressure shift N2 in <b> Hz/Pa </b>: */
  Numeric Delta_N2() const { return mdelta_n2; }
  
  /** F Pressure shift O2 in <b> Hz/Pa </b>: */
  Numeric Delta_O2() const { return mdelta_o2; }
  
  /** F Pressure shift H2O in <b> Hz/Pa </b>: */
  Numeric Delta_H2O() const { return mdelta_h2o; }
  
  /** F Pressure shift CO2 in <b> Hz/Pa </b>: */
  Numeric Delta_CO2() const { return mdelta_co2; }
  
  /** F Pressure shift H2 in <b> Hz/Pa </b>: */
  Numeric Delta_H2() const { return mdelta_h2; }
  
  /** F Pressure shift He in <b> Hz/Pa </b>: */
  Numeric Delta_He() const { return mdelta_he; }
  
  /** Read one line from a stream associated with a HITRAN 1986-2001 file. The
    HITRAN format is as follows (directly from the HITRAN documentation):

    \verbatim
    Each line consists of 100
    bytes of ASCII text data, followed by a line feed (ASCII 10) and
    carriage return (ASCII 13) character, for a total of 102 bytes per line.
    Each line can be read using the following READ and FORMAT statement pair
    (for a FORTRAN sequential access read):

          READ(3,800) MO,ISO,V,S,R,AGAM,SGAM,E,N,d,V1,V2,Q1,Q2,IERF,IERS,
         *  IERH,IREFF,IREFS,IREFH
    800   FORMAT(I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2I3,2A9,3I1,3I2)

    Each item is defined below, with its format shown in parenthesis.

      MO  (I2)  = molecule number
      ISO (I1)  = isotope number (1 = most abundant, 2 = second, etc)
      V (F12.6) = frequency of transition in wavenumbers (cm-1)
      S (E10.3) = intensity in cm-1/(molec * cm-2) at 296 Kelvin
      R (E10.3) = transition probability squared in Debyes**2
      AGAM (F5.4) = air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      SGAM (F5.4) = self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
      E (F10.4) = lower state energy in wavenumbers (cm-1)
      N (F4.2) = coefficient of temperature dependence of air-broadened halfwidth
      d (F8.6) = shift of transition due to pressure (cm-1)
      V1 (I3) = upper state global quanta index
      V2 (I3) = lower state global quanta index
      Q1 (A9) = upper state local quanta
      Q2 (A9) = lower state local quanta
      IERF (I1) = accuracy index for frequency reference
      IERS (I1) = accuracy index for intensity reference
      IERH (I1) = accuracy index for halfwidth reference
      IREFF (I2) = lookup index for frequency
      IREFS (I2) = lookup index for intensity
      IREFH (I2) = lookup index for halfwidth

    The molecule numbers are encoded as shown in the table below:

      0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
      6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
     12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
     18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
     24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
     30=  SF6   31=  H2S   32=HCOOH
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped. It is checked if the data record has the right
    number of characters. If not, a runtime_error is thrown.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \author Stefan Buehler */
  bool ReadFromHitranStream(istream& is, const Verbosity& verbosity);



  /** Read one line from a stream associated with a HITRAN 2004 file. The
    HITRAN format is as follows:

    \verbatim
    Each line consists of 160 ASCII characters, followed by a line feed (ASCII 10)
    and carriage return (ASCII 13) character, for a total of 162 bytes per line.

    Each item is defined below, with its Fortran format shown in parenthesis.

    (I2)     molecule number
    (I1)     isotopologue number (1 = most abundant, 2 = second, etc)
    (F12.6)  vacuum wavenumbers (cm-1)
    (E10.3)  intensity in cm-1/(molec * cm-2) at 296 Kelvin
    (E10.3)  Einstein-A coefficient (s-1)
    (F5.4)   air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    (F5.4)   self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    (F10.4)  lower state energy (cm-1)
    (F4.2)   coefficient of temperature dependence of air-broadened halfwidth
    (F8.6)   air-broadened pressure shift of line transition at 296 K (cm-1)
    (A15)    upper state global quanta
    (A15)    lower state global quanta
    (A15)    upper state local quanta
    (A15)    lower state local quanta
    (I1)     uncertainty index for wavenumber
    (I1)     uncertainty index for intensity
    (I1)     uncertainty index for air-broadened half-width
    (I1)     uncertainty index for self-broadened half-width
    (I1)     uncertainty index for temperature dependence
    (I1)     uncertainty index for pressure shift
    (I2)     index for table of references correspond. to wavenumber
    (I2)     index for table of references correspond. to intensity
    (I2)     index for table of references correspond. to air-broadened half-width
    (I2)     index for table of references correspond. to self-broadened half-width
    (I2)     index for table of references correspond. to temperature dependence
    (I2)     index for table of references correspond. to pressure shift
    (A1)     flag (*) for lines supplied with line-coupling algorithm
    (F7.1)   upper state statistical weight
    (F7.1)   lower state statistical weight

     The molecule numbers are encoded as shown in the table below:

      0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=    CO
      6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=   NH3
     12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=    HI
     18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=   HCN
     24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29=  COF2
     30=  SF6   31=  H2S   32=HCOOH   33=  HO2   34=    O   35=ClONO2
     36=  NO+   37= HOBr   38= C2H4
    \endverbatim

    CH3OH is not included in ARTS because its total internal partition
    sum is not known yet.

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped. It is checked if the data record has the right number
    of characters (comment lines are ignored). If not, a runtime_error is thrown.
    If the molecule is unknown to ARTS, a warning is prompted but the
    program continues (ignoring this record). For CH3OH this
    warning will be issued even when using the regular Hitran 2004 data base
    (see above).

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \author Stefan Buehler, Hermann Berg */
  bool ReadFromHitran2004Stream(istream& is, const Verbosity& verbosity);




  /** Read one line from a stream associated with a MYTRAN2 file. The MYTRAN2
    format is as follows (directly taken from the abs_my.c documentation):

    \verbatim
    The MYTRAN format is as follows (FORTRAN notation):
    FORMAT(I2,I1,F13.4,1PE10.3,0P2F5.2,F10.4,2F4.2,F8.6,F6.4,2I3,2A9,4I1,3I2)
   
    Each item is defined below, with its FORMAT String shown in
    parenthesis.
   
       MO  (I2)      = molecule number
       ISO (I1)      = isotope number (1 = most abundant, 2 = second, etc)
    *  F (F13.4)     = frequency of transition in MHz
    *  errf (F8.4)   = error in f in MHz
       S (E10.3)     = intensity in cm-1/(molec * cm-2) at 296 K
    *  AGAM (F5.4)   = air-broadened halfwidth (HWHM) in MHz/Torr at Tref
    *  SGAM (F5.4)   = self-broadened halfwidth (HWHM) in MHz/Torr at Tref
       E (F10.4)     = lower state energy in wavenumbers (cm-1)
       N (F4.2)      = coefficient of temperature dependence of 
                       air-broadened halfwidth
    *  N_self (F4.2) = coefficient of temperature dependence of 
                       self-broadened halfwidth
    *  Tref (F7.2)   = reference temperature for AGAM and SGAM 
    *  d (F9.7)      = shift of transition due to pressure (MHz/Torr)
       V1 (I3)       = upper state global quanta index
       V2 (I3)       = lower state global quanta index
       Q1 (A9)       = upper state local quanta
       Q2 (A9)       = lower state local quanta
       IERS (I1)     = accuracy index for S
       IERH (I1)     = accuracy index for AGAM
    *  IERN (I1)     = accuracy index for N

   
    The asterisks mark entries that are different from HITRAN.

    Note that AGAM and SGAM are for the temperature Tref, while S is
    still for 296 K!
   
    The molecule numbers are encoded as shown in the table below:
   
     0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
     6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
    30=  SF6   31=  H2S   32=HCOOH   33= HO2    34=    O   35= CLONO2
    36=  NO+   37= Null   38= Null   39= Null   40=H2O_L   41= Null
    42= Null   43= OCLO   44= Null   45= Null   46=BRO     47= Null
    48= H2SO4  49=CL2O2

    All molecule numbers are from HITRAN, except for species with id's
    greater or equals 40, which are not included in HITRAN.
    (E.g.: For BrO, iso=1 is Br-79-O,iso=2 is  Br-81-O.)
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped.

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 31.10.00
    \author Axel von Engeln 
  */
  bool ReadFromMytran2Stream(istream& is, const Verbosity& verbosity);


  /** Read one line from a stream associated with a JPL file. The JPL
    format is as follows (directly taken from the jpl documentation):

    \verbatim 
    The catalog line files are composed of 80-character lines, with one
    line entry per spectral line.  The format of each line is:

    \label{lfmt}
    \begin{tabular}{@{}lccccccccr@{}}
    FREQ, & ERR, & LGINT, & DR, & ELO, & GUP, & TAG, & QNFMT, & QN${'}$, & QN${''}$\\ 
    (F13.4, & F8.4, & F8.4, & I2, & F10.4, & I3, & I7, & I4, & 6I2, & 6I2)\\
    \end{tabular}

    \begin{tabular}{lp{4.5in}} 
    FREQ: & Frequency of the line in MHz.\\ 
    ERR: & Estimated or experimental error of FREQ in MHz.\\ 
    LGINT: &Base 10 logarithm of the integrated intensity 
    in units of \linebreak nm$^2$$\cdot$MHz at 300 K. (See Section 3 for 
    conversions to other units.)\\ 
    DR: & Degrees of freedom in the rotational partition 
    function (0 for atoms, 2 for linear molecules, and 3 for nonlinear 
    molecules).\\ 
    ELO: &Lower state energy in cm$^{-1}$ relative to the lowest energy 
    spin--rotation level in ground vibronic state.\\ 
    GUP: & Upper state degeneracy.\\ 
    TAG: & Species tag or molecular identifier. 
    A negative value flags that the line frequency has 
    been measured in the laboratory.  The absolute value of TAG is then the 
    species tag and ERR is the reported experimental error.  The three most 
    significant digits of the species tag are coded as the mass number of the 
    species, as explained above.\\ 
    QNFMT: &Identifies the format of the quantum numbers 
    given in the field QN. These quantum number formats are given in Section 5 
    and are different from those in the first two editions of the catalog.\\ 
    QN${'}$: & Quantum numbers for the upper state coded 
    according to QNFMT.\\ 
    QN${''}$: & Quantum numbers for the lower state.\\
    \end{tabular} 
    \endverbatim

    The function attempts to read a line of data from the
    catalogue. It returns false if it succeeds. Otherwise, if eof is
    reached, it returns true. If an error occurs, a runtime_error is
    thrown. When the function looks for a data line, comment lines are
    automatically skipped (unused in jpl).

    \param is Stream from which to read
    \exception runtime_error Some error occured during the read
    \return false=ok (data returned), true=eof (no data returned)

    \date 01.11.00
    \author Axel von Engeln */
  bool ReadFromJplStream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with an ARTSCAT-3 file.

      Format: see Documentation of class LineRecord

      The function attempts to read a line of data from the
      catalogue. It returns false if it succeeds. Otherwise, if eof is
      reached, it returns true. If an error occurs, a runtime_error is
      thrown. When the function looks for a data line, comment lines are
      automatically skipped.

      \param is Stream from which to read
      \exception runtime_error Some error occured during the read
      \return false=ok (data returned), true=eof (no data returned)

      \date   2000-12-15
      \author Oliver Lemke

      
      \date   2001-06-20
      \author Stefan Buehler
*/
  bool ReadFromArtscat3Stream(istream& is, const Verbosity& verbosity);

  /** Read one line from a stream associated with an ARTSCAT-4 file.
   
   Format: see Documentation of class LineRecord
   
   The function attempts to read a line of data from the
   catalogue. It returns false if it succeeds. Otherwise, if eof is
   reached, it returns true. If an error occurs, a runtime_error is
   thrown. When the function looks for a data line, comment lines are
   automatically skipped.
   
   \param is Stream from which to read
   \exception runtime_error Some error occured during the read
   \return false=ok (data returned), true=eof (no data returned)
   
   \date   2012-02-10
   \author Oliver Lemke
   */
  bool ReadFromArtscat4Stream(istream& is, const Verbosity& verbosity);
  

private:
  // Version number:
  Index mversion;
  // Molecular species index: 
  Index mspecies;
  // Isotopic species index:
  Index misotope;
  // The line center frequency in Hz:
  Numeric mf;
  // The pressure shift parameter in Hz/Pa:
  Numeric mpsf;
  // The line intensity in m^2/Hz:
  Numeric mi0;
  // Reference temperature for I0 in K:
  Numeric mti0;
  // Lower state energy in cm^-1:
  Numeric melow;
  // Air broadened width in Hz/Pa:
  Numeric magam;
  // Self broadened width in Hz/Pa:
  Numeric msgam;
  // AGAM temperature exponent (dimensionless):
  Numeric mnair;
  // SGAM temperature exponent (dimensionless):
  Numeric mnself;
  // Reference temperature for AGAM and SGAM in K:
  Numeric mtgam;
  // Array to hold auxiliary parameters:
  ArrayOfNumeric maux;
  //
  // Fields for the spectroscopic parameters accuracies
  //
  // Accuracy for line center frequency in Hz:
  Numeric mdf;
  // Accuracy for line intensity in %:
  Numeric mdi0;
  // Accuracy for air broadened width in %:
  Numeric mdagam;
  // Accuracy for self broadened width in %:
  Numeric mdsgam;
  // Accuracy for AGAM temperature exponent in %:
  Numeric mdnair;
  //  Accuracy for SGAM temperature exponent in %:
  Numeric mdnself;
 //  Accuracy for pressure shift in %:
  Numeric mdpsf;
  
  //// New fields in ARTSCAT-4
  
  // Einstein A-coefficient in 1/s:
  Numeric ma;
  // Upper state stat. weight:
  Numeric mgupper;
  // Lower state stat. weight:
  Numeric mglower;
  
  // Broadening parameter self in Hz/Pa:
  Numeric mgamma_self;
  // Broadening parameter N2 in Hz/Pa:
  Numeric mgamma_n2;
  // Broadening parameter O2 in Hz/Pa:
  Numeric mgamma_o2;
  // Broadening parameter H2O in Hz/Pa:
  Numeric mgamma_h2o;
  // Broadening parameter CO2 in Hz/Pa:
  Numeric mgamma_co2;
  // Broadening parameter H2 in Hz/Pa:
  Numeric mgamma_h2;
  // Broadening parameter He in Hz/Pa:
  Numeric mgamma_he;

  // GAM temp. exponent self:
  Numeric mn_self;
  // GAM temp. exponent N2:
  Numeric mn_n2;
  // GAM temp. exponent O2:
  Numeric mn_o2;
  // GAM temp. exponent H2O:
  Numeric mn_h2o;
  // GAM temp. exponent CO2:
  Numeric mn_co2;
  // GAM temp. exponent H2:
  Numeric mn_h2;
  // GAM temp. exponent He:
  Numeric mn_he;
  
  // F Pressure shift N2 in Hz/Pa:
  Numeric mdelta_n2;
  // F Pressure shift O2 in Hz/Pa:
  Numeric mdelta_o2;
  // F Pressure shift H2O in Hz/Pa:
  Numeric mdelta_h2o;
  // F Pressure shift CO2 in Hz/Pa:
  Numeric mdelta_co2;
  // F Pressure shift H2 in Hz/Pa:
  Numeric mdelta_h2;
  // F Pressure shift He in Hz/Pa:
  Numeric mdelta_he;
  
};

// is needed to map jpl tags/arts identifier to the species/isotope data within arts
class SpecIsoMap{
public:
  SpecIsoMap():mspeciesindex(0), misotopeindex(0){}
  SpecIsoMap(const Index& speciesindex,
                const Index& isotopeindex)
    : mspeciesindex(speciesindex),
      misotopeindex(isotopeindex) 
  {}

  // Return the index to the species 
  const Index& Speciesindex() const { return mspeciesindex; }
  // Return the index to the isotope
  const Index& Isotopeindex() const { return misotopeindex; }

private:
  Index mspeciesindex;
  Index misotopeindex;
};



/** Holds a list of spectral line data.
    \author Stefan Buehler */
typedef Array<LineRecord> ArrayOfLineRecord;

/** Holds a lists of spectral line data for each tag group.
    Dimensions: (tag_groups.nelem()) (number of lines for this tag)
    \author Stefan Buehler */
typedef Array< Array<LineRecord> > ArrayOfArrayOfLineRecord;



/** Output operator for LineRecord. The result should look like a
    catalogue line.

    \author Stefan Buehler */
ostream& operator<< (ostream& os, const LineRecord& lr);



/** Define the species data map.

    \author Stefan Buehler  */
void define_species_map();


// Doc header in absorption.cc
void write_lines_to_stream(ostream& os,
                           const ArrayOfLineRecord& abs_lines);


void xsec_species( MatrixView               xsec,
                   ConstVectorView          f_mono,
                   ConstVectorView          abs_p,
                   ConstVectorView          abs_t,           
                   ConstVectorView          h2o_abs,           
                   ConstVectorView          vmr,
                   const ArrayOfLineRecord& abs_lines,
                   const Index              ind_ls,
                   const Index              ind_lsn,
                   const Numeric            cutoff,
                   const Verbosity&         verbosity );


// A helper function for energy conversion:
Numeric wavenumber_to_joule(Numeric e);



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

#endif // absorption_h
