/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

/** @file
    Declarations required for the calculation of absorption ciefficients.

    @author Stefan Buehler
*/

// This file contains 

#ifndef absorption_h
#define absorption_h

#include <iostream>
#include "vecmat.h"

/** Contains the lookup data for one isotope.
    @author Stefan Buehler */
class IsotopeRecord{
public:
  /** Default constructor. Needed by make_array. */
  IsotopeRecord() { /* Nothing to do here */ }
  /** Constructor that sets the values. */
  IsotopeRecord(const string&  	  name,
		const Numeric& 	  abundance,
		const Numeric& 	  mass,
		const int&     	  mytrantag,
		const int&     	  hitrantag,
		const ARRAY<int>& jpltags)
    : mname(name),
      mabundance(abundance),
      mmass(mass),
      mmytrantag(mytrantag),
      mhitrantag(hitrantag),
      mjpltags(jpltags)
  {
    // Some consistency checks whether the given data makes sense.
#ifndef NDEBUG
      {
	/* 1. All the tags must be positive or -1 */
	assert( (0<mmytrantag) || (-1==mmytrantag) );
	assert( (0<mhitrantag) || (-1==mhitrantag) );
	for ( size_t i=0; i<mjpltags.size(); ++i )
	  assert( (0<mjpltags[i]) || (-1==mjpltags[i]) );


	/* 2. Check, whether the isotope name (e.g., 686) is consistent with the mass.
	      (And with JPL tag numbers.) */

	/* Convert string to number: */
	int name_as_number;
	istringstream is(mname);
	is >> name_as_number;
	//	cout << "name_as_number = " << name_as_number << endl;

	// Collect the sum of the individual digits:
	int sum = 0;
	while (0 < name_as_number)
	  {
	    sum += name_as_number % 10;
	    name_as_number /= 10;
	  }

	//	cout << "sum = " << sum << endl;
	sum %= 10;

	// Make sure that modulo 10, sum and rounded mass are the
	// same.  This is so complicated, because we have to convert
	// the mass to an integer first.
	Numeric rm  = floor(mmass+.5);
	int    irm = static_cast<int>(rm); 
	//	assert( sum == irm%10 );
	// SAB&AvE 08.08.2000 Disabled this check, since not
	// always all atoms are explicitly listed in the isotope name.

	// Same test for JPL tag numbers:
	for ( size_t i=0; i<mjpltags.size(); ++i )
	  {
	    //	    assert( sum == (mjpltags[i] / 1000) % 10 );
	    // SAB&AvE 08.08.2000 Disabled this check, since not
	    // always all atoms are explicitly listed in the isotope name.
	  }
      }
#endif // ifndef NDEBUG
  }

  const string&       Name()         const { return mname;  }
  const Numeric&      Abundance()    const { return mabundance; }
  const Numeric&      Mass()         const { return mmass;    }
  const int&          MytranTag()    const { return mmytrantag;    }
  const int&          HitranTag()    const { return mhitrantag;    }
  const ARRAY<int>&   JplTags()      const { return mjpltags;      }
  
private:
  /** Isotope names. */
  string mname;
  /** Normal abundance ( = isotopic ratio). (Absolute number.) */
  Numeric mabundance;
  /** Mass of the isotope. (In unified atomic mass units u)
      If I understand this correctly this is the same as g/mol. */
  Numeric mmass;
  /** MYTRAN2 tag numers for all isotopes. -1 means not included. */
  int mmytrantag;
  /** HITRAN-96 tag numers for all isotopes. -1 means not included. */
  int mhitrantag;
  /** JPL tag numbers for all isotopes. Empty array means not included. There
      can be more than one JPL tag for an isotopic species, because in
      JPL different vibrational states have different tags. */
  ARRAY<int> mjpltags;
};


/** Contains the lookup data for one species.

    @author Stefan Buehler  */
class SpeciesRecord{
public:
  /** The constructor used in define_species_data. */
  SpeciesRecord(const char name[],
		const int degfr,
		const ARRAY<IsotopeRecord>& isotope)
    : mname(name),
      mdegfr(degfr),
      misotope(isotope)
  {
#ifndef NDEBUG
      {
	/* Check that the isotopes are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
	  {
	    assert( misotope[i].Abundance() >= misotope[i+1].Abundance() );
	  }

	/* Check that the Mytran tags are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
	  {
	    if ( (0<misotope[i].MytranTag()) && (0<misotope[i+1].MytranTag()) )
	      {
		assert( misotope[i].MytranTag() < misotope[i+1].MytranTag() );
	    
		// Also check that the tags have the same base number:
		assert( misotope[i].MytranTag()/10 == misotope[i].MytranTag()/10 );
	      }
	  }

	/* Check that the Hitran tags are correctly sorted. */
	for ( size_t i=0; i<misotope.size()-1; ++i )
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

  const string&               Name()     const { return mname;     }   
  int                         Degfr()    const { return mdegfr;    }
  const ARRAY<IsotopeRecord>& Isotope()  const { return misotope;  }
  
private:
  /** Species name. */
  string mname;
  /** Degrees of freedom. */
  int mdegfr;
  /** Isotope data. */
  ARRAY<IsotopeRecord> misotope;
};


/** Spectral line catalog data. Here is a description of the ARTS
    catalogue format, largely taken from the Bredbeck book, except
    that some units are slightly changed:

    The line catalogue should not have any fixed column widths because the
    precision of the parameters should not be limited by the format.  The
    catalogue can then be stored as binary or ASCII. In the ASCII version
    the columns are separated by one or more blanks. The line format is
    then specified by only the order and the units of the columns. As
    the catalogue entry for each transition can be quite long, it can be
    broken across lines in the ASCII file. Each new transition is marked
    with a `@' character.

    The first column will contain the species and isotope, following
    the naming scheme described below.  Scientific notation is allowed,
    e.g. 501.12345e9.  The transitions of different isotopes should be
    kept in separate files (this is for the catalogue, the forward model
    should also be able to use a single line file). The suggested line
    format is:

    \verbatim
    Col  Variable                Label    Unit     Comment
    ------------------------------------------------------------------      
    0   `@'                      ENTRY       -     marks start of entry
    1   name                     NAME        -     e.g. O3-666
    2   center frequency            F       Hz     e.g. 501.12345e9 
    3   pressure shift of F       PSF    Hz/Pa    
    4   line intensity             I0   m^2/Hz 
    5   reference temp. for I0   T_I0        K
    6   lower state energy       ELOW     cm-1
    7   air broadened width      AGAM    Hz/Pa     values around 2
    8   self broadened width     SGAM    Hz/Pa
    9   AGAM temp. exponent      NAIR        -
    10   SGAM temp. exponent     NSELF       - 
    11   ref. temp. for AGAM, SGAM T_GAM     K
    12   number of aux. parameters N_AUX     -
    13   auxiliary parameter       AUX1      -
    14   ...
    15   error for F                DF      Hz
    16   error for AGAM          DAGAM       %
    17   error for SGAM          DSGAM       %
    18   error for NAIR          DNAIR       %
    19   error for NSELF        DNSELF       %
    20   error for PSF            DPSF       %
    21   quantum number code                       string or number ?
    22   lower state quanta                        string inside quotes
    23   upper state quanta                        string inside quotes
    24   information source of F
    25   information source of I0
    26   information source of line width variables
    27   information source of overlap constants
    \endverbatim

    One line could be:
    {\small
    \verbatim
    @ O3-666 110.83604e9 0 12.43453e-22 300 17.5973 1.52 2.03 0.73 0.73 296 ....   
    \endverbatim}

    The format used in the line file used by the FM
    can be a truncated version of the full line format.

    Some species need special parameters that are not needed by other
    species (for example overlap coefficients for O2). In the case of
    oxygen two parameters are sufficient to describe the overlap, but
    other species, e.g., methane, may need more coefficients. The
    default for \texttt{N\_AUX} is zero. In that case, no further
    \texttt{AUX} fields are present.

    The names of the private members and public access functions of
    this data structure follow the above table. The only difference is
    that underscores are omited and only the first letter of each name
    is capitalized. This is for consistency with the notation
    elsewhere in the program.

    @author Stefan Buehler */
class LineRecord {
public:

  /** Default constructor. Initialize to default values. The indices
      are initialized to large numbers, so that we at least get range
      errors when we try to used un-initialized data. */
  LineRecord()
    : mspecies (1000000),
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
      maux     (       )
 { /* Nothing to do here. */ }

  /** Constructor that sets all data elements explicitly. If
      assertions are not disabled (i.e., if NDEBUG is not #defined),
      assert statements check that the species and isotope data
      exists. */
  LineRecord( size_t  	     	    species,
	      size_t  	     	    isotope,
	      Numeric 	     	    f,
	      Numeric 	     	    psf,
	      Numeric 	     	    i0,
	      Numeric 	     	    ti0,
	      Numeric 	     	    elow,
	      Numeric 	     	    agam,
	      Numeric 	     	    sgam,
	      Numeric 	     	    nair,
	      Numeric 	     	    nself,
	      Numeric 	     	    tgam,
	      const ARRAY<Numeric>& aux       )
    : mspecies (species),
      misotope (isotope),
      mf       (f      ),
      mpsf     (psf    ),
      mi0      (i0     ),
      mti0     (ti0    ),
      melow    (elow   ),
      magam    (agam   ),
      msgam    (sgam   ),
      mnair    (nair   ),
      mnself   (nself  ),
      mtgam    (tgam   ),
      maux     (aux    )
  {
    // Check if this species is legal, i.e., if species and isotope
    // data exists.
    extern const ARRAY<SpeciesRecord> species_data;
    assert( mspecies < species_data.size() );
    assert( misotope < species_data[mspecies].Isotope().size() );
  }

  /** The index of the molecular species that this line belongs
      to. The species data can be accessed by species_data[Species()]. */
  size_t Species() const { return mspecies; }

  /** The index of the isotopic species that this line belongs
      to. The isotopic species data can be accessed by
      species_data[Species()].Isotope()[Isotope()].  */
  size_t Isotope() const { return misotope; }

  /** The full name of the species and isotope. E.g., `O3-666'. The
      name is found by looking up the information in species_data,
      using the species and isotope index. */
  string Name() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    const SpeciesRecord& sr = species_data[mspecies];
    return sr.Name() + "-" + sr.Isotope()[misotope].Name();
  }

  /** The matching SpeciesRecord from species_data. To get at the
      species data of a LineRecord lr, you can use:
      \begin{enumerate}
      \item species_data[lr.Species()]
      \item lr.SpeciesData()
      \end{enumerate}
      The only advantages of the latter are that the notation is
      slightly nicer and that you don't have to declare the external
      variable species_data. */
  const SpeciesRecord& SpeciesData() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    return species_data[mspecies];
  }

  /** The matching IsotopeRecord from species_data. The IsotopeRecord
      is a subset of the SpeciesRecord. To get at the isotope data of
      a LineRecord lr, you can use:
      \begin{enumerate}
      \item species_data[lr.Species()].Isotope()[lr.Isotope()]
      \item lr.SpeciesData().Isotope()[lr.Isotope()]
      \item lr.IsotopeData()
      \end{enumerate}
      The last option is clearly the shortest, and has the advantage
      that you don't have to declare the external variable
      species_data. */
  const IsotopeRecord& IsotopeData() const {
    // The species lookup data:
    extern const ARRAY<SpeciesRecord> species_data;
    return species_data[mspecies].Isotope()[misotope];
  }

  /** The line center frequency in #Hz#. */
  Numeric F() const     { return mf; }

  /** The pressure shift parameter in #Hz/Pa#. */
  Numeric Psf() const   { return mpsf; }

  /** The line intensity in #m^2/Hz#. */
  Numeric I0() const    { return mi0; }

  /** Reference temperature for I0 in #K#: */
  Numeric Ti0() const   { return mti0; }

  /** Lower state energy in #cm^-1#: */
  Numeric Elow() const  { return melow; }

  /** Air broadened width in #Hz/Pa#: */
  Numeric Agam() const  { return magam; }

  /** Self broadened width in #Hz/Pa#: */
  Numeric Sgam() const  { return msgam; }

  /** AGAM temperature exponent (dimensionless): */
  Numeric Nair() const  { return mnair; }

  /** SGAM temperature exponent (dimensionless): */
  Numeric Nself() const { return mnself; }

  /** Reference temperature for AGAM and SGAM in #K#: */
  Numeric Tgam() const  { return mtgam; }

  /** Number of auxiliary parameters. This function is actually
      redundant, since the number of auxiliary parameters can also be
      obtained directly with Aux.size(). I just added the function in
      order to have consistency of the interface with the catalgue
      format. */
  size_t Naux() const   { return maux.size(); }

  /** Auxiliary parameters. */
  const ARRAY<Numeric>& Aux() const { return maux; }

  /** Read one line from a stream associated with a HITRAN file. The HITRAN
    format is as follows (directly from the HITRAN documentation):

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
    automatically skipped.

    @param is Stream from which to read
    @exception runtime_error Some error occured during the read
    @return false=ok (data returned), true=eof (no data returned)

    @author Stefan Buehler */
  bool ReadFromHitranStream(istream& is);

private:
  // Molecular species index: 
  size_t mspecies;
  // Isotopic species index:
  size_t misotope;
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
  ARRAY<Numeric> maux;
};



/** Holds a list of spectral line data.
    @author Stefan Buehler */
typedef ARRAY<LineRecord> ARRAYofLineRecord;

/** Holds a lists of spectral line data for each tag group.
    Dimensions: (tag_groups.dim()) (number of lines for this tag)
    @author Stefan Buehler */
typedef ARRAY< ARRAY<LineRecord> > ARRAYofARRAYofLineRecord;



/** Output operator for LineRecord. The result should look like a
    catalogue line.

    @author Stefan Buehler */
ostream& operator << (ostream& os, const LineRecord& lr);


/** Define species lookup data.

    Molecular masses: From looking at the table in forward_4_96,
    glob_def.c, the difference between the actual mass and the mass
    simply estimated from the Atom number is only 0.001. This seems
    not worth the trouble. Anyway, the field for the mass is
    there. Should anybody feel like adding the true numbers, just go
    ahead.

    @author Stefan Buehler  */
void define_species_data();

/** Define the species data map.

    @author Stefan Buehler  */
void define_species_map();


//------------------------------< Tag Group Stuff >------------------------------

/** A tag group can consist of the sum of several of these.

    @author Stefan Buehler */
class OneTag {
public:
  /** Default constructor. */
  OneTag() { /* Nothing to be done here. */ }

  /** Constructor from a tag definition string (Bredbeck
      notation). For examples see member function Name(). 

      @exception runtime_error The given string could not be mapped to
      a sensible tag description. */
  OneTag(string def); 

  /** Return the full name of this tag according to Bredbeck
      convention. Examples:
      \verbatim
      O3-*-*-*         : All O3 lines
      O3-666-*-*       : All O3-666 lines
      O3-*-500e9-501e9 : All O3 lines between 500 and 501 GHz.
      \endverbatim */
  string Name() const;
    
  /** Molecular species index. */
  size_t Species() const { return mspecies; }

  /** Isotopic species index.
      If this is equal to the number of isotopes (one more than
      allowed) it means all isotopes of this species. */ 
  size_t Isotope() const { return misotope; }

  /** The lower line center frequency in Hz.
      If this is <0 it means no lower limit. */
  Numeric Lf() const { return mlf; }

  /** The upper line center frequency in Hz:
      If this is <0 it means no upper limit. */
  Numeric Uf() const { return muf; }

private:
  // Molecular species index: 
  size_t mspecies;
  // Isotopic species index.
  // If this is equal to the number of isotopes (one more than
  // allowed) it means all isotopes of this species.
  size_t misotope;
  // The lower line center frequency in Hz.
  // If this is <0 it means no lower limit. 
  Numeric mlf;
  // The upper line center frequency in Hz:
  // If this is <0 it means no upper limit. 
  Numeric muf;
};


/** Output operator for OneTag. 

    @author Stefan Buehler */
ostream& operator << (ostream& os, const OneTag& ot);


/** Contains the available tag groups. Contrary to the Bredbeck
    definition, tag groups may only consist of tags belonging to the
    same species. The reason for this is that there is one VMR profile
    associated with each tag group.

    @author Stefan Buehler */
typedef  ARRAY< ARRAY<OneTag> > TagGroups;


/** A helper function that writes lines in a line list to a stream. 

    @param os       Output. The stream to write to.
    @param lines    The line list to write.

    @author Stefan Buehler 12.06.2000 */
void write_lines_to_stream(ostream& os,
			   const ARRAYofLineRecord& lines);


/** Calculate absorption coefficients for one tag group. All lines in
    the line list must belong to the same species. This must be
    ensured by lines_per_tgCreateFromLines, so it is only verified
    with assert. Also, the input vectors p_abs, t_abs, and vmr must
    all have the same dimension.

    This is a strongly simplified routine which seves mainly
    the purpose of demonstration.

    @param abs Output. Absorption coefficients.
    @param f_abs Frequency grid.
    @param p_abs Pressure grid.
    @param t_abs Temperatures associated with p_abs.
    @param vmrs  Volume mixing ratios of the species.
    @param lines The spectroscopic line list.

    @author Stefan Buehler 16.06.2000. */
void abs_species( MATRIX&                  abs,
		  const VECTOR&  	   f_abs,
		  const VECTOR&  	   p_abs,
		  const VECTOR&  	   t_abs,           
		  const VECTOR&            vmr,
		  const ARRAYofLineRecord& lines );

#endif // absorption_h
