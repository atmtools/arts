// This file contains declarations required for the calculation of
// absorption ciefficients.  

#ifndef absorption_h
#define absorption_h

#include "vecmat.h"


/** Contains the lookup data for one isotope.
    @author Stefan Buehler */
class IsotopeRecord{
public:
  /** Default constructor. Needed by make_array. */
  IsotopeRecord() { /* Nothing to do here */ };
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
	strstream is;
	is << mname;
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
	//	cout << "sum = " << sum << "  irm%10 = " << irm%10 << endl;
	assert( sum == irm%10 );

	// Same test for JPL tag numbers:
	for ( size_t i=0; i<mjpltags.size(); ++i )
	  {
	    //	    cout << "sum = " << sum << "  (mjpltags[i] / 1000) % 10 = " << (mjpltags[i] / 1000) % 10 << endl;
	    assert( sum == (mjpltags[i] / 1000) % 10 );
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
	    assert( misotope[i].Abundance() > misotope[i+1].Abundance() );
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
    catalogue format, taken from the Bredbeck book:

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

    \begin{verbatim}
    Col  Variable                Label    Unit     Comment
    ------------------------------------------------------------------      
    0   `@'                      ENTRY      -     marks start of entry
    1   name                     NAME       -     e.g. O3-666
    2   center frequency            F      Hz     e.g. 501.12345e9 
    3   pressure shift of F       PSF  MHz/hPa    
    4   line intensity             I0  m^2/Hz 
    5   reference temp. for I0   T_I0       K
    6   lower state energy       ELOW    cm-1
    7   air broadened width      AGAM  MHz/hPa     values around 2
    8   self broadened width     SGAM  MHz/hPa
    9   AGAM temp. exponent      NAIR       -
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
    \end{verbatim}

    One line could be:
    {\small
    \begin{verbatim}
    @ O3-666 110.83604e9 0 12.43453e-22 300 17.5973 1.52 2.03 0.73 0.73 296 ....   
    \end{verbatim}}

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

  /** The pressure shift parameter in #MHz/hPa#. */
  Numeric Psf() const   { return mpsf; }

  /** The line intensity in #m^2/Hz#. */
  Numeric I0() const    { return mi0; }

  /** Reference temperature for I0 in #K#: */
  Numeric Ti0() const   { return mti0; }

  /** Lower state energy in #cm^-1#: */
  Numeric Elow() const  { return melow; }

  /** Air broadened width in #MHz/hPa#: */
  Numeric Agam() const  { return magam; }

  /** Self broadened width in #MHz/hPa#: */
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

    \begin{verbatim}
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
    \end{verbatim}

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
  // The pressure shift parameter in MHz/hPa:
  Numeric mpsf;
  // The line intensity in m^2/Hz:
  Numeric mi0;
  // Reference temperature for I0 in K:
  Numeric mti0;
  // Lower state energy in cm^-1:
  Numeric melow;
  // Air broadened width in MHz/hPa:
  Numeric magam;
  // Self broadened width in MHz/hPa:
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





#endif // absorption_h
