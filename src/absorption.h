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
	// FIXME: Strange tricks, to circumvent non-ansi behaviour of
	// EGCS. Change this when it's no longer necessary to:
	// std::istrstream(mname) >> name_as_number  (should be sufficient)]
	istrstream is(&mname[0],mname.size());
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
