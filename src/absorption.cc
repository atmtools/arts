// Physical absorption routines. The absorption workspace methods are
// in file m_abs.cc 

#include "arts.h"
#include "make_array.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DEGFR(x) x
#define ISOTOPES     make_array<IsotopeRecord>
#define REC	     IsotopeRecord
#define TAGS	     make_array<int>



void define_species_data()
{
  extern ARRAY<SpeciesRecord> species_data;

  // Initialize to zero, just in case:
  species_data.clear();

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( ""	,		,	,	,	,TAGS() ),
	 ) ) );

  For good looks, keep the commas on the marks!

  */


  // H2O
  // Source: 
  // Isotope ratios: HITRAN92 tables.3
  //		     ISOTOPERATIO for double isotopic species (not in
  //		     HITRAN) are products of individual ratios.	 
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "161"	,0.997317	,18.	,11	,11	,TAGS(18003, 18005) ),
	 REC( "181"	,0.00199983	,20.	,12	,12	,TAGS(20003) ),
	 REC( "171"	,0.000372	,19.	,13	,13	,TAGS(19003) ),
	 REC( "162"	,0.00031069	,19.	,14	,14	,TAGS(19002) ),
	 REC( "182"	,621.327e-9	,21.	,-1	,-1	,TAGS(21001) ),
	 REC( "262"	,96.528e-9	,20.	,-1	,-1	,TAGS(20001) )
	 ) ) );

  // CO2 
  // (missing in JPL)
  // Source: 
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // Isotope ratios: HITRAN92 tables.3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO2"),
	DEGFR(2),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "626"	,0.98420	,44.	,21	,21	,TAGS() ),
	 REC( "636"	,0.01106	,45.	,22	,22	,TAGS() ),
	 REC( "628"	,0.0039471	,46.	,23	,23	,TAGS() ),
	 REC( "627"	,0.000734	,45.	,24	,24	,TAGS() ),
	 REC( "638"	,0.00004434	,47.	,25	,25	,TAGS() ),
	 REC( "637"	,0.00000825	,46.	,26	,26	,TAGS() ),
	 REC( "828"	,0.0000039573	,48.	,27	,27	,TAGS() ),
	 REC( "728"	,0.00000147	,47.	,28	,28	,TAGS() )
	 ) ) );
  

  // O3
  // Source: forward_4_96, glob_def.c
  // FIXME: I just copied the number without checking them. This should be checked carefully.
  //	    Do all these HITRAN tags really exist?
  //	    Are JPL tags missing?
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HI-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("666"	,0.9929		,48.	,31	,31	,TAGS(48004)),
	 REC("668"	,4.0E-3		,50.	,32	,32	,TAGS(50004)),
	 REC("686"	,2.0E-3		,50.	,33	,33	,TAGS(50003)),
	 REC("667"	,7.4E-4		,49.	,34	,34	,TAGS(49002)),
	 REC("676"	,3.7E-4		,49.	,35	,35	,TAGS(49001))
	 ) ) );


}


void define_species_map()
{
  extern const ARRAY<SpeciesRecord> species_data;
  extern std::map<string, size_t> SpeciesMap;

  for ( size_t i=0 ; i<species_data.size() ; ++i)
    {
      SpeciesMap[species_data[i].Name()] = i;
    }
}


ostream& operator << (ostream& os, const LineRecord& lr)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  int precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; exit(1);
  }

  os << "@"
     << " " << lr.Name  ()
     << " "
     << setprecision(precision)
     <<        lr.F     ()
     << " " << lr.Psf   ()
     << " " << lr.I0    ()
     << " " << lr.Ti0   ()
     << " " << lr.Elow  ()
     << " " << lr.Agam  ()
     << " " << lr.Sgam  ()
     << " " << lr.Nair  ()
     << " " << lr.Nself ()
     << " " << lr.Tgam  ()
     << " " << lr.Naux  ();

  for ( size_t i=0; i<lr.Naux(); ++i )
    os << " " << lr.Aux()[i];

  return os;
}


/** Extract something from a HITRAN line. This is just a small helper
    function to safe some typing. 

    @param x Output.    What was extracted from the beginning of the line.
    @param line Output. What was extracted is also cut away from line.
    @param n            The width of the stuff to extract.

    @author Stefan Buehler */
template<class T>
void extract(T&      x,
	     string& line,
	     size_t  n)
{
  // This will contain the short substring with the item to extract.
  // Make it a string stream, for easy parsing,
  // extracting substring of width n from line:
  istringstream item( line.substr(0,n) );

//  cout << "line = '" << line << "'\n";
//   cout << "line.substr(0,n) = " << line.substr(0,n) << endl;
//   cout << "item = " << item.str() << endl;

  // Shorten line by n:
  line.erase(0,n);
//  cout << "line = " << line << endl;

  // Convert with the aid of string stream item:
  item >> x;
}


bool LineRecord::ReadFromHitranStream(istream& is)
{
  // Global species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. This is based on the assumption that there are a
  // lot more species than isotopes. Anyway, we also add 100, because
  // at the moment we have only 3 species and more than 3 isotopes.
  const size_t missing = species_data.size() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>]. The value of
  // missing means that we don't have this species.
  //
  // Allow for up to 100 species in HITRAN in the future.
  static ARRAY< size_t >        hspec(100,missing);	

  // This is  an array of arrays for each hitran tag. It contains the
  // ARTS indices of the HITRAN isotopes. 
  static ARRAY< ARRAY<size_t> > hiso(100);

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  // Remember, about which missing species we have already issued a
  // warning: 
  static ARRAY<size_t> warned_missing;

  if ( !hinit )
    {
      for ( size_t i=0; i<species_data.size(); ++i )
	{
	  const SpeciesRecord& sr = species_data[i];

	  // The HITRAN tags are stored as species plus isotope tags
	  // (MO and ISO)
	  // in the Isotope() part of the species record.
	  // We can extract the MO part from any of the isotope tags,
	  // so we use the first one. We do this by taking an integer
	  // division by 10.
	  size_t mo = sr.Isotope()[0].HitranTag() / 10;
//	  cout << "mo = " << mo << endl;
	  hspec[mo] = i; 
	  
	  // Get a nicer to handle array of HITRAN iso tags:
	  size_t n_iso = sr.Isotope().size();
	  ARRAY<int> iso_tags;
	  iso_tags.resize(n_iso);
	  for ( size_t j=0; j<n_iso; ++j )
	    {
	      iso_tags[j] = sr.Isotope()[j].HitranTag();
	    }

	  // Reserve elements for the isotope tags. How much do we
	  // need? This depends on the largest HITRAN tag that we know
	  // about!
	  // Also initialize the tags to missing.
// 	  cout << "iso_tags = " << iso_tags << endl;
// 	  cout << "static_cast<size_t>(max(iso_tags))%10 + 1 = "
// 	       << static_cast<size_t>(max(iso_tags))%10 + 1 << endl;
	  hiso[mo].resize( static_cast<size_t>(max(iso_tags))%10 + 1,
			   missing );
//	  cout << "hiso[mo].size() = " << hiso[mo].size() << endl;

	  // Set the isotope tags:
	  for ( size_t j=0; j<n_iso; ++j )
	    {
	      if ( 0 < iso_tags[j] )				  // ignore -1 elements
		{
		  // To get the iso tags from HitranTag() we also have to take
		  // modulo 10 to get rid of mo.
//		  cout << "iso_tags[j] % 10 = " << iso_tags[j] % 10 << endl;
		  hiso[mo][iso_tags[j] % 10] = j;
		}
	    }
	}
      
//      cout << "hiso = " << hiso << endl << "***********" << endl;


      // Print the generated data structures (for debugging):
      out3 << "  HITRAN index table:\n";
      for ( size_t i=0; i<hspec.size(); ++i )
	{
	  if ( missing != hspec[i] )
	    {
	      // The explicit conversion of Name to a c-string is
	      // necessary, because setw does not work correctly for
	      // stl strings.
	      out3 << "  mo = " << i << "   Species = "
		   << setw(10) << setiosflags(ios::left)
		   << species_data[hspec[i]].Name().c_str()
		   << "iso = ";
	      for ( size_t j=1; j<hiso[i].size(); ++j )
		{
		  if ( missing==hiso[i][j] )
		    out3 << " " << "m";
		  else
		    out3 << " " << species_data[hspec[i]].Isotope()[hiso[i][j]].Name();
		}
	      out3 << "\n";
	    }
	}

      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  string line;

  // The first item is the molecule number:
  size_t mo;

  // Look for more comments?
  bool comment = true;

  while (comment)
    {
      // Return true if eof is reached:
      if (is.eof()) return true;

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract molecule number:
      mo = 0;
      // Initialization of mo is important, because mo stays the same
      // if line is empty.
      extract(mo,line,2);
      //      cout << "mo = " << mo << endl;
  
      // If mo == 0 this is just a comment line:
      if ( 0 != mo )
	{
	  // See if we know this species. At the moment, because I
	  // only implemented 3 species, we will only give a warning
	  // if a species is unknown, and not an error. FIXME: This should
	  // probably be changed in the future.
	  if ( missing != hspec[mo] )	    comment = false ;
	  else
	    {
	      // See if this is already in warned_missing, use
	      // std::count for that:
	      if ( 0 == std::count(warned_missing.begin(),
				   warned_missing.end(),
				   mo) )
		{
		  out0 << "Warning: HITRAN mo = " << mo << " is not "
		       << "known to ARTS.\n";
		  warned_missing.push_back(mo);
		}
	    }
	}
    }

  // Ok, we seem to have a valid species here.

  // Set mspecies from my cool index table:
  mspecies = hspec[mo];

  // Extract isotope:
  size_t iso;				  
  extract(iso,line,1);
  //  cout << "iso = " << iso << endl;


  // Set misotope from the other cool index table.
  // We have to be careful to issue an error for unknown iso tags. Iso
  // could be either larger than the size of hiso[mo], or set
  // explicitly to missing. Unfortunately we have to test both cases. 
  misotope = missing;
  if ( iso < hiso[mo].size() )
    if ( missing != hiso[mo][iso] )
      misotope = hiso[mo][iso];

  // Issue error message if misotope is still missing:
  if (missing == misotope)
    {
      ostringstream os;
      os << "Species: " << species_data[mspecies].Name()
	 << ", isotope iso = " << iso
	 << " is unknown.";
      throw runtime_error(os.str());
    }

  
  // Position.
  {
    // HITRAN position in wavenumbers (cm^-1):
    Numeric v;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a line
    // position in wavenumber (cm^-1) by this constant, you get the
    // frequency in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 100.;

    // Extract HITRAN postion:
    extract(v,line,12);

    // ARTS position in Hz:
    mf = v * w2Hz;
//    cout << "mf = " << mf << endl;
  }

  // Intensity.
  {
    // HITRAN intensity in cm-1/(molec * cm-2) at 296 Kelvin.
    // Includes isotpic ratio!
    Numeric s;
    // Conversion from HITRAN intensities to JPL intensities 
    // (nm^2 MHz). This factor is given in the JPL documentation. 
    const Numeric hi2jpl = 2.99792458e18;
    // Because we want SI units m and Hz instead of nm and MHz, we
    // need to shift the decimals a bit.
    // More importantly, we have to take out the isotope ratio, which
    // is included in the HITRAN intensity.
    const Numeric hi2arts = hi2jpl * 1e-12 / IsotopeData().Abundance();

    // Extract HITRAN intensity:
    extract(s,line,10);

    // ARTS intensity in m^2/Hz:
    mi0 = s * hi2arts;
//    cout << "mi0 = " << mi0 << endl;
  }

  
  // Skip transition probability:
  {
    Numeric r;
    extract(r,line,10);
  }
  

  // Air broadening parameters.
  {
    // HITRAN parameter is in cm-1/atm at 296 Kelvin
    // All parameters are HWHM (I hope this is true!)
    Numeric gam;
    // External constant from constants.cc: Converts atm to
    // hPa. Multiply value in atm by this number to get value in hPa. 
    extern const Numeric ATM2HPA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to MHz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in MHz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e-4;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2HPA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in MHz/hPa:
    magam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in MHz/hPa:
    msgam = gam * hi2arts;

//    cout << "agam, sgam = " << magam << ", " << msgam << endl;
  }


  // Lower state energy.
  {
    // HITRAN parameter is in wavenumbers (cm^-1).
    // The same unit as in ARTS, we can extract directly!
    extract(melow,line,10);
//    cout << "mf, melow = " << mf << ", " << setprecision(20) << melow << endl;
  }

  
  // Temperature coefficient of broadening parameters.
  {
    // This is dimensionless, we can also extract directly.
    extract(mnair,line,4);

    // Set self broadening temperature coefficient to the same value:
    mnself = mnair;
//    cout << "mnair = " << mnair << endl;
  }


  // Pressure shift.
  {
    // HITRAN value in cm^-1
    // (I guess they really mean cm^-1 / atm. Anyway, that's what I'll
    // assume. So the conversion goes exactly as for the broadening
    // parameters. 
    Numeric d;
    // External constant from constants.cc: Converts atm to
    // hPa. Multiply value in atm by this number to get value in hPa. 
    extern const Numeric ATM2HPA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to MHz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in MHz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e-4;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2HPA;

    // Extract HITRAN value:
    extract(d,line,8);

    // ARTS value in MHz/hPa
    mpsf = d * hi2arts;
  }


  // These were all the parameters that we can extract from
  // HITRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for HITRAN)
  mti0 = 296.0;					  

  // Reference temperature for AGAM and SGAM in K.
  // (This is also fix for HITRAN)
  mtgam = 296.0;


  // That's it!
  return false;
}

OneTag::OneTag(string def) 
{
  // Species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;
  // The species map. This is used to find the species id.
  extern std::map<string, size_t> SpeciesMap;
  // Name of species and isotope (aux variables):
  string name, isoname;
  // Aux index:
  size_t n;


  // Set frequency limits to default values (no limits):
  mlf = -1;
  muf = -1;

  // We cannot set a default value for the isotope, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species. 
    

  // Extract the species name:
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      name = def.substr(0,n);      // Extract before '-'
      def.erase(0,n+1);		    // Remove from def
    }
  else
    {
      // n==def.size means that def does not contain a '-'. In that
      // case we assume that it contains just a species name and
      // nothing else
      name = def;
      def  = "";
    }


//  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<string, size_t>::const_iterator mi = SpeciesMap.find(name);
  if ( mi != SpeciesMap.end() )
    {
      // Ok, we've found the species. Set mspecies.
      mspecies = mi->second;
    }
  else
    {
      ostringstream os;
      os << "Species " << name << " is not a valid species.";
      throw runtime_error(os.str());
    }

  // Get a reference to the relevant Species Record:
  const SpeciesRecord& spr = species_data[mspecies];

  if ( 0 == def.size() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopes and no frequency limits.
      // Frequency defaults are already set. Set isotope defaults:
      misotope = spr.Isotope().size();
      // This means all isotopes.

      return;
    }
    
  // Extract the isotope name:
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      isoname = def.substr(0,n);          // Extract before '-'
      def.erase(0,n+1);		    // Remove from def
    }
  else
    {
      // n==def.size means that def does not contain a '-'. In that
      // case we assume that it contains just the isotope name and
      // nothing else.
      isoname = def;
      def  = "";
    }
    
  // Check for joker:
  if ( "*" == isoname )
    {
      // The user wants all isotopes. Set this accordingly:
      misotope = spr.Isotope().size();
    }
  else
    {
      // Make an array containing the isotope names:
      ARRAY<string> ins;
      for ( size_t i=0; i<spr.Isotope().size(); ++i )
	ins.push_back( spr.Isotope()[i].Name() );

      // Use the find algorithm from the STL to find the isotope ID. It
      // returns an iterator, so to get the index we take the
      // difference to the begin() iterator.
      misotope = find( ins.begin(),
		       ins.end(),
		       isoname ) - ins.begin();

      // Check if we found a matching isotope:
      if ( misotope >= ins.size() ) 
	{
	  ostringstream os;
	  os << "Isotope " << isoname << " is not a valid isotope for "
	     << "species " << name << ".\n"
	     << "Valid isotopes are:";
	  for ( size_t i=0; i<ins.size(); ++i )
	    os << " " << ins[i];
	  throw runtime_error(os.str());
	}
    }

  if ( 0 == def.size() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants no frequency limits.  Frequency defaults are
      // already set, so we can return directly.

      return;
    }


  // Look for the two frequency limits:
  
  // Extract first frequency
  n    = def.find('-');    // find the '-'
  if (n < def.size() )
    {
      // Frequency as a string:
      string fname;
      fname = def.substr(0,n);              // Extract before '-'
      def.erase(0,n+1);		      // Remove from def

      // Check for joker:
      if ( "*" == fname )
	{
	  // The default for mlf is already -1, meaning `ALL'.
	  // So there is nothing to do here.
	}
      else
	{
	  // Convert to Numeric:
	  istringstream is(fname);
	  is >> mlf;
	}
    }
  else
    {
      // n==def.size means that def does not contain a '-'. In this
      // case that is not allowed!
      throw runtime_error("You must either speciefy both frequency limits\n"
			  "(at least with jokers), or none.");
    }


  // Now there should only be the upper frequency left in def.
  // Check for joker:
  if ( "*" == def )
    {
      // The default for muf is already -1, meaning `ALL'.
      // So there is nothing to do here.
    }
  else
    {
      // Convert to Numeric:
      istringstream is(def);
      is >> muf;
    }
}


string OneTag::Name() const 
{
  // Species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;
  // A reference to the relevant record of the species data:
  const  SpeciesRecord& spr = species_data[mspecies];
  // For return value:
  ostringstream os;

  // First the species name:
  os << spr.Name() << "-";

  // Now the isotope. Can be a single isotope or ALL.
  if ( misotope == spr.Isotope().size() )
    {
      // One larger than allowed means all isotopes!
      os << "*-";
    }
  else
    {
      os << spr.Isotope()[misotope].Name() << "-";
    }

  // Now the frequency limits, if there are any. For this we first
  // need to determine the floating point precision.

  // Determine the precision, depending on whether Numeric is double
  // or float:  
  int precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; exit(1);
  }

  if ( 0 > mlf )
    {
      // mlf < 0 means no lower limit.
      os << "*-";
    }
  else
    {
      os << setprecision(precision);
      os << mlf << "-";
    }

  if ( 0 > muf )
    {
      // muf < 0 means no upper limit.
      os << "*";
    }
  else
    {
      os << setprecision(precision);
      os << muf;
    }

  return os.str();
}

ostream& operator << (ostream& os, const OneTag& ot)
{
  return os << ot.Name();
}


void write_lines_to_stream(ostream& os,
			   const ARRAYofLineRecord& lines)
{
  for ( size_t i=0; i<lines.size(); ++i )
    {
      os << lines[i] << '\n';
    }
}
