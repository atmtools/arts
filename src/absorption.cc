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

/*!
  \file   absorption.cc
  \brief  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  \author Stefan Buehler
*/


#include "arts.h"
#include "absorption.h"
#include "math_funcs.h"
#include "messages.h"


/*! The map associated with species_data. */
std::map<string, size_t> SpeciesMap;


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

	  // We have to be careful and check for the case that all
	  // HITRAN isotope tags are -1 (this species is missing in HITRAN).

	  if ( 0 < sr.Isotope()[0].HitranTag() )
	    {
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
		  out0 << "Error: HITRAN mo = " << mo << " is not "
		       << "known to ARTS.\n";
		  warned_missing.push_back(mo);
		  exit(1);
		  // SAB 08.08.2000 If you want to make the program
		  // continue anyway, just comment out the exit
		  // line.
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
    // Pa. Multiply value in atm by this number to get value in Pa. 
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    magam = gam * hi2arts;

    // Extract HITRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
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
    // Pa. Multiply value in atm by this number to get value in Pa. 
    extern const Numeric ATM2PA;
    // External constant from constants.cc:
    extern const Numeric SPEED_OF_LIGHT;
    // Conversion from wavenumber to Hz. If you multiply a value in
    // wavenumber (cm^-1) by this constant, you get the value in Hz.
    const Numeric w2Hz = SPEED_OF_LIGHT * 1e2;
    // Ok, put together the end-to-end conversion that we need:
    const Numeric hi2arts = w2Hz / ATM2PA;

    // Extract HITRAN value:
    extract(d,line,8);

    // ARTS value in Hz/Pa
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

/** The Lorentz line shape. This is a quick and dirty implementation.

    @param F     Output: The shape function.
    @param f0    Line center frequency.
    @param f_abs The frequency grid.
    @param gamma The pressure broadening parameter.

    @author Stefan Buehler 16.06.2000 */
void line_shape_lorentz(VECTOR&       ls,
			Numeric	      f0,
			Numeric       gamma,
			const VECTOR& f_abs)
{
  // FIXME: Maybe try if call by reference is faster for f0 and gamma?

  // 1/PI:
  extern const Numeric PI;
  static const Numeric invPI = 1. / PI;

  assert( ls.dim() == f_abs.dim() );

  for ( size_t i=0; i<f_abs.dim(); ++i )
    {
      ls[i] = invPI * gamma / ( pow( f_abs[i]-f0, 2) + pow(gamma,2) );
    }
}

void abs_species( MATRIX&                  abs,
		  const VECTOR&  	   f_abs,
		  const VECTOR&  	   p_abs,
		  const VECTOR&  	   t_abs,           
		  const VECTOR&            vmr,
		  const ARRAYofLineRecord& lines )
{
  // Define the vector for the line shape function here, so that we
  // don't need so many free store allocations.
  VECTOR ls(f_abs.dim());

  // Check that p_abs, t_abs, and vmr all have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error. 

  if ( t_abs.dim() != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable t_abs must have the same dimension as p_abs.\n"
	 << "t_abs.dim() = " << t_abs.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  if ( vmr.dim() != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable vmr must have the same dimension as p_abs.\n"
	 << "vmr.dim() = " << vmr.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  // Check that the dimension of abs is indeed [f_abs.dim(),
  // p_abs.dim()]:
  if ( abs.dim(1) != f_abs.dim() || abs.dim(2) != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable abs must have dimensions [f_abs.dim(),p_abs.dim()].\n"
	 << "[abs.dim(1),abs.dim(2)] = [" << abs.dim(1)
	 << ", " << abs.dim(2) << "]\n"
	 << "f_abs.dim() = " << f_abs.dim() << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  // Loop all pressures:
  for ( size_t i=0; i<p_abs.dim(); ++i )
  {
    out3 << "  p = " << p_abs[i] << " Pa\n";
    // Loop all lines:
    for ( size_t l=0; l<lines.dim(); ++l )
      {
	// 1. Convert intensity to units of [ m^-1 * Hz / Pa ]
	Numeric intensity = lines[l].I0();
	// FIXME: Move the constants here to constants.cc?
	  {
	    // I have no time to find out about the correct intensity
	    // conversion. So, as a quick hack, we will convert to
	    // MYTRAN intensiy (HITRAN, except isotopic ratio is taken
	    // out, and then use the old conversion from abs_my.c. See
	    // LineRecord::ReadFromHitranStream for the below two
	    // conversions. 

	    // Conversion from HITRAN intensities to JPL intensities 
	    // (nm^2 MHz). This factor is given in the JPL documentation. 
	    const Numeric hi2jpl = 2.99792458e18;
	    // Because we want SI units m and Hz instead of nm and MHz, we
	    // need to shift the decimals a bit.
	    // More importantly, we have to take out the isotope ratio, which
	    // is included in the HITRAN intensity.
	    const Numeric hi2arts = hi2jpl * 1e-12;

	    intensity = intensity / hi2arts;

	    // intensity is now MYTRAN intensity.

	    // From abs_my.c:
	    const Numeric T_0_C = 273.15;  	       /* temp. of 0 Celsius in [K]  */
	    const Numeric p_0   = 101300.25; 	       /* standard p in [Pa]        */
	    const Numeric n_0   = 2.686763E19;         /* Loschmidt constant [cm^-3] */
	    const Numeric c     = 29.9792458;          /*  	     [ GHz / cm^-1 ] */

	    intensity = intensity * n_0 * T_0_C / (p_0*t_abs[i]);
	    /* n = n0*T0/p0 * p/T, ideal gas law
	       ==> [intensity] = [cm^2 / Pa] */

	    intensity = intensity * c; 			/* [cm^-1 * GHz / Pa] */
	    intensity = intensity * 1E5; 		/* [km^-1 * GHz / Pa] */
	    intensity = intensity * 1e6; 		/* [m^-1 * Hz / Pa] */
	  }

	// FIXME: Are the units of the intensity correct?
	
	// FIXME: We should calculate the temperature dependence of the
	// intensities here.

	// Scale with isotopic ratio. FIXME: This is inefficient. The
	// scaled intensities could be stored in lines.
	intensity = intensity * lines[l].IsotopeData().Abundance();

	// 2. Get pressure broadened line width:
	// (Agam is in Hz/Pa, p_abs is in Pa, f_abs is in Hz,
	// gamma is in Hz/Pa)
	// FIXME: This is inefficient, the scaling to Hz could be
	// stored in lines.

	const Numeric theta = lines[l].Tgam() / t_abs[i];

	const Numeric p_partial = p_abs[i] * vmr[i];

	Numeric gamma
	  = lines[l].Agam() * pow(theta, lines[l].Nair())  * (p_abs[i] - p_partial)
	  + lines[l].Sgam() * pow(theta, lines[l].Nself()) * p_partial;

	// Ignore Doppler broadening for now. FIXME: Add this.

	// Get line shape:
	line_shape_lorentz(ls,
			   lines[l].F(),
			   gamma,
			   f_abs);

	// Add line to abs:
	for ( size_t j=0; j<f_abs.dim(); ++j )
	  {
	    abs[j][i] = abs[j][i]
	      + p_partial * intensity * ls[j];
	  }
      }
  }
}
