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

/**
  \file   absorption.cc

  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  \author Stefan Buehler
*/


#include "arts.h"
#include "absorption.h"
#include "math_funcs.h"
#include "md.h"
#include "messages.h"


/** The map associated with species_data. */
std::map<string, size_t> SpeciesMap;


// member fct of isotoperecord, calculates the partition fct at the
// given temperature  from the partition fct coefficients (3rd order
// polynomial in temperature)
Numeric IsotopeRecord::CalculatePartitionFctAtTemp( Numeric
						    temperature ) const
{
  Numeric result = 0.;
  Numeric exponent = 1.;

  ARRAY<Numeric>::const_iterator it;

  for (it=mqcoeff.begin(); it != mqcoeff.end(); it++)
    {
      result += *it * exponent;
      exponent *= temperature;
    }
  return result;
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


/** Extract something from a catalogue line. This is just a small helper
    function to safe some typing. 

    \retval x    What was extracted from the beginning of the line.
    \retval line What was extracted is also cut away from line.
    \param n     The width of the stuff to extract.

    \author Stefan Buehler */
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
  extern ARRAY<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotope.
  const size_t missing = species_data.size() + 100;

  // We need a species index sorted by HITRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<HITRAN tag>]. 
  //
  // Allow for up to 100 species in HITRAN in the future.
  static ARRAY< size_t >        hspec(100);

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
      // Initialize hspec.
      // The value of missing means that we don't have this species.
      mtl::set(hspec,missing);

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
	      hiso[mo] = ARRAY<size_t>( max(iso_tags)%10 + 1 );
	      mtl::set(hiso[mo], missing);

	      // Set the isotope tags:
	      for ( size_t j=0; j<n_iso; ++j )
		{
		  if ( 0 < iso_tags[j] )				  // ignore -1 elements
		    {
		      // To get the iso tags from HitranTag() we also have to take
		      // modulo 10 to get rid of mo.
		      hiso[mo][iso_tags[j] % 10] = j;
		    }
		}
	    }
	}


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
	  // See if we know this species. Exit with an error if the species is unknown. 
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
    extern const Numeric SPEED_OF_LIGHT; // in [m/s]

    // HITRAN intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // It already includes the isotpic ratio.
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function). 
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
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

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );

  // Reference temperature for AGAM and SGAM in K.
  // (This is also fix for HITRAN)
  mtgam = 296.0;


  // That's it!
  return false;
}


bool LineRecord::ReadFromMytran2Stream(istream& is)
{
  // Global species lookup data:
  extern ARRAY<SpeciesRecord> species_data;

  // This value is used to flag missing data both in species and
  // isotope lists. Could be any number, it just has to be made sure
  // that it is neither the index of a species nor of an isotope.
  const size_t missing = species_data.size() + 100;

  // We need a species index sorted by MYTRAN tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is hind[<MYTRAN tag>]. The value of
  // missing means that we don't have this species.
  //
  // Allow for up to 100 species in MYTRAN in the future.
  static ARRAY< size_t >        hspec(100,missing);	

  // This is  an array of arrays for each mytran tag. It contains the
  // ARTS indices of the MYTRAN isotopes. 
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
	  // MYTRAN isotope tags are -1 (this species is missing in MYTRAN).

	  if ( 0 < sr.Isotope()[0].MytranTag() )
	    {
	      // The MYTRAN tags are stored as species plus isotope tags
	      // (MO and ISO)
	      // in the Isotope() part of the species record.
	      // We can extract the MO part from any of the isotope tags,
	      // so we use the first one. We do this by taking an integer
	      // division by 10.
	  
	      size_t mo = sr.Isotope()[0].MytranTag() / 10;
	      //	  cout << "mo = " << mo << endl;
	      hspec[mo] = i; 
	  
	      // Get a nicer to handle array of MYTRAN iso tags:
	      size_t n_iso = sr.Isotope().size();
	      ARRAY<int> iso_tags;
	      iso_tags.resize(n_iso);
	      for ( size_t j=0; j<n_iso; ++j )
		{
		  iso_tags[j] = sr.Isotope()[j].MytranTag();
		}

	      // Reserve elements for the isotope tags. How much do we
	      // need? This depends on the largest MYTRAN tag that we know
	      // about!
	      // Also initialize the tags to missing.
	      // 	  cout << "iso_tags = " << iso_tags << endl;
	      // 	  cout << "static_cast<size_t>(max(iso_tags))%10 + 1 = "
	      // 	       << static_cast<size_t>(max(iso_tags))%10 + 1 << endl;
	      hiso[mo] = ARRAY<size_t>( max(iso_tags)%10 + 1 );
	      mtl::set(hiso[mo], missing);

	      // Set the isotope tags:
	      for ( size_t j=0; j<n_iso; ++j )
		{
		  if ( 0 < iso_tags[j] )				  // ignore -1 elements
		    {
		      // To get the iso tags from MytranTag() we also have to take
		      // modulo 10 to get rid of mo.
		      //		  cout << "iso_tags[j] % 10 = " << iso_tags[j] % 10 << endl;
		      hiso[mo][iso_tags[j] % 10] = j;
		    }
		}
	    }
	}
      
//      cout << "hiso = " << hiso << endl << "***********" << endl;


      // Print the generated data structures (for debugging):
      out3 << "  MYTRAN index table:\n";
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
      if (is.eof()) {
	//	cout << "Eof" << endl;
	return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);
      //      cout << line << endl;

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract molecule number:
      mo = 0;
      // Initialization of mo is important, because mo stays the same
      // if line is empty.
      extract(mo,line,2);
      //           cout << "mo = " << mo << endl;
  
      // If mo == 0 this is just a comment line:
      if ( 0 != mo )
	{
	  // See if we know this species. We will give a warning
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
		  out0 << "Error: MYTRAN mo = " << mo << " is not "
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
    // MYTRAN position in MHz:
    Numeric v;

    // Extract MYTRAN postion:
    extract(v,line,13);

    // ARTS position in Hz:
    mf = v * 1E6;
//    cout << "mf = " << mf << endl;
  }

  // skip transition frequency error
  {
    Numeric r;
    extract(r,line,8);
  } 
  
  // Intensity.
  {
    extern const Numeric SPEED_OF_LIGHT; // in [m/s]

    // MYTRAN2 intensity is in cm-1/(molec * cm-2) at 296 Kelvin.
    // (just like HITRAN, only isotopic ratio is not included)
    // The first cm-1 is the frequency unit (it cancels with the
    // 1/frequency unit of the line shape function). 
    //
    // We need to do the following:
    // 1. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 2. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)
    // 3. Multiply with the isotopic ratio

    const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;

    Numeric s;

    // Extract HITRAN intensity:
    extract(s,line,10);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s * hi2arts;

    // multiply with the isotopic ratio
    mi0 *= species_data[mspecies].Isotope()[misotope].Abundance();
  }

  
  // Air broadening parameters.
  {
    // MYTRAN parameter is in MHz/Torr at reference temperature
    // All parameters are HWHM
    Numeric gam;
    // External constant from constants.cc: Converts torr to
    // Pa. Multiply value in torr by this number to get value in Pa. 
    extern const Numeric TORR2PA;

    // Extract HITRAN AGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    magam = gam * 1E6 / TORR2PA;

    // Extract MYTRAN SGAM value:
    extract(gam,line,5);

    // ARTS parameter in Hz/Pa:
    msgam = gam * 1E6 / TORR2PA;

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

    // Extract the self broadening parameter:
    extract(mnself,line,4);
//    cout << "mnair = " << mnair << endl;
  }

  
  // Reference temperature for broadening parameter in K:
  {
    // correct units, extract directly
    extract(mtgam,line,7);
  }


  // Pressure shift.
  {
    // MYTRAN value in MHz/Torr
    Numeric d;
    // External constant from constants.cc: Converts torr to
    // Pa. Multiply value in torr by this number to get value in Pa. 
    extern const Numeric TORR2PA;

    // Extract MYTRAN value:
    extract(d,line,8);

    // ARTS value in Hz/Pa
    mpsf = d * 1E6 / TORR2PA;
  }


  // These were all the parameters that we can extract from
  // MYTRAN. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for MYTRAN2)
  mti0 = 296.0;

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );



  // That's it!
  return false;
}


bool LineRecord::ReadFromJplStream(istream& is)
{
  // Global species lookup data:
  extern ARRAY<SpeciesRecord> species_data;

  // We need a species index sorted by JPL tag. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is JplMap[<JPL tag>]. We need int in this map,
  // because the tag array within the species data is an int array.
  static map<int, SpecIsoMap> JplMap;

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  if ( !hinit )
    {

      out3 << "  JPL index table:\n";

      for ( size_t i=0; i<species_data.size(); ++i )
	{
	  const SpeciesRecord& sr = species_data[i];


	  for ( size_t j=0; j<sr.Isotope().size(); ++j)
	    {
	      
	      for (size_t k=0; k<sr.Isotope()[j].JplTags().size(); ++k)
		{

		  SpecIsoMap indicies(i,j);

		  JplMap[sr.Isotope()[j].JplTags()[k]] = indicies;

		  // Print the generated data structures (for debugging):
		  // The explicit conversion of Name to a c-string is
		  // necessary, because setw does not work correctly for
		  // stl strings.
		  const size_t& i1 = JplMap[sr.Isotope()[j].JplTags()[k]].Speciesindex();
		  const size_t& i2 = JplMap[sr.Isotope()[j].JplTags()[k]].Isotopeindex();
					 
		  out3 << "  JPL TAG = " << sr.Isotope()[j].JplTags()[k] << "   Species = "
		       << setw(10) << setiosflags(ios::left)
		       << species_data[i1].Name().c_str()
		       << "iso = " 
		       << species_data[i1].Isotope()[i2].Name().c_str()
		       << "\n";
		}

	    }
	}
      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  string line;







  // Look for more comments?
  bool comment = true;

  while (comment)
    {
      // Return true if eof is reached:
      if (is.eof()) {
	//	cout << "Eof" << endl;
	return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);

      // Because of the fixed FORTRAN format, we need to break up the line
      // explicitly in apropriate pieces. Not elegant, but works!

      // Extract center frequency:
      // Initialization of v is important, because v stays the same
      // if line is empty.
      // JPL position in MHz:
      Numeric v = 0.0;
      
      // Extract JPL position:
      extract(v,line,13);
	
      // check for empty line
      if (v != 0.0)
	{
	  // ARTS position in Hz:
	  mf = v * 1E6;
	  
	  comment = false;
	}
    }


  // skip transition frequency error
  {
    Numeric r;
    extract(r,line,8);
  } 
  
  // Intensity.
  {
    // JPL has log (10) of intensity in nm2 MHz at 300 Kelvin.
    //
    // We need to do the following:
    // 1. take 10^intensity 
    // 2. convert to cm-1/(molecule * cm-2): devide by c * 1e10
    // 3. Convert frequency from wavenumber to Hz (factor 1e2 * c)
    // 4. Convert [molec * cm-2] to [molec * m-2] (factor 1e-4)
    // 5. Multiply with the isotopic ratio (done later after the tag 
    //    is extracted)

    Numeric s;

    // Extract JPL intensity:
    extract(s,line,8);

    // remove log
    s = pow(10,s);

    // Convert to ARTS units (Hz / (molec * m-2) ), or shorter: Hz*m^2
    mi0 = s / 1E12;
  }

  // Degrees of freedom
  {
    size_t dr;

    // Extract degrees of freedom
    extract(dr,line,2);
  }

  // Lower state energy.
  {
    // JPL parameter is in wavenumbers (cm^-1).
    // The same unit as in ARTS, we can extract directly!
    extract(melow,line,10);
  }

  // Upper state degeneracy
  {
    size_t gup;

    // Extract upper state degeneracy
    extract(gup,line,3);
  }

  // Tag number
  int tag;
  {
    // Extract Tag number
    extract(tag,line,7);

    // make sure tag is not negative (damned jpl cat):
    tag = abs(tag);
  }

  // ok, now for the cool index map:

  // is this tag valid?
  const map<int, SpecIsoMap>::const_iterator i = JplMap.find(tag);
  if ( i == JplMap.end() )
    {
      ostringstream os;
      os << "JPL Tag: " << tag << " is unknown.";
      throw runtime_error(os.str());
    }

  SpecIsoMap id = i->second;


  // Set mspecies:
  mspecies = id.Speciesindex();

  // Set misotope:
  misotope = id.Isotopeindex();

  // multiply intensity with the isotopic ratio
  mi0 *= species_data[mspecies].Isotope()[misotope].Abundance();

  // Air broadening parameters: unknown to jpl, use old iup forward
  // model default values, which is mostly set to 0.0025 GHz/hPa, even
  // though for some lines the pressure broadening is given explicitly
  // in the program code. The explicitly given values are ignored and
  // only the default value is set. Self broadening was in general not
  // considered in the old forward model.
  {
    // ARTS parameter in Hz/Pa:
    magam = 2.5E4;

    // ARTS parameter in Hz/Pa:
    msgam = 0.0;
  }


  // Temperature coefficient of broadening parameters. Was set to 0.75
  // in old forward model, even though for some lines the parameter is
  // given explicitly in the program code. The explicitly given values
  // are ignored and only the default value is set. Self broadening
  // not considered.
  {
    mnair = 0.75;
    mnself = 0.0;
  }

  
  // Reference temperature for broadening parameter in K, was
  // generally set to 300 K in old forward model, with the exceptions
  // as already mentioned above:
  {
    mtgam = 300.0;
  }


  // Pressure shift: not given in JPL, set to 0
  {
    mpsf = 0.0;
  }


  // These were all the parameters that we can extract from
  // JPL. However, we still have to set the reference temperatures
  // to the appropriate value:

  // Reference temperature for Intensity in K.
  // (This is fix for JPL)
  mti0 = 300.0;

  // calculate the partition fct at the reference temperature. This is
  // done for all lines read from the catalogue, even the ones never
  // used in calculation. FIXME: are there better ways to implement this?
  species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );



  // That's it!
  return false;
}


bool LineRecord::ReadFromArtsStream(istream& is)
{
  // Global species lookup data:
  extern ARRAY<SpeciesRecord> species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts string>]. 
  static map<string, SpecIsoMap> ArtsMap;

  // Remeber if this stuff has already been initialized:
  static bool hinit = false;

  if ( !hinit )
    {

      out3 << "  ARTS index table:\n";

      for ( size_t i=0; i<species_data.size(); ++i )
	{
	  const SpeciesRecord& sr = species_data[i];


	  for ( size_t j=0; j<sr.Isotope().size(); ++j)
	    {
	      
	      SpecIsoMap indicies(i,j);
	      string buf = sr.Name()+"-"+sr.Isotope()[j].Name();
	      
	      ArtsMap[buf] = indicies;
	      
	      // Print the generated data structures (for debugging):
	      // The explicit conversion of Name to a c-string is
	      // necessary, because setw does not work correctly for
	      // stl strings.
	      const size_t& i1 = ArtsMap[buf].Speciesindex();
	      const size_t& i2 = ArtsMap[buf].Isotopeindex();
					 
	      out3 << "  Arts Identifier = " << buf << "   Species = "
		   << setw(10) << setiosflags(ios::left)
		   << species_data[i1].Name().c_str()
		   << "iso = " 
		   << species_data[i1].Isotope()[i2].Name().c_str()
		   << "\n";

	    }
	}
      hinit = true;
    }


  // This contains the rest of the line to parse. At the beginning the
  // entire line. Line gets shorter and shorter as we continue to
  // extract stuff from the beginning.
  string line;


  // Look for more comments?
  bool comment = true;

  while (comment)
    {
      // Return true if eof is reached:
      if (is.eof()) {
	//	cout << "Eof" << endl;
	return true;
      }

      // Throw runtime_error if stream is bad:
      if (!is) throw runtime_error ("Stream bad.");

      // Read line from file into linebuffer:
      getline(is,line);

      // @ as first character marks catalogue entry
      char c;
      extract(c,line,1);
      
      // check for empty line
      if (c == '@')
	{
	  comment = false;
	}
    }


  // read the arts identifier string
  istringstream icecream(line);

  string artsid;
  icecream >> artsid;

  if (artsid.length() != 0)
    {

      // ok, now for the cool index map:
      // is this arts identifier valid?
      const map<string, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
      if ( i == ArtsMap.end() )
	{
	  ostringstream os;
	  os << "ARTS Tag: " << artsid << " is unknown.";
	  throw runtime_error(os.str());
	}

      SpecIsoMap id = i->second;


      // Set mspecies:
      mspecies = id.Speciesindex();

      // Set misotope:
      misotope = id.Isotopeindex();


      // Extract center frequency:
      icecream >> mf;


      // Extract pressure shift:
      icecream >> mpsf;
  
      // Extract intensity:
      icecream >> mi0;


      // Extract reference temperature for Intensity in K:
      icecream >> mti0;


      // Extract lower state energy:
      icecream >> melow;


      // Extract air broadening parameters:
      icecream >> magam;
      icecream >> msgam;

      // Extract temperature coefficient of broadening parameters:
      icecream >> mnair;
      icecream >> mnself;

  
      // Extract reference temperature for broadening parameter in K:
      icecream >> mtgam;

      // Extract the aux parameters:
      size_t naux;
      icecream >> naux;

      // resize the aux array and read it
      maux.resize(naux);

      for (size_t i = 0; i<naux; i++)
	{
	  icecream >> maux[i];
	}


      // calculate the partition fct at the reference temperature. This is
      // done for all lines read from the catalogue, even the ones never
      // used in calculation. FIXME: are there better ways to implement this?
      species_data[mspecies].Isotope()[misotope].CalculatePartitionFctAtRefTemp( mti0 );

    }

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



/** 
   Returns the index among some tag groups for an array of tag strings. 
   
   \begin{verbatim}
   For example, if tags1 correspond to the definition
     ["O3","H2O-161,H2O-162"]
   and the tag strings are
     ["H2O-161,H2O-162","O3"]
   the tags1_index becomes
     [2,1]
   \end{verbatim}

   @exception runtime_error  Some string is not a valid tag item.
   @exception runtime_error  Not all strings are not found among the tags.

   \retval tags1_index     Index in tags1 for tags2_strings
   \param  tags1           The tags to search in.
   \param  tags2_strings   The tag strings for which indeces shall be found.

   \author Patrick Eriksson 
   \date 2000-12-06
*/
void get_tagindex_for_strings( 
              ARRAYofsizet&   tags1_index, 
        const TagGroups&      tags1, 
        const ARRAYofstring&  tags2_strings )
{
  const size_t   n1 = tags1.size();
  const size_t   n2 = tags2_strings.size();
     TagGroups   tags2;                // Internal tag names for tag_strings
        size_t   i1, i2, nj, j, found, ok;

  tags1_index.resize(n2);
  tag_groupsDefine( tags2, tags2_strings );

  for ( i2=0; i2<n2; i2++ )
  {
    found = 0;
    for ( i1=0; (i1<n1) && !found; i1++ )
    {
      nj = tags2[i2].size(); 
      if ( nj  == tags1[i1].size() )
      {
        ok = 1;
        for ( j=0; j<nj; j++ )
	{
          if ( tags2[i2][j].Name() != tags1[i1][j].Name() )
            ok = 0;
	}
        if ( ok )
	{
          found = 1;
          tags1_index[i2] = i1;
	}
      }
    } 
    if ( !found )
    {
      ostringstream os;
      os << "The tag string \"" << tags2_strings[i2] << 
            "\" does not match any of the given tags.\n";
      throw runtime_error(os.str());
    }
  }
}


/** A helper function that writes lines in a line list to a stream. 

    \retval os      The stream to write to.
    \param  lines   The line list to write.

    \author Stefan Buehler 
    \date 2000-06-12 */
void write_lines_to_stream(ostream& os,
			   const ARRAYofLineRecord& lines)
{
  for ( size_t i=0; i<lines.size(); ++i )
    {
      os << lines[i] << '\n';
    }
}


/** Calculate absorption coefficients for one tag group. All lines in
    the line list must belong to the same species. This must be
    ensured by lines_per_tgCreateFromLines, so it is only verified
    with assert. Also, the input vectors p_abs, t_abs, and vmr must
    all have the same dimension.

    \retval abs    Absorption coefficients.
    \param f_mono  Frequency grid.
    \param p_abs   Pressure grid.
    \param t_abs   Temperatures associated with p_abs.
    \param vmrs    Volume mixing ratios of the species.
    \param lines   The spectroscopic line list.
    \param ind_ls  Index to used lineshape function.
    \param ind_lsn Index to used lineshape normalization function.

    \author Stefan Buehler 16.06.2000. */
void abs_species( MATRIX&                  abs,
		  const VECTOR&  	   f_mono,
		  const VECTOR&  	   p_abs,
		  const VECTOR&  	   t_abs,           
		  const VECTOR&            vmr,
		  const ARRAYofLineRecord& lines,
		  const size_t             ind_ls,
		  const size_t             ind_lsn)
{
  // Make lineshape and species lookup data visible:
  extern const ARRAY<LineshapeRecord> lineshape_data;
  extern const ARRAY<LineshapeNormRecord> lineshape_norm_data;

  // speed of light constant
  extern const Numeric SPEED_OF_LIGHT;

  // Boltzmann constant
  extern const Numeric BOLTZMAN_CONST;

  // Avogadros constant
  extern const Numeric AVOGADROS_NUMB;

  // Planck constant
  extern const Numeric PLANCK_CONST;

  // sqrt(ln(2))
  // extern const Numeric SQRT_NAT_LOG_2;

  // PI
  extern const Numeric PI;

  // constant sqrt(1/pi)
  const Numeric sqrt_invPI =  sqrt(1/PI);

  // Constant within the Doppler Broadening calculation:
  const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST *
				      AVOGADROS_NUMB) / SPEED_OF_LIGHT; 

  // dimension of f_mono
  size_t nf = f_mono.dim();

  // Define the vector for the line shape function and the
  // normalization factor of the lineshape here, so that we don't need
  // so many free store allocations.
  VECTOR ls(nf);
  VECTOR fac(nf);

  // Check that p_abs, t_abs, and vmr all have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error. 

  if ( t_abs.size() != p_abs.size() )
    {
      ostringstream os;
      os << "Variable t_abs must have the same dimension as p_abs.\n"
	 << "t_abs.size() = " << t_abs.size() << '\n'
	 << "p_abs.size() = " << p_abs.size();
      throw runtime_error(os.str());
    }

  if ( vmr.size() != p_abs.size() )
    {
      ostringstream os;
      os << "Variable vmr must have the same dimension as p_abs.\n"
	 << "vmr.size() = " << vmr.size() << '\n'
	 << "p_abs.size() = " << p_abs.size();
      throw runtime_error(os.str());
    }

  // Check that the dimension of abs is indeed [f_mono.dim(),
  // p_abs.dim()]:
  if ( abs.dim(1) != nf || abs.dim(2) != p_abs.dim() )
    {
      ostringstream os;
      os << "Variable abs must have dimensions [f_mono.size(),p_abs.size()].\n"
	 << "[abs.dim(1),abs.dim(2)] = [" << abs.dim(1)
	 << ", " << abs.dim(2) << "]\n"
	 << "f_mono.dim() = " << nf << '\n'
	 << "p_abs.dim() = " << p_abs.dim();
      throw runtime_error(os.str());
    }

  // Loop all pressures:
  for ( size_t i=1; i<=p_abs.size(); ++i )
  {

    // store variables p_abs(i) and t_abs(i),
    // this is slightly faster
    Numeric p_i = p_abs(i);
    Numeric t_i = t_abs(i);

    
    //out3 << "  p = " << p_i << " Pa\n";

    // Calculate total number density from pressure and temperature. n
    // = n0*T0/p0 * p/T, ideal gas law
    Numeric n;
    {
      // FIXME: Should these be moved to constants.cc?
      // The number density cab ne calculated as n  = p/KB/t. No new 
      // constants are needed (PE 001215). 
      const Numeric T_0_C = 273.15;  	       /* temp. of 0 Celsius in [K]  */
      const Numeric p_0   = 101300.25; 	       /* standard p in [Pa]        */
      const Numeric n_0   = 2.686763E25;         /* Loschmidt constant [m^-3] */
      const Numeric fac   = n_0 * T_0_C / p_0;
      n = fac * p_i / t_i;
    }

    // For the pressure broadening, we also need the partial pressure:
    const Numeric p_partial = p_i * vmr(i);


    // Loop all lines:
    for ( size_t l=0; l<lines.size(); ++l )
      {

	// lines[l] is used several times, this construct should be
	// faster (Oliver Lemke)
	LineRecord l_l = lines[l];

	// Intensity is already in the right units (Hz*m^2). It also
	// includes already the isotopic ratio. Needs
	// only to be multiplied by total number density, VMR, and lineshape. 
	Numeric intensity = l_l.I0();

	// Lower state energy is in wavenumbers
	Numeric e_lower = l_l.Elow() *  PLANCK_CONST * SPEED_OF_LIGHT
	  * 1E2;

	// Upper state energy
	Numeric e_upper = e_lower + l_l.F() * PLANCK_CONST;

	// get the ratio of the partition function
	Numeric part_fct_ratio = 0.0;
	try 
	  {
	    part_fct_ratio =
	      l_l.IsotopeData().CalculatePartitionFctRatio( t_i );
	  }

	catch (runtime_error e)
	  {
	    ostringstream os;
	    os << "Partition function of "
	       << "Species-Isotope = " << l_l.Name()
	       << "is unknown." << endl;
	  }

	// Boltzmann factors
	Numeric nom = exp(- e_lower / (BOLTZMAN_CONST * t_i ) ) - 
       	              exp(- e_upper /( BOLTZMAN_CONST * t_i ) );

	Numeric denom = exp(- e_lower / (BOLTZMAN_CONST * l_l.Ti0() ) ) - 
       	                exp(- e_upper /( BOLTZMAN_CONST * l_l.Ti0() ) );


	// intensity at temperature
	intensity *= part_fct_ratio * nom / denom;


	// 2. Get pressure broadened line width:
	// (Agam is in Hz/Pa, p_abs is in Pa, gamma is in Hz)
	const Numeric theta = l_l.Tgam() / t_i;
	const Numeric theta_Nair = pow(theta, l_l.Nair());

	Numeric gamma
	  = l_l.Agam() * theta_Nair  * (p_i - p_partial)
	  + l_l.Sgam() * pow(theta, l_l.Nself()) * p_partial;

	// 3. Doppler broadening without the sqrt(ln(2)) factor, which
	// seems to be redundant FIXME: verify .
	Numeric sigma = l_l.F() * doppler_const * 
	  sqrt( t_i / l_l.IsotopeData().Mass());


	// Calculate the line shape:
	lineshape_data[ind_ls].Function()(ls,
					  l_l.F(),
					  gamma,
					  sigma,
					  f_mono,
					  nf);


	// Calculate the chosen normalization factor:
 	lineshape_norm_data[ind_lsn].Function()(fac,
						l_l.F(),
						f_mono,
						nf);



	// Add line to abs:

	// little speeding up factor
	Numeric sfac = 1.0 / sigma * sqrt_invPI;
	for ( size_t j=1; j<=nf; ++j )
	  {
	    abs(j,i) = abs(j,i)
	      + n * vmr(i) * intensity * ls(j) * sfac * fac(j);

	    // if (i == 20)
	    //  {
	    //    cout << f_mono(j) << ' ' << n * intensity << ' ' <<
	    //      abs(j,i) << ' ' <<  ls(j) * factor * f_mono(j) *
	    //      f_mono(j)  << endl;
	    //  }


	  }
      }
  }
}
