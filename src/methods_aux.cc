/* Copyright (C) 2000, 2001, 2002
   Stefan Buehler <sbuehler@uni-bremen.de>

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
  \file   methods_aux.cc
  \brief  Auxiliary material for the workspace
          methods, which used to be in methods.cc. 

  The reason for the separation is that the stuff here hardly ever
  should be changed, whereas methods.cc has to be edited each time a
  new method is added. See methods.h for more documentation.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "arts.h"
#include <map>
#include "make_array.h"
#include "auto_wsv.h"
#include "methods.h"
#include "wsv_aux.h"


//! Initializing constructor for MdRecord
/*!
  This is the only non-trivial constructor, which sets all the
  fields. The flag for supergenericity is not set directly, but
  inferred from the presence of Any_ in the generic input or output
  list. 
*/
MdRecord::MdRecord(const char 		        name[],
		   const char 		        description[],
		   const MakeArray<Index>&      output,
		   const MakeArray<Index>&      input,   
		   const MakeArray<Index>&      goutput,
		   const MakeArray<Index>&      ginput,   
		   const MakeArray<String>&     keywords,
		   const MakeArray<TokValType>& types,
		   bool                         agenda_method) :
    mname(          name            	  ),
    mdescription(   description     	  ),    
    moutput(        output       	  ),  
    minput(         input        	  ),   
    mgoutput(       goutput      	  ),  
    mginput(        ginput       	  ),   
    mkeywords(      keywords     	  ),
    mtypes(         types        	  ),
    magenda_method( agenda_method         ),
    msupergeneric(  false                 )
    { 
      // Initializing the various arrays with input data should now
      // work correctly.  

      // Keywords and type must have the same number of
      // elements. (Types specifies the types associated with each
      // keyword.)
      assert( mkeywords.nelem() == mtypes.nelem() );

      // Find out if this method is supergeneric, and set the flag if
      // yes:
      for ( Index j=0; j<mgoutput.nelem(); ++j )
	if ( Any_ == mgoutput[j] )
	  msupergeneric = true;
      for ( Index j=0; j<mginput.nelem(); ++j )
	if ( Any_ == mginput[j] )
	  msupergeneric = true;
    }


// //! Find out if the method is supergeneric.
// /*! 
//    We will need special treatment for those.

//    \param mdd The record for the method to check.

//    \return True if this is a supergeneric method.
// */
// bool is_supergeneric(const MdRecord& mdd)
// {
//   bool is_supergeneric = false;

//   const ArrayOfIndex&  vgo=mdd.GOutput();   // Generic Output 
//   const ArrayOfIndex&  vgi=mdd.GInput();    // Generic Input

//   for ( Index j=0; j<vgo.nelem(); ++j )
//     if ( Any_ == vgo[j] )		is_supergeneric = true;
//   for ( Index j=0; j<vgi.nelem(); ++j )
//     if ( Any_ == vgi[j] )		is_supergeneric = true;

//   return is_supergeneric;
// }

//! Expand supergeneric record for given group.
/*! 
  This function will substitute any occurance of Any_ in the goutput
  and ginput list of MdRecord by group g. 

  It also adds the group to the name like this: Copy becomes
  Copy_sg_Vector, Copy_sg_Matrix, etc..

  \param g The group for which to expand.
*/
void MdRecord::subst_any_with_group( Index g )
{
  // The group names, we need them for the expansion:
  extern const ArrayOfString wsv_group_names;

  // Make sure they are initialized:
  assert( 0!=wsv_group_names.nelem() );

  // Make sure that this really is a supergeneric method:
  assert( Supergeneric() );

  // Modify the name:
  {
    ostringstream os;
    os << mname << "_sg_" << wsv_group_names[g];
    mname = os.str();
  }
  
  for ( Index j=0; j<mgoutput.nelem(); ++j )
    if ( Any_ == mgoutput[j] )		mgoutput[j] = g;
  for ( Index j=0; j<mginput.nelem(); ++j )
    if ( Any_ == mginput[j] )		mginput[j] = g;

  // After supstitution, the record is no longer supergeneric, it is
  // just a normal method data record:
  msupergeneric = false;
}

//! Expand supergeneric methods.
/*!
  This creates md_data from md_data_raw, by explicitly expanding
  supergeneric methods for all groups. That means, e.g., instead of
  supergeneric method Copy(Any,Any) there will be
  Copy_sg_Vector(Vector,Vector), Copy_sg_Matrix(Matrix,Matrix), etc..

  Not only the goutput and ginput lists are manipulated, also the
  method name.
*/
void expand_md_data_raw_to_md_data()
{
  extern const Array<MdRecord> md_data_raw;
  extern Array<MdRecord>       md_data;

  // The group names, we need them for the expansion:
  extern const ArrayOfString wsv_group_names;

  // Make sure that they have been initialized:
  assert ( 0!=wsv_group_names.nelem() );

  // Reset md_data, just in case:
  md_data.resize(0);

  for ( Index i=0; i<md_data_raw.nelem(); ++i )
    {
      const MdRecord& mdd = md_data_raw[i];

      if ( !mdd.Supergeneric() )
	{
	  md_data.push_back(mdd);
	}
      else
	{
	  // Special treatment for supergeneric methods:

	  // We have to create method records for all groups!
	  for ( Index j=0; j<wsv_group_names.nelem(); ++j )
	    {
	      // Any_ itself is also a group, but we don't want to
	      // create a record for Any_!
	      if ( Any_ != j )
		{
		  // Not a reference but a copy this time, since we
		  // have to manipulate this.
		  MdRecord mdlocal = mdd;
		      
		  mdlocal.subst_any_with_group( j );

		  md_data.push_back(mdlocal);
		}
	    }

	}
    }
}

//! Define MdMap. 
/*!
  MdMap can be used to find method data by method name.
*/
void define_md_map()
{
  extern const Array<MdRecord> md_data;
  extern std::map<String, Index> MdMap;

  for ( Index i=0 ; i<md_data.nelem() ; ++i)
    {
      MdMap[md_data[i].Name()] = i;
    }
}

//! Define MdRawMap. 
/*!
  MdRawMap can be used to find method data by method name. In the
  md_data_raw lookup table. This is the method table before expansion
  of supergeneric methods.
*/
void define_md_raw_map()
{
  extern const Array<MdRecord> md_data_raw;
  extern std::map<String, Index> MdRawMap;

  for ( Index i=0 ; i<md_data_raw.nelem() ; ++i)
    {
      MdRawMap[md_data_raw[i].Name()] = i;
    }
}


ostream& MdRecord::PrintTemplate(ostream& os,
				 bool show_description) const
{
  extern const  ArrayOfString wsv_group_names;

  if (show_description)
    {
      // FIXME: Print description String!
    }
  
  os << Name();

  // Is this a generic method? -- Then we need round braces.
  if ( 0 != GOutput().nelem()+GInput().nelem() )
    {
      // First entry needs to comma before:
      bool first=true;

      os << '(';

      for (Index i=0; i<GOutput().nelem(); ++i)
	{
	  if (first) first=false;
	  else os << ",\n";

	  os << wsv_group_names[GOutput()[i]];
	}

      for (Index i=0; i<GInput().nelem(); ++i)
	{
	  if (first) first=false;
	  else os << ",\n";

	  os << wsv_group_names[GInput()[i]];
	}

      os << ')';
    }

  // Now the keywords:

  os << '{';

  // Determine the length of the longest keyword:
  Index maxsize = 0;
  for (Index i=0; i<Keywords().nelem(); ++i)
    if ( Keywords()[i].nelem() > maxsize )
      maxsize = Keywords()[i].nelem();


  for (Index i=0; i<Keywords().nelem(); ++i)
    {
      os << "\t" << setw(maxsize)
	 << Keywords()[i] << " = \n";
    }

  os << '}';

  return os;
}

//! Output operator for MdRecord.
ostream& operator<<(ostream& os, const MdRecord& mdr)
{
  extern const Array<WsvRecord> wsv_data;
  extern const ArrayOfString wsv_group_names;
  extern const String TokValTypeName[];
  bool first;

  os << "\n*-------------------------------------------------------------------*\n"
     << "Workspace method = " << mdr.Name() << 
        "\n---------------------------------------------------------------------\n"
     << "\n" << mdr.Description() << "\n" << 
        "\n---------------------------------------------------------------------\n";

  //  os << "\n-----\nName = " << mdr.Name() << '\n\n'
  //     << "Description =\n" << mdr.Description() << "\n\n";

  // Output:
  first = true;
  os << "Output = ";
  for ( Index i=0; i<mdr.Output().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Output()[i]].Name();
    }
  os << '\n';

  // Input:
  first = true;
  os << "Input = ";
  for ( Index i=0; i<mdr.Input().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Input()[i]].Name();
    }
  os << '\n';
      
  // GOutput:
  first = true;
  os << "GOutput = ";
  for ( Index i=0; i<mdr.GOutput().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GOutput()[i]];
    }
  os << '\n';

  // GInput:
  first = true;
  os << "GInput = ";
  for ( Index i=0; i<mdr.GInput().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GInput()[i]];
    }
  os << '\n';

  // Keywords:
  first = true;
  os << "Keywords = ";
  for ( Index i=0; i<mdr.Keywords().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << mdr.Keywords()[i];
    }
  os << '\n';

  // Types:
  first = true;
  os << "Types = ";
  for ( Index i=0; i<mdr.Types().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << TokValTypeName[mdr.Types()[i]];
    }
  os << "\n*-------------------------------------------------------------------*\n";

  return os;
}

