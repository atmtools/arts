/* Copyright (C) 2000-2008 Stefan Buehler <sbuehler@ltu.se>

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
#include "methods.h"
#include "wsv_aux.h"
#include "workspace_ng.h"

//! Initializing constructor for MdRecord
/*!
  This is the only non-trivial constructor, which sets all the
  fields. The flag for supergenericity is not set directly, but
  inferred from the presence of Any_ in the generic input or output
  list. 
*/
MdRecord::MdRecord(const char                   name[],
                   const char                   description[],
                   const MakeArray<String>&     authors,
                   const MakeArray<String>&     output,
                   const MakeArray<String>&     input,   
                   const MakeArray<String>&     goutput,
                   const MakeArray<String>&     ginput,   
                   const MakeArray<String>&     keywords,
                   const MakeArray<String>&     defaults,
                   const MakeArray<String>&     types,
                   bool                         agenda_method,
                   bool                         suppress_header,
                   bool                         pass_workspace,
                   bool                         pass_wsv_names ) :
    mname(          name                  ),
    mdescription(   description           ),    
    mauthors(       authors               ),
    moutput(        0                     ),  
    minput(         0                     ),   
    mgoutput(       0                     ),  
    mginput(        0                     ),   
    mkeywords(      keywords              ),
    mdefaults(      defaults              ),
    mtypes(         0                     ),
    magenda_method(   agenda_method       ),
    msupergeneric(    false               ),
    msuppress_header( suppress_header     ),
    mpass_workspace( pass_workspace       ),
    mpass_wsv_names( pass_wsv_names       ),
    mactual_group( -1 )
    { 
      // Initializing the various arrays with input data should now
      // work correctly.  

      // Keywords and type must have the same number of
      // elements. (Types specifies the types associated with each
      // keyword.)
      mtypes.resize(types.nelem());
      assert( mkeywords.nelem() == mtypes.nelem() );

      // Keywords and Defaults must have the same number of
      // elements. (Defaults specifies the default values associated with each
      // keyword.)
      assert( mkeywords.nelem() == mdefaults.nelem() );

      // Map the WSV names to indexes
      moutput.resize(output.nelem());
      for ( Index j=0; j<output.nelem(); ++j )
        {
          moutput[j] = get_wsv_id (output[j]);
          if (moutput[j] == -1)
            {
              ostringstream os;

              os << "Unknown WSV " << output[j] << " for output "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      minput.resize(input.nelem());
      for ( Index j=0; j<input.nelem(); ++j )
        {
          minput[j] = get_wsv_id (input[j]);
          if (minput[j] == -1)
            {
              ostringstream os;

              os << "Unknown WSV " << input[j] << " for input "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      // Map the group names to groups' indexes
      mgoutput.resize(goutput.nelem());
      for ( Index j=0; j<goutput.nelem(); ++j )
        {
          mgoutput[j] = get_wsv_group_id (goutput[j]);
          if (mgoutput[j] == -1)
            {
              ostringstream os;

              os << "Unknown WSV Group " << goutput[j] << " for generic output "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      mginput.resize(ginput.nelem());
      for ( Index j=0; j<ginput.nelem(); ++j )
        {
          mginput[j] = get_wsv_group_id (ginput[j]);
          if (mginput[j] == -1)
            {
              ostringstream os;

              os << "Unknown WSV Group " << ginput[j] << " for generic input "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }


      // Map the keyword types to group indexes
      for ( Index j=0; j<types.nelem(); ++j )
        {
          mtypes[j] = get_wsv_group_id (types[j]);
          if (mtypes[j] == -1)
            {
              ostringstream os;

              os << "Unknown Group " << types[j] << " for keyword"
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      // Find out if this method is supergeneric, and set the flag if
      // yes:
      for ( Index j=0; j<mgoutput.nelem(); ++j )
        if ( get_wsv_group_id("Any") == mgoutput[j] )
          msupergeneric = true;
      for ( Index j=0; j<mginput.nelem(); ++j )
        if ( get_wsv_group_id("Any") == mginput[j] )
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
//     if ( Any_ == vgo[j] )            is_supergeneric = true;
//   for ( Index j=0; j<vgi.nelem(); ++j )
//     if ( Any_ == vgi[j] )            is_supergeneric = true;

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
  const Index wsv_group_id_Any = get_wsv_group_id("Any");
  // The group names, we need them for the expansion:
  DEBUG_ONLY (extern const ArrayOfString wsv_group_names);

  // Make sure they are initialized:
  assert( 0 != wsv_group_names.nelem() );

  // Make sure that g is in the allowed range, which means
  // 0<=g<wsv_group_names.nelem() and g != Any_
  assert( 0 <= g );
  assert( wsv_group_id_Any != g );
  assert( g < wsv_group_names.nelem() );

  // Make sure that this really is a supergeneric method:
  assert( Supergeneric() );

  // Modify the name:
//   {
//     ostringstream os;
//     os << mname << "_sg_" << wsv_group_names[g];
//     mname = os.str();
//   }
  
  for ( Index j=0; j<mgoutput.nelem(); ++j )
    if ( wsv_group_id_Any == mgoutput[j] )          mgoutput[j] = g;
  for ( Index j=0; j<mginput.nelem(); ++j )
    if ( wsv_group_id_Any == mginput[j] )           mginput[j] = g;

  // Set the field for the actual group:
  mactual_group = g;
}


//! Get list of input only WSVs.
/*!
  This function returns an array with the indexes of WSVs which are
  only input variables but not output.

  \param[out] inonly Index array of input only WSVs.

  \author Oliver Lemke
  \date   2008-02-27
*/
void MdRecord::input_only(ArrayOfIndex& inonly) const
{
  inonly = minput;    // Input
  for (ArrayOfIndex::const_iterator j=moutput.begin(); j<moutput.end(); ++j)
    for (ArrayOfIndex::iterator k=inonly.begin(); k<inonly.end(); ++k)
      if ( *j == *k )
        {
          //              erase_vector_element(vi,k);
          k = inonly.erase(k) - 1;
          // We need the -1 here, otherwise due to the
          // following increment we would miss the element
          // behind the erased one, which is now at the
          // position of the erased one.
        }
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

  const Index wsv_group_id_Any = get_wsv_group_id("Any");

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
              if ( wsv_group_id_Any != j )
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
  // md_data is constant here and should never be changed
  extern Array<MdRecord> md_data;
  extern map<String, Index> MdMap;
  extern const ArrayOfString wsv_group_names;

  // Check that md_data and wsv_group_names have already be defined:
  assert( 0 != md_data.nelem() );
  assert( 0 != wsv_group_names.nelem() );

  for ( Index i=0 ; i<md_data.nelem() ; ++i)
    {
      const MdRecord& mdd = md_data[i];

//       cout << "mdd.ActualGroup() = "
//         << mdd.ActualGroup() << "\n";

//       cout << "wsv_group_names[mdd.ActualGroup()] = "
//         << wsv_group_names[mdd.ActualGroup()] << "\n";

      // For supergeneric methods, add group to method name
      String methodname;
      ostringstream os;
      if ( mdd.Supergeneric() )
        {
          os << mdd.Name() << "_sg_"
             << wsv_group_names[mdd.ActualGroup()];
        }
      else
        {
          os << mdd.Name();
        }
      methodname = os.str();

      MdMap[methodname] = i;
    }
}

//! Define MdRawMap. 
/*!
  MdRawMap can be used to find method data by method name. In the
  md_data_raw lookup table. This is the method table before expansion
  of supergeneric methods.

  We add the _sg_Type string to the methodname here, so that
  supergeneric methods can be picked out for the right type.
*/
void define_md_raw_map()
{
  extern const Array<MdRecord> md_data_raw;
  extern map<String, Index> MdRawMap;

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
  extern const ArrayOfString wsv_group_names;
  extern const String TokValTypeName[];
  bool first;

  os << "\n*-------------------------------------------------------------------*\n"
     << "Workspace method = " << mdr.Name() << 
        "\n---------------------------------------------------------------------\n"
     << "\n" << mdr.Description() << "\n";

  if (mdr.Description()[mdr.Description().nelem() - 1] != '\n')
    {
      os << "\n";
    }

  // Print the method's synopsis
  {
    bool separate_lines = false;
    String indent = "";

    // If the method has more than 4 arguments, put them on
    // separate lines for better readability
    if (mdr.Output().nelem() + mdr.GOutput().nelem() + mdr.Input().nelem()
        + mdr.GInput().nelem() + mdr.Keywords().nelem() > 4)
      {
        separate_lines = true;
        indent = "\n";
        for (size_t i = 0; i < mdr.Name().length() + 2; i++)
          {
            indent += ' ';
          }
      }

    os << "Synopsis (Arts2 Syntax):\n\n";
    os << mdr.Name() << "( ";
    first = true;
    for ( Index i=0; i<mdr.Output().nelem(); ++i )
      {
        if (first) first=false; else os << ", " << indent;
        os << Workspace::wsv_data[mdr.Output()[i]].Name();
      }

    for ( Index i=0; i<mdr.GOutput().nelem(); ++i )
      {
        if (first) first=false; else os << ", " << indent;
        os << wsv_group_names[mdr.GOutput()[i]];
      }

    ArrayOfIndex inonly;
    mdr.input_only (inonly);
    for ( Index i=0; i<inonly.nelem(); ++i )
      {
        if (first) first=false; else os << ", " << indent;
        os << Workspace::wsv_data[inonly[i]].Name();
      }

    for ( Index i=0; i<mdr.GInput().nelem(); ++i )
      {
        if (first) first=false; else os << ", " << indent;
        os << wsv_group_names[mdr.GInput()[i]];
      }

    for ( Index i=0; i<mdr.Keywords().nelem(); ++i )
      {
        if (first) first=false; else os << ", " << indent;
        os << '"' << mdr.Keywords()[i]<< '"';
      }
    os << " )\n\n";
  }

  {
    bool is_first_author = true;
    for (Index i = 0; i < mdr.Authors().nelem(); i++)
      {
        if (is_first_author)
          {
            os << "Authors: ";
            is_first_author = false;
          }
        else
          os << ", ";

        os << mdr.Authors()[i];
      }
    os << "\n";
  }

  os << "\n---------------------------------------------------------------------\n";

  //  os << "\n-----\nName = " << mdr.Name() << '\n\n'
  //     << "Description =\n" << mdr.Description() << "\n\n";

  // Output:
  first = true;
  os << "Output = ";
  for ( Index i=0; i<mdr.Output().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << Workspace::wsv_data[mdr.Output()[i]].Name();
    }
  os << '\n';

  // Input:
  first = true;
  os << "Input = ";
  for ( Index i=0; i<mdr.Input().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << Workspace::wsv_data[mdr.Input()[i]].Name();
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

  // Defaults:
  first = true;
  os << "Defaults = ";
  for ( Index i=0; i<mdr.Defaults().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      if (mdr.Defaults()[i] != NODEF)
        os << mdr.Defaults()[i];
      else
        os << "none";
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

