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


void limit_line_length( ostream& os, ostringstream& curline,
                        ostringstream& token, const String& indent,
                        size_t linelen);

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
                   const MakeArray<String>&     gout,
                   const MakeArray<String>&     gouttype,
                   const MakeArray<String>&     goutdesc,
                   const MakeArray<String>&     input,   
                   const MakeArray<String>&     gin,   
                   const MakeArray<String>&     gintype,   
                   const MakeArray<String>&     gindefault,
                   const MakeArray<String>&     gindesc,
                   bool                         set_method,
                   bool                         agenda_method,
                   bool                         suppress_header,
                   bool                         pass_workspace,
                   bool                         pass_wsv_names ) :
    mname(            name                ),
    mdescription(     description         ),    
    mauthors(         authors             ),
    moutput(          0                   ),  
    mgout(            gout                ),  
    mgouttype(        0                   ),  
    mgoutdesc(        goutdesc            ),
    minput(           0                   ),   
    mgin(             gin                 ),   
    mgintype(         0                   ),   
    mgindefault(      gindefault          ),
    mgindesc(         gindesc             ),
    mset_method(      set_method          ),
    magenda_method(   agenda_method       ),
    msupergeneric(    false               ),
    msuppress_header( suppress_header     ),
    mpass_workspace(  pass_workspace      ),
    mpass_wsv_names(  pass_wsv_names      ),
    mactual_group(    -1                  )
    { 
      // Initializing the various arrays with input data should now
      // work correctly.  

      // Generic variable names, types and defaults must have the same number of
      // elements. (Defaults specifies the default values associated with each
      // generic input.)
      assert( mgout.nelem() == gouttype.nelem() );
      assert( mgout.nelem() == goutdesc.nelem() );
      assert( mgin.nelem() == mgindefault.nelem() );
      assert( mgin.nelem() == gintype.nelem() );
      assert( mgin.nelem() == gindesc.nelem() );

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
      mgouttype.resize(gouttype.nelem());
      for ( Index j=0; j<gouttype.nelem(); ++j )
        {
          mgouttype[j] = get_wsv_group_id (gouttype[j]);
          if (mgouttype[j] == -1)
            {
              ostringstream os;
              os << "Unknown WSV Group " << gouttype[j] << " for generic output "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      mgintype.resize(gintype.nelem());
      for ( Index j=0; j<gintype.nelem(); ++j )
        {
          mgintype[j] = get_wsv_group_id (gintype[j]);
          if (mgintype[j] == -1)
            {
              ostringstream os;
              os << "Unknown WSV Group " << gintype[j] << " for generic input "
                << "in WSM " << mname;
              throw runtime_error( os.str() );
            }
        }

      // Find out if this method is supergeneric, and set the flag if
      // yes:
      for ( Index j=0; j<mgouttype.nelem(); ++j )
        if ( get_wsv_group_id("Any") == mgouttype[j] )
          msupergeneric = true;
      for ( Index j=0; j<mgintype.nelem(); ++j )
        if ( get_wsv_group_id("Any") == mgintype[j] )
          msupergeneric = true;
    }


//! Expand supergeneric record for given group.
/*! 
  This function will substitute any occurance of Any_ in the GOutType
  and GInType list of MdRecord by group g. 

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
  
  for ( Index j=0; j<mgouttype.nelem(); ++j )
    if ( wsv_group_id_Any == mgouttype[j] )          mgouttype[j] = g;
  for ( Index j=0; j<mgintype.nelem(); ++j )
    if ( wsv_group_id_Any == mgintype[j] )           mgintype[j] = g;

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

  Not only the GOutType and GInType lists are manipulated, also the
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


void format_paragraph (String &s, const String &indent, const size_t linelen)
{
  String out = indent;
  String token;
  size_t lastreturn = 0;
  for (size_t i = 0; i < s.length(); i++)
    {
      if (s[i] == '\n') s[i] = ' ';
      token += s[i];
      if (s[i] == ' ')
        {
          if (out.length() - lastreturn + token.length() > linelen)
            {
              out += '\n' + indent;
              lastreturn = i;
            }
          out += token;
          token = "";
        }
    }

  if (token.length())
    {
      if (out.length() - lastreturn + token.length() > linelen)
        out += '\n' + indent;
      out += token;
    }
  s = out;
}


void get_short_wsv_description(String &s, const String &desc)
{
  Index pos;
  Index pos2;

  // Find the end of the short WSV description
  pos = desc.find(".\n");
  pos2 = desc.find(". ");
  if (pos == String::npos || pos2 != String::npos && pos2 < pos) pos = pos2;
  if (pos == String::npos) pos = desc.find("\n");
  if (pos != String::npos)
    s = desc.substr(0, pos+1);
  else
    s = desc;

  // Replace any newlines inside the description with spaces
  while ((pos = s.find("\n")) != String::npos)
    {
      s.replace (pos, 1, " ");
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
  if ( 0 != GOutType().nelem()+GInType().nelem() )
    {
      // First entry needs to comma before:
      bool first=true;

      os << '(';

      for (Index i=0; i<GOutType().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",\n";

          os << wsv_group_names[GOutType()[i]];
        }

      for (Index i=0; i<GInType().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",\n";

          os << wsv_group_names[GInType()[i]];
        }

      os << ')';
    }

  // Now the keywords:

  os << '{';

  // Determine the length of the longest keyword:
  Index maxsize = 0;
  for (Index i=0; i<GIn().nelem(); ++i)
    if ( GIn()[i].nelem() > maxsize )
      maxsize = GIn()[i].nelem();


  for (Index i=0; i<GIn().nelem(); ++i)
    {
      os << "\t" << setw(maxsize)
         << GIn()[i] << " = \n";
    }

  os << '}';

  return os;
}


//! Limit length of output
/*! Automatically inserts linebreaks at certain length.
  
  \author Oliver Lemke
  \date   2008-09-03
*/
void limit_line_length( ostream& os,
                        ostringstream& curline,
                        ostringstream& token,
                        const String& indent,
                        size_t linelen)
{
  if (indent.length() + curline.str().length() + token.str().length() > linelen)
    {
      os << curline.str() << endl << indent;
      curline.str("");
    }
  curline << token.str();
  token.str("");
}


//! Output operator for MdRecord.
ostream& operator<<(ostream& os, const MdRecord& mdr)
{
  extern const ArrayOfString wsv_group_names;
  bool first;
  ostringstream buf;
  ostringstream param;
  String indent = "";
  const size_t linelen = 68;

  os << "\n*-------------------------------------------------------------------*\n"
     << "Workspace method = " << mdr.Name() << 
        "\n---------------------------------------------------------------------\n"
     << "\n" << mdr.Description() << "\n";

  if (mdr.Description()[mdr.Description().nelem() - 1] != '\n')
    {
      os << "\n";
    }

  // Print the method's synopsis
  while (indent.length() < mdr.Name().length() + 2) indent += ' ';

  os << "\nSynopsis:\n\n";
  buf << mdr.Name() << "( ";
  first = true;
  for ( Index i=0; i<mdr.Out().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      param << Workspace::wsv_data[mdr.Out()[i]].Name();

      limit_line_length( os, buf, param, indent, linelen );
    }

  for ( Index i=0; i<mdr.GOutType().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      if (mdr.GOut()[i].length())
        param << mdr.GOut()[i];
      else
        param << "gout" << i;

      limit_line_length( os, buf, param, indent, linelen );
    }

  ArrayOfIndex inonly;
  mdr.input_only (inonly);
  for ( Index i=0; i<inonly.nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      param << Workspace::wsv_data[inonly[i]].Name();

      limit_line_length( os, buf, param, indent, linelen );
    }

  for ( Index i=0; i<mdr.GInType().nelem(); ++i )
    {
      if (first) first=false; else buf << ", ";
      if (mdr.GIn()[i].length())
        {
          param << mdr.GIn()[i];
        }
      else
        {
          param << "gin" << i;
        }

      limit_line_length( os, buf, param, indent, linelen );
    }
  if (buf.str().length()) os << buf.str();

  os << " )\n\n\n";

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

  os << "\n\nVariables:\n\n";

  //  os << "\n-----\nName = " << mdr.Name() << '\n\n'
  //     << "Description =\n" << mdr.Description() << "\n\n";

  // Out:
  first = true;
  indent = "      ";
  String desc;
  for ( Index i=0; i<mdr.Out().nelem(); ++i )
    {
      buf.str("");
      buf <<    "OUT   ";

      buf << Workspace::wsv_data[mdr.Out()[i]].Name();
      buf << " (";
      buf << wsv_group_names[Workspace::wsv_data[mdr.Out()[i]].Group()];
      buf << "): ";

      get_short_wsv_description(
        desc, Workspace::wsv_data[mdr.Out()[i]].Description());

      if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen);
          buf << endl << desc;
        }
      else
        {
          buf << desc;
        }

      os << buf.str() << endl;
    }

  for ( Index i=0; i<mdr.In().nelem(); ++i )
    {
      buf.str("");
      buf <<    "IN    ";

      buf << Workspace::wsv_data[mdr.In()[i]].Name();
      buf << " (";
      buf << wsv_group_names[Workspace::wsv_data[mdr.In()[i]].Group()];
      buf << "): ";

      get_short_wsv_description(
        desc, Workspace::wsv_data[mdr.In()[i]].Description());

      if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen);
          buf << endl << desc;
        }
      else
        {
          buf << desc;
        }

      os << buf.str() << endl;
    }

  for ( Index i=0; i<mdr.GOut().nelem(); ++i )
    {
      buf.str("");
      buf <<    "GOUT  ";

      if (mdr.GOut()[i].length())
        buf << mdr.GOut()[i];
      else
        buf << "gout" << i;
      buf << " (";
      buf << wsv_group_names[mdr.GOutType()[i]];

      buf << "): ";
      desc = mdr.GOutDescription()[i];
      if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen);
          buf << endl << desc;
        }
      else
        {
          buf << desc;
        }

      os << buf.str() << endl;
    }

  for ( Index i=0; i<mdr.GIn().nelem(); ++i )
    {
      buf.str("");
      buf <<    "GIN   ";

      if (mdr.GIn()[i].length())
        buf << mdr.GIn()[i];
      else
        buf << "gin" << i;
      buf << " (";
      buf << wsv_group_names[mdr.GInType()[i]];

      if (mdr.GInDefault()[i] != NODEF)
        {
          buf << ", Default: ";
          if (mdr.GInType()[i] == get_wsv_group_id ("String"))
            {
              buf << "\"" << mdr.GInDefault()[i] << "\"";
            }
          else
            {
              buf << mdr.GInDefault()[i];
            }

        }

      buf << "): ";
      desc = mdr.GInDescription()[i];
      if (buf.str().length() + desc.length() > linelen)
        {
          format_paragraph (desc, indent, linelen);
          buf << endl << desc;
        }
      else
        {
          buf << desc;
        }

      os << buf.str() << endl;
    }

  os << "\n*-------------------------------------------------------------------*\n";

  return os;
}

