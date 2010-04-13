/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_compound_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.
*/

#include "arts.h"
#include "xml_io_private.h"
#include "xml_io_types.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "gridded_fields.h"
#include "jacobian.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////


//=== GasAbsLookup ===========================================================

//! Reads GasAbsLookup from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gal     GasAbsLookup return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GasAbsLookup& gal,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GasAbsLookup");

  xml_read_from_stream (is_xml, gal.species, pbifs);
  xml_read_from_stream (is_xml, gal.nonlinear_species, pbifs);
  xml_read_from_stream (is_xml, gal.f_grid, pbifs);
  xml_read_from_stream (is_xml, gal.p_grid, pbifs);
  xml_read_from_stream (is_xml, gal.vmrs_ref, pbifs);
  xml_read_from_stream (is_xml, gal.t_ref, pbifs);
  xml_read_from_stream (is_xml, gal.t_pert, pbifs);
  xml_read_from_stream (is_xml, gal.nls_pert, pbifs);
  xml_read_from_stream (is_xml, gal.xsec, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GasAbsLookup");
}


//! Writes GasAbsLookup to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gal     GasAbsLookup
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GasAbsLookup& gal,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GasAbsLookup");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, gal.species, pbofs);
  xml_write_to_stream (os_xml, gal.nonlinear_species, pbofs,
                       "NonlinearSpecies");
  xml_write_to_stream (os_xml, gal.f_grid, pbofs, "FrequencyGrid");
  xml_write_to_stream (os_xml, gal.p_grid, pbofs, "PressureGrid");
  xml_write_to_stream (os_xml, gal.vmrs_ref, pbofs, "ReferenceVmrProfiles");
  xml_write_to_stream (os_xml, gal.t_ref, pbofs, "ReferenceTemperatureProfile");
  xml_write_to_stream (os_xml, gal.t_pert, pbofs, "TemperaturePertubations");
  xml_write_to_stream (os_xml, gal.nls_pert, pbofs,
                       "NonlinearSpeciesVmrPertubations");
  xml_write_to_stream (os_xml, gal.xsec, pbofs, "AbsorptionCrossSections");

  close_tag.set_name ("/GasAbsLookup");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== GriddedField ===========================================================

//! Reads the grids for gridded fields from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GriddedField& gfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  for(Index i=0; i < gfield.get_dim(); i++)
    {
      tag.read_from_stream (is_xml);
      if (tag.get_name() == "Vector")
        {
          String s;
          tag.get_attribute_value ("name", s);
          if (s.length())
            gfield.set_grid_name (i, s);

          Vector v;
          xml_parse_from_stream (is_xml, v, pbifs, tag);
          gfield.set_grid (i, v);
          tag.read_from_stream (is_xml);
          tag.check_name ("/Vector");
        }
      else if (tag.get_name() == "Array")
        {
          String s;
          tag.get_attribute_value ("name", s);
          if (s.length())
            gfield.set_grid_name (i, s);

          tag.get_attribute_value ("type", s);
          if (s == "String")
            {
              ArrayOfString as;
              xml_parse_from_stream (is_xml, as, pbifs, tag);
              gfield.set_grid (i, as);
              tag.read_from_stream (is_xml);
              tag.check_name ("/Array");

            }
          else
            {
              xml_parse_error ("Grids must be of type <Vector> or <ArrayOfString> but <ArrayOf"
                               + s + "> found.");
            }
        }
      else
        {
          xml_parse_error ("Grids must be of type <Vector> or <ArrayOfString> but <"
                           + tag.get_name() + "> found.");
        }
    }
}


//! Writes the grids for gridded fields to an XML input stream
/*!
  \param os_xml  XML output stream
  \param gfield  GriddedField with the grids
  \param pbofs   Pointer to binary output stream. NULL in case of ASCII file.
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GriddedField& gfield,
                     bofstream *pbofs)
{
  for (Index i = 0; i < gfield.get_dim(); i++)
    {
      switch (gfield.get_grid_type(i))
        {
        case GRIDTYPE_NUMERIC:
          xml_write_to_stream (os_xml, gfield.get_numeric_grid(i),
                               pbofs, gfield.get_grid_name(i));
          break;
        case GRIDTYPE_STRING:
          xml_write_to_stream (os_xml, gfield.get_string_grid(i),
                               pbofs, gfield.get_grid_name(i));
          break;
        }
    }
}


//=== GriddedField1 ===========================================================

//! Reads GriddedField1 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField1 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GriddedField1& gfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GriddedField1");

  String s;
  tag.get_attribute_value ("name", s);
  if (s.length())
    gfield.set_name (s);

  xml_read_from_stream (is_xml, (GriddedField&)gfield, pbifs);
  xml_read_from_stream (is_xml, (Vector&)gfield, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GriddedField1");

  if (!gfield.checksize())
    throw runtime_error ("Dimension mismatch between the grids and the data");
}


//! Writes GriddedField1 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField1
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GriddedField1& gfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GriddedField1");
  if (!name.length () && (gfield.get_name().length ()))
    open_tag.add_attribute ("name", gfield.get_name());
  else if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream (os_xml, (GriddedField&)gfield, pbofs);
  xml_write_to_stream (os_xml, (Vector&)gfield, pbofs, "Data");

  close_tag.set_name ("/GriddedField1");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== GriddedField2 ===========================================================

//! Reads GriddedField2 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField2 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GriddedField2& gfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GriddedField3");

  String s;
  tag.get_attribute_value ("name", s);
  if (s.length())
    gfield.set_name (s);

  xml_read_from_stream (is_xml, (GriddedField&)gfield, pbifs);
  xml_read_from_stream (is_xml, (Matrix&)gfield, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GriddedField3");

  if (!gfield.checksize())
    throw runtime_error ("Dimension mismatch between the grids and the data");
}


//! Writes GriddedField2 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField3
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GriddedField2& gfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GriddedField2");
  if (!name.length () && (gfield.get_name().length ()))
    open_tag.add_attribute ("name", gfield.get_name());
  else if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream (os_xml, (GriddedField&)gfield, pbofs);
  xml_write_to_stream (os_xml, (Matrix&)gfield, pbofs, "Data");

  close_tag.set_name ("/GriddedField2");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== GriddedField3 ===========================================================

//! Reads GriddedField3 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField3 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GriddedField3& gfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GriddedField3");

  String s;
  tag.get_attribute_value ("name", s);
  if (s.length())
    gfield.set_name (s);

  xml_read_from_stream (is_xml, (GriddedField&)gfield, pbifs);
  xml_read_from_stream (is_xml, (Tensor3&)gfield, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GriddedField3");

  if (!gfield.checksize())
    throw runtime_error ("Dimension mismatch between the grids and the data");
}


//! Writes GriddedField3 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField3
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GriddedField3& gfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GriddedField3");
  if (!name.length () && (gfield.get_name().length ()))
    open_tag.add_attribute ("name", gfield.get_name());
  else if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream (os_xml, (GriddedField&)gfield, pbofs);
  xml_write_to_stream (os_xml, (Tensor3&)gfield, pbofs, "Data");

  close_tag.set_name ("/GriddedField3");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== GriddedField4 ===========================================================

//! Reads GriddedField4 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField4 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GriddedField4& gfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GriddedField4");

  String s;
  tag.get_attribute_value ("name", s);
  if (s.length())
    gfield.set_name (s);

  xml_read_from_stream (is_xml, (GriddedField&)gfield, pbifs);
  xml_read_from_stream (is_xml, (Tensor4&)gfield, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GriddedField4");
}


//! Writes GriddedField4 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField3
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GriddedField4& gfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GriddedField3");
  if (!name.length () && (gfield.get_name().length ()))
    open_tag.add_attribute ("name", gfield.get_name());
  else if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, (GriddedField&)gfield, pbofs);
  xml_write_to_stream (os_xml, (Tensor4&)gfield, pbofs, "Data");

  close_tag.set_name ("/GriddedField3");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== GridPos =====================================================

//! Reads GridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gpos    GridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      GridPos& gpos,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("GridPos");

  xml_read_from_stream (is_xml, gpos.idx, pbifs);
  xml_read_from_stream (is_xml, gpos.fd[0], pbifs);
  xml_read_from_stream (is_xml, gpos.fd[1], pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/GridPos");
}

//! Writes GridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gpos    GridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const GridPos& gpos,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GridPos");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, gpos.idx, pbofs,
                       "OriginalGridIndexBelowInterpolationPoint");
  xml_write_to_stream (os_xml, gpos.fd[0], pbofs,
                       "FractionalDistanceToNextPoint_1");
  xml_write_to_stream (os_xml, gpos.fd[1], pbofs,
                       "FractionalDistanceToNextPoint_2");

  close_tag.set_name ("/GridPos");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== IsotopeRecord ================================================

//! Reads IsotopeRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param irecord  SpeciesRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      IsotopeRecord& irecord,
                      bifstream *pbifs)
{
  ArtsXMLTag    tag;
  String        name;
  Numeric       abundance;
  Numeric       mass;
  Index         mytrantag;
  Index         hitrantag;
  ArrayOfIndex  jpltags;

  tag.read_from_stream (is_xml);
  tag.check_name ("IsotopeRecord");

  xml_read_from_stream (is_xml, name, pbifs);
  xml_read_from_stream (is_xml, abundance, pbifs);
  xml_read_from_stream (is_xml, mass, pbifs);
  xml_read_from_stream (is_xml, mytrantag, pbifs);
  xml_read_from_stream (is_xml, hitrantag, pbifs);
  xml_read_from_stream (is_xml, jpltags, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/IsotopeRecord");

  irecord = IsotopeRecord (name, abundance, mass, mytrantag, hitrantag,
                           jpltags);
}


//! Writes IsotopeRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param irecord  IsotopeRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const IsotopeRecord& irecord,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("IsotopeRecord");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream (os_xml, irecord.Name (), pbofs, "Name");
  xml_write_to_stream (os_xml, irecord.Abundance (), pbofs, "Abundance");
  xml_write_to_stream (os_xml, irecord.Mass (), pbofs, "Mass");
  xml_write_to_stream (os_xml, irecord.MytranTag (), pbofs, "MytranTag");
  xml_write_to_stream (os_xml, irecord.HitranTag (), pbofs, "HitranTag");
  xml_write_to_stream (os_xml, irecord.JplTags (), pbofs, "JplTags");

  close_tag.set_name ("/IsotopeRecord");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}



//=== Ppath =====================================================

//! Reads Ppath from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ppath   Ppath return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      Ppath& ppath,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is_xml);
  tag.check_name ("Ppath");

  xml_read_from_stream (is_xml, ppath.dim, pbifs);
  xml_read_from_stream (is_xml, ppath.np, pbifs);
  xml_read_from_stream (is_xml, ppath.refraction, pbifs);
  xml_read_from_stream (is_xml, ppath.method, pbifs);
  xml_read_from_stream (is_xml, ppath.constant, pbifs);
  xml_read_from_stream (is_xml, ppath.pos, pbifs);
  xml_read_from_stream (is_xml, ppath.z, pbifs);
  xml_read_from_stream (is_xml, ppath.l_step, pbifs);
  xml_read_from_stream (is_xml, ppath.gp_p, pbifs);
  xml_read_from_stream (is_xml, ppath.gp_lat, pbifs);
  xml_read_from_stream (is_xml, ppath.gp_lon, pbifs);
  xml_read_from_stream (is_xml, ppath.los, pbifs);
  xml_read_from_stream (is_xml, ppath.background, pbifs);
  xml_read_from_stream (is_xml, ppath.tan_pos, pbifs);
  xml_read_from_stream (is_xml, ppath.geom_tan_pos, pbifs);
  xml_read_from_stream (is_xml, ppath.p, pbifs);
  xml_read_from_stream (is_xml, ppath.t, pbifs);
  xml_read_from_stream (is_xml, ppath.vmr, pbifs);
  xml_read_from_stream (is_xml, ppath.next_parts, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/Ppath");
}

//! Writes Ppath to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ppath   Ppath
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const Ppath& ppath,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Ppath");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, ppath.dim, pbofs, "AtmosphericDimensionality");
  xml_write_to_stream (os_xml, ppath.np, pbofs,
                       "NumberOfPositionInPropagationPath");
  xml_write_to_stream (os_xml, ppath.refraction, pbofs, "RefractionFlag");
  xml_write_to_stream (os_xml, ppath.method, pbofs, "CalculationApproach");
  xml_write_to_stream (os_xml, ppath.constant, pbofs,
                       "PropagationPathConstant");
  xml_write_to_stream (os_xml, ppath.pos, pbofs,
                       "PropagationPathPointPositions");
  xml_write_to_stream (os_xml, ppath.z, pbofs, "GeometricalAltitudes");
  xml_write_to_stream (os_xml, ppath.l_step, pbofs,
                       "PropagationPathPositionLength");
  xml_write_to_stream (os_xml, ppath.gp_p, pbofs, "PressureGridIndexPosition");
  xml_write_to_stream (os_xml, ppath.gp_lat, pbofs,
                       "LatitudeGridIndexPosition");
  xml_write_to_stream (os_xml, ppath.gp_lon, pbofs,
                       "LongitudeGridIndexPosition");
  xml_write_to_stream (os_xml, ppath.los, pbofs, "LineOfSight");
  xml_write_to_stream (os_xml, ppath.background, pbofs, "RadiativeBackground");
  xml_write_to_stream (os_xml, ppath.tan_pos, pbofs, "TangentPointPosition");
  xml_write_to_stream (os_xml, ppath.geom_tan_pos, pbofs,
                       "GeometricalTangentPointPosition");
  xml_write_to_stream (os_xml, ppath.p, pbofs, "Pressures"); 
  xml_write_to_stream (os_xml, ppath.t, pbofs, "Temperatures");
  xml_write_to_stream (os_xml, ppath.vmr, pbofs, "VMRs");
  xml_write_to_stream (os_xml, ppath.next_parts, pbofs, "NextPpathParts");

  close_tag.set_name ("/Ppath");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== RetrievalQuantity =========================================

//! Reads RetrievalQuantity from XML input stream
/*!
  \param is_xml  XML Input stream
  \param rq      RetrievalQuantity return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      RetrievalQuantity& rq,
                      bifstream *pbifs)
{
  ArtsXMLTag     tag;
  String         maintag;
  String         subtag;
  String         mode;
  Index          analytical;
  Numeric        perturbation;
  ArrayOfVector  grids;

  tag.read_from_stream (is_xml);
  tag.check_name ("RetrievalQuantity");

  xml_read_from_stream (is_xml, maintag, pbifs);
  xml_read_from_stream (is_xml, subtag, pbifs);
  xml_read_from_stream (is_xml, mode, pbifs);
  xml_read_from_stream (is_xml, analytical, pbifs);
  xml_read_from_stream (is_xml, perturbation, pbifs);
  xml_read_from_stream (is_xml, grids, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/RetrievalQuantity");

  rq = RetrievalQuantity( maintag, subtag, mode, analytical, perturbation,
                          grids );
}

//! Writes RetrievalQuantity to XML output stream
/*!
  \param os_xml  XML Output stream
  \param rq      RetrievalQuantity
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const RetrievalQuantity& rq,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("RetrievalQuantity");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, rq.MainTag(), pbofs, "MainTag");
  xml_write_to_stream (os_xml, rq.Subtag(), pbofs, "Subtag");
  xml_write_to_stream (os_xml, rq.Mode(), pbofs, "Mode");
  xml_write_to_stream (os_xml, rq.Analytical(), pbofs, "Analytical");
  xml_write_to_stream (os_xml, rq.Perturbation(), pbofs, "Perturbation");
  xml_write_to_stream (os_xml, rq.Grids(), pbofs, "Grids");

  close_tag.set_name ("/RetrievalQuantity");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== SingleScatteringData ======================================

//! Reads SingleScatteringData from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ssdata  SingleScatteringData return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      SingleScatteringData& ssdata,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index ptype;

  tag.read_from_stream (is_xml);
  tag.check_name ("SingleScatteringData");

  xml_read_from_stream (is_xml, ptype, pbifs);
  ssdata.ptype = PType (ptype);
  xml_read_from_stream (is_xml, ssdata.description, pbifs);
  xml_read_from_stream (is_xml, ssdata.f_grid, pbifs);
  xml_read_from_stream (is_xml, ssdata.T_grid, pbifs);
  xml_read_from_stream (is_xml, ssdata.za_grid, pbifs);
  xml_read_from_stream (is_xml, ssdata.aa_grid, pbifs);

  xml_read_from_stream (is_xml, ssdata.pha_mat_data, pbifs);
  if (ssdata.pha_mat_data.nlibraries () != ssdata.f_grid.nelem ())
    {
      throw runtime_error("Number of frequencies in f_grid and pha_mat_data "
                          "not matching!!!");
    }

  xml_read_from_stream (is_xml, ssdata.ext_mat_data, pbifs);
  xml_read_from_stream (is_xml, ssdata.abs_vec_data, pbifs);

  tag.read_from_stream (is_xml);
  tag.check_name ("/SingleScatteringData");
}


//! Writes SingleScatteringData to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ssdata  SingleScatteringData
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const SingleScatteringData& ssdata,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("SingleScatteringData");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, Index (ssdata.ptype), pbofs);
  xml_write_to_stream (os_xml, ssdata.description, pbofs);
  xml_write_to_stream (os_xml, ssdata.f_grid, pbofs);
  xml_write_to_stream (os_xml, ssdata.T_grid, pbofs);
  xml_write_to_stream (os_xml, ssdata.za_grid, pbofs);
  xml_write_to_stream (os_xml, ssdata.aa_grid, pbofs);
  xml_write_to_stream (os_xml, ssdata.pha_mat_data, pbofs);
  xml_write_to_stream (os_xml, ssdata.ext_mat_data, pbofs);
  xml_write_to_stream (os_xml, ssdata.abs_vec_data, pbofs);

  close_tag.set_name ("/SingleScatteringData");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== SLIData2 =====================================================
//! Reads SLIData2 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param slidata  SLIData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/

void
xml_read_from_stream (istream& is_xml,
                      SLIData2& slidata,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  
  tag.read_from_stream (is_xml);
  tag.check_name ("SLIData2");

  xml_read_from_stream (is_xml, slidata.x1a, pbifs);
  xml_read_from_stream (is_xml, slidata.x2a, pbifs);
  xml_read_from_stream (is_xml, slidata.ya, pbifs);
  
  tag.read_from_stream (is_xml);
  tag.check_name ("/SLIData2");
}


void
xml_write_to_stream (ostream& os_xml,
                     const SLIData2& slidata,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("SLIData2");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  xml_write_to_stream (os_xml, slidata.x1a, pbofs);
  xml_write_to_stream (os_xml, slidata.x2a, pbofs);
  xml_write_to_stream (os_xml, slidata.ya, pbofs);
  
  close_tag.set_name ("/SLIData2");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}

//=== SpeciesRecord ================================================

//! Reads SpeciesRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param srecord  SpeciesRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      SpeciesRecord& srecord,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  String sname;
  Index degfr;
  Array<IsotopeRecord> airecord;

  tag.read_from_stream (is_xml);
  tag.check_name ("SpeciesRecord");

  xml_read_from_stream (is_xml, sname, pbifs);
  xml_read_from_stream (is_xml, degfr, pbifs);
  xml_read_from_stream (is_xml, airecord, pbifs);

  srecord = SpeciesRecord (sname.c_str (), degfr, airecord);

  tag.read_from_stream (is_xml);
  tag.check_name ("/SpeciesRecord");
}


//! Writes SpeciesRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param srecord  SpeciesRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const SpeciesRecord& srecord,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("SpeciesRecord");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  xml_write_to_stream (os_xml, srecord.Name (), pbofs);
  xml_write_to_stream (os_xml, srecord.Degfr (), pbofs);
  xml_write_to_stream (os_xml, srecord.Isotope (), pbofs);

  close_tag.set_name ("/SpeciesRecord");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


//=== SpeciesTag ================================================

//! Reads SpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param stag    SpeciesTag return value
*/
/*  param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
                 Ignored because SpeciesTag is always stored in ASCII format.*/
void
xml_read_from_stream (istream& is_xml,
                      SpeciesTag& stag,
                      bifstream * /* pbifs */)
{
  ArtsXMLTag tag;
  stringbuf  strbuf;
  char dummy;

  tag.read_from_stream (is_xml);
  tag.check_name ("SpeciesTag");

  // Skip whitespaces
  bool string_starts_with_quotes = true;
  do
    {
      is_xml >> dummy;
      switch (dummy)
        {
      case ' ':
      case '\"':
      case '\n':
      case '\r':
      case '\t':
        break;
      default:
        string_starts_with_quotes = false;
        }
    } while (is_xml.good () && dummy != '"' && string_starts_with_quotes);

  // Throw exception if first char after whitespaces is not a quote
  if (!string_starts_with_quotes)
    {
      xml_parse_error ("SpeciesTag must begin with \"");
    }

  is_xml.get (strbuf, '"');
  if (is_xml.fail ())
    {
      xml_parse_error ("SpeciesTag must end with \"");
    }

  stag = SpeciesTag (strbuf.str ());

  // Ignore quote
  is_xml >> dummy;

  tag.read_from_stream (is_xml);
  tag.check_name ("/SpeciesTag");
}


//! Writes SpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param stag    SpeciesTag
  \param name    Optional name attribute
*/
/*  param pbofs   Pointer to binary file stream. NULL for ASCII output.
                 Ignore because SpeciesTag is always stored in ASCII format. */
void
xml_write_to_stream (ostream& os_xml,
                     const SpeciesTag& stag,
                     bofstream * /* pbofs */,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("SpeciesTag");
  if (name.length ())
    open_tag.add_attribute ("name", name);
  open_tag.write_to_stream (os_xml);

  os_xml << '\"' << stag.Name () << '\"';

  close_tag.set_name ("/SpeciesTag");
  close_tag.write_to_stream (os_xml);
  os_xml << '\n';
}


////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// FIXME: These should be implemented, sooner or later...

void
xml_read_from_stream (istream&,
                      Agenda&,
                      bifstream * /* pbifs */)
{
  throw runtime_error("Method not implemented!");
}

void
xml_write_to_stream (ostream&,
                     const Agenda&,
                     bofstream * /* pbofs */,
                     const String& /* name */)
{
  throw runtime_error("Method not implemented!");
}

//=== MCAntenna ================================================

void
xml_read_from_stream (istream&,
                      MCAntenna&,
                      bifstream * /* pbifs */)
{
  throw runtime_error("Method not implemented!");
}

void
xml_write_to_stream (ostream&,
                     const MCAntenna&,
                     bofstream * /* pbofs */,
                     const String& /* name */)
{
  throw runtime_error("Method not implemented!");
}

