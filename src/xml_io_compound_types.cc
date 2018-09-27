/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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
#include "global_data.h"
#include "cloudbox.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////


//=== CIARecord ================================================

//! Reads CIARecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param irecord  SpeciesRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          CIARecord& cr,
                          bifstream* pbifs,
                          const Verbosity& verbosity)
{
    ArtsXMLTag    tag(verbosity);
    String        name;
    String molecule1;
    String molecule2;
    Index species1;
    Index species2;

    tag.read_from_stream(is_xml);
    tag.check_name("CIARecord");
    tag.get_attribute_value("molecule1", molecule1);
    tag.get_attribute_value("molecule2", molecule2);

    species1 = species_index_from_species_name(molecule1);
    species2 = species_index_from_species_name(molecule2);

    if (species1 == -1)
    {
        ostringstream os;
        os << "Unknown species (1st molecule) in CIARecord: " << molecule1;
        throw runtime_error(os.str());
    }
    if (species2 == -1)
    {
        ostringstream os;
        os << "Unknown species (2nd molecule) in CIARecord: " << molecule2;
        throw runtime_error(os.str());
    }

    cr.SetSpecies(species1, species2);

    xml_read_from_stream(is_xml, cr.mdata, pbifs, verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("/CIARecord");
}


//! Writes CIARecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param irecord  CIARecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const CIARecord& cr,
                         bofstream* pbofs,
                         const String& name _U_, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("CIARecord");
  open_tag.add_attribute("molecule1", cr.MoleculeName(0));
  open_tag.add_attribute("molecule2", cr.MoleculeName(1));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, cr.Data(), pbofs, "", verbosity);

  close_tag.set_name("/CIARecord");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//! Reads CovarianceMatrix from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    CovarianceMatrix
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
  \param verbosity
*/
void xml_read_from_stream(istream& is_xml,
                          CovarianceMatrix& covmat,
                          bifstream* pbifs,
                          const Verbosity& verbosity)
{
    ArtsXMLTag    tag(verbosity);
    String        name, type;
    Index n_blocks, row_start, row_extent, column_start, column_extent, row_index, column_index,
        is_inverse;

    tag.read_from_stream(is_xml);
    tag.check_name("CovarianceMatrix");
    tag.get_attribute_value("n_blocks", n_blocks);

    for (Index i = 0; i < n_blocks; i++) {
        tag.read_from_stream(is_xml);
        tag.check_name("Block");

        tag.get_attribute_value("row_index",  row_index);
        tag.get_attribute_value("column_index",  column_index);
        tag.get_attribute_value("row_start",  row_start);
        tag.get_attribute_value("row_extent", row_extent);
        tag.get_attribute_value("column_start",  column_start);
        tag.get_attribute_value("column_extent", column_extent);
        tag.get_attribute_value("type", type);
        tag.get_attribute_value("is_inverse", is_inverse);

        Range row_range(row_start, row_extent);
        Range column_range(column_start, column_extent);
        if (type == "Matrix") {
            std::shared_ptr<Matrix> M = std::make_shared<Matrix>(row_extent, column_extent);
            xml_read_from_stream(is_xml, *M, pbifs, verbosity);
            if (!is_inverse) {
                covmat.correlations_.emplace_back(row_range, column_range,
                                                  std::make_pair(row_index, column_index), M);
            } else {
                covmat.inverses_.emplace_back(row_range, column_range,
                                           std::make_pair(row_index, column_index), M);
            }
        } else if (type == "Sparse") {
            std::shared_ptr<Sparse> M = std::make_shared<Sparse>(row_extent, column_extent);
            xml_read_from_stream(is_xml, *M, pbifs, verbosity);
            if (!is_inverse) {
                covmat.correlations_.emplace_back(row_range, column_range,
                                                  std::make_pair(row_index, column_index), M);
            } else {
                covmat.inverses_.emplace_back(row_range, column_range,
                                              std::make_pair(row_index, column_index), M);
            }
        }
        tag.read_from_stream(is_xml);
        tag.check_name("/Block");
    }
    tag.read_from_stream(is_xml);
    tag.check_name("/CovarianceMatrix");
}

//=== CovarianceMatrix =========================================================

//! Write CovarianceMatrix to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    CovarianceMatrix
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
  \param verbosity
*/
void xml_write_to_stream(ostream& os_xml,
                         const CovarianceMatrix& covmat,
                         bofstream* pbofs,
                         const String& name _U_, const Verbosity& verbosity)
{
    ArtsXMLTag covmat_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);


    covmat_tag.set_name("CovarianceMatrix");
    covmat_tag.add_attribute("n_blocks", covmat.correlations_.size() + covmat.inverses_.size());
    covmat_tag.write_to_stream(os_xml);
    os_xml << '\n';
    for (const Block &c : covmat.correlations_) {
        ArtsXMLTag block_tag(verbosity);
        block_tag.set_name("Block");

        Index i,j;
        std::tie(i,j) = c.get_indices();
        block_tag.add_attribute("row_index", i);
        block_tag.add_attribute("column_index", j);

        Range row_range = c.get_row_range();
        Range column_range = c.get_column_range();
        block_tag.add_attribute("row_start", row_range.get_start());
        block_tag.add_attribute("row_extent", row_range.get_extent());
        block_tag.add_attribute("column_start", column_range.get_start());
        block_tag.add_attribute("column_extent", column_range.get_extent());
        block_tag.add_attribute("is_inverse", 0);
        if (c.get_matrix_type() == Block::MatrixType::dense) {
            block_tag.add_attribute("type", "Matrix");
            block_tag.write_to_stream(os_xml);
            os_xml << '\n';
            xml_write_to_stream(os_xml, c.get_dense(), pbofs, name, verbosity);
        } else {
            block_tag.add_attribute("type", "Sparse");
            block_tag.write_to_stream(os_xml);
            os_xml << '\n';
            xml_write_to_stream(os_xml, c.get_sparse(), pbofs, name, verbosity);
        }
        close_tag.set_name("/Block");
        close_tag.write_to_stream(os_xml);
        os_xml << '\n';
    }
    for (const Block &c : covmat.inverses_) {
        ArtsXMLTag block_tag(verbosity);
        block_tag.set_name("Block");

        Index i,j;
        std::tie(i,j) = c.get_indices();
        block_tag.add_attribute("row_index", i);
        block_tag.add_attribute("column_index", j);

        Range row_range = c.get_row_range();
        Range column_range = c.get_column_range();
        block_tag.add_attribute("row_start", row_range.get_start());
        block_tag.add_attribute("row_extent", row_range.get_extent());
        block_tag.add_attribute("column_start", column_range.get_start());
        block_tag.add_attribute("column_extent", column_range.get_extent());
        block_tag.add_attribute("is_inverse", 1);
        if (c.get_matrix_type() == Block::MatrixType::dense) {
            block_tag.add_attribute("type", "Matrix");
            block_tag.write_to_stream(os_xml);
            os_xml << '\n';
            xml_write_to_stream(os_xml, c.get_dense(), pbofs, name, verbosity);
        } else {
            block_tag.add_attribute("type", "Sparse");
            block_tag.write_to_stream(os_xml);
            os_xml << '\n';
            xml_write_to_stream(os_xml, c.get_sparse(), pbofs, name, verbosity);
        }
        close_tag.set_name("/Block");
        close_tag.write_to_stream(os_xml);
        os_xml << '\n';
    }
    os_xml << '\n';
    close_tag.set_name("/CovarianceMatrix");
    close_tag.write_to_stream(os_xml);
}

//=== GasAbsLookup ===========================================================

//! Reads GasAbsLookup from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gal     GasAbsLookup return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GasAbsLookup& gal,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GasAbsLookup");

  xml_read_from_stream(is_xml, gal.species, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.nonlinear_species, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.f_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.p_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.vmrs_ref, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.t_ref, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.t_pert, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.nls_pert, pbifs, verbosity);
  xml_read_from_stream(is_xml, gal.xsec, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GasAbsLookup");
}


//! Writes GasAbsLookup to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gal     GasAbsLookup
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GasAbsLookup& gal,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GasAbsLookup");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, gal.species, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gal.nonlinear_species, pbofs,
                      "NonlinearSpecies", verbosity);
  xml_write_to_stream(os_xml, gal.f_grid, pbofs, "FrequencyGrid", verbosity);
  xml_write_to_stream(os_xml, gal.p_grid, pbofs, "PressureGrid", verbosity);
  xml_write_to_stream(os_xml, gal.vmrs_ref, pbofs, "ReferenceVmrProfiles", verbosity);
  xml_write_to_stream(os_xml, gal.t_ref, pbofs, "ReferenceTemperatureProfile", verbosity);
  xml_write_to_stream(os_xml, gal.t_pert, pbofs, "TemperaturePerturbations", verbosity);
  xml_write_to_stream(os_xml, gal.nls_pert, pbofs,
                      "NonlinearSpeciesVmrPerturbations", verbosity);
  xml_write_to_stream(os_xml, gal.xsec, pbofs, "AbsorptionCrossSections", verbosity);

  close_tag.set_name("/GasAbsLookup");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField ===========================================================

//! Reads the grids for gridded fields from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  for (Index i = 0; i < gfield.get_dim(); i++)
    {
      tag.read_from_stream(is_xml);
      if (tag.get_name() == "Vector")
        {
          String s;
          tag.get_attribute_value("name", s);
          if (s.length())
            gfield.set_grid_name(i, s);

          Vector v;
          xml_parse_from_stream(is_xml, v, pbifs, tag, verbosity);
          gfield.set_grid(i, v);
          tag.read_from_stream(is_xml);
          tag.check_name("/Vector");
        }
      else if (tag.get_name() == "Array")
        {
          String s;
          tag.get_attribute_value("name", s);
          if (s.length())
            gfield.set_grid_name(i, s);

          tag.get_attribute_value("type", s);
          if (s == "String")
            {
              ArrayOfString as;
              xml_parse_from_stream(is_xml, as, pbifs, tag, verbosity);
              gfield.set_grid(i, as);
              tag.read_from_stream(is_xml);
              tag.check_name("/Array");

            }
          else
            {
              xml_parse_error("Grids must be of type *Vector* or *ArrayOfString*\n"
                              "but *ArrayOf" + s + "* found.");
            }
        }
      else
        {
          ostringstream os;
          os <<"Grids must be of type *Vector* or *ArrayOfString*\n"
            << "but tag <"+ tag.get_name() + "> found.";
          if (tag.get_name() == "ArrayOfString")
            os << "\nCorrect XML tag for *ArrayOfString* is <Array type=\"String\" ...>.";
          xml_parse_error(os.str());
        }
    }
}


//! Writes the grids for gridded fields to an XML input stream
/*!
  \param os_xml  XML output stream
  \param gfield  GriddedField with the grids
  \param pbofs   Pointer to binary output stream. NULL in case of ASCII file.
*/
void xml_write_to_stream(ostream&            os_xml,
                         const GriddedField& gfield,
                         bofstream*          pbofs,
                         const               String& /* name */,
                         const Verbosity&    verbosity)
{
  for (Index i = 0; i < gfield.get_dim(); i++)
    {
      switch (gfield.get_grid_type(i))
        {
        case GRID_TYPE_NUMERIC:
           xml_write_to_stream(os_xml, gfield.get_numeric_grid(i),
                               pbofs, gfield.get_grid_name(i), verbosity);
          break;
        case GRID_TYPE_STRING:
          xml_write_to_stream(os_xml, gfield.get_string_grid(i),
                              pbofs, gfield.get_grid_name(i), verbosity);
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
void xml_read_from_stream(istream& is_xml,
                          GriddedField1& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField1");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField1");

  gfield.checksize_strict();
}


//! Writes GriddedField1 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField1
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField1& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField1");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField1");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField2 ===========================================================

//! Reads GriddedField2 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField2 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField2& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField2");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField2");

  gfield.checksize_strict();
}


//! Writes GriddedField2 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField2
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField2& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField2");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField2");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField3 ===========================================================

//! Reads GriddedField3 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField3 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField3& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField3");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField3");

  gfield.checksize_strict();
}


//! Writes GriddedField3 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField3
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField3& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField3");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField3");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField4 ===========================================================

//! Reads GriddedField4 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField4 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField4& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField4");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField4");

  gfield.checksize_strict();
}


//! Writes GriddedField4 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField4
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField4& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField4");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField4");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField5 ===========================================================

//! Reads GriddedField5 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField5 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField5& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField5");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField5");

  gfield.checksize_strict();
}


//! Writes GriddedField5 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField5
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField5& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField5");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField5");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GriddedField6 ===========================================================

//! Reads GriddedField6 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField6 return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GriddedField6& gfield,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField6");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length())
    gfield.set_name(s);

  xml_read_from_stream(is_xml, *((GriddedField*)&gfield), pbifs, verbosity);
  xml_read_from_stream(is_xml, gfield.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GriddedField6");

  gfield.checksize_strict();
}


//! Writes GriddedField6 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gfield  GriddedField6
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField6& gfield,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GriddedField6");
  if (!name.length() && (gfield.get_name().length()))
    open_tag.add_attribute("name", gfield.get_name());
  else if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, *((GriddedField*)&gfield), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data", verbosity);

  close_tag.set_name("/GriddedField6");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== GridPos =====================================================

//! Reads GridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gpos    GridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          GridPos& gpos,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GridPos");

  xml_read_from_stream(is_xml, gpos.idx, pbifs, verbosity);
  xml_read_from_stream(is_xml, gpos.fd[0], pbifs, verbosity);
  xml_read_from_stream(is_xml, gpos.fd[1], pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/GridPos");
}


//! Writes GridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gpos    GridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const GridPos& gpos,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GridPos");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, gpos.idx, pbofs,
                      "OriginalGridIndexBelowInterpolationPoint", verbosity);
  xml_write_to_stream(os_xml, gpos.fd[0], pbofs,
                      "FractionalDistanceToNextPoint_1", verbosity);
  xml_write_to_stream(os_xml, gpos.fd[1], pbofs,
                      "FractionalDistanceToNextPoint_2", verbosity);

  close_tag.set_name("/GridPos");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== IsotopologueRecord ================================================

//! Reads IsotopologueRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param irecord  SpeciesRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          IsotopologueRecord& irecord,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag    tag(verbosity);
  String        name;
  Numeric       abundance;
  Numeric       mass;
  Index         mytrantag;
  Index         hitrantag;
  ArrayOfIndex  jpltags;

  tag.read_from_stream(is_xml);
  tag.check_name("IsotopologueRecord");

  xml_read_from_stream(is_xml, name, pbifs, verbosity);
  xml_read_from_stream(is_xml, abundance, pbifs, verbosity);
  xml_read_from_stream(is_xml, mass, pbifs, verbosity);
  xml_read_from_stream(is_xml, mytrantag, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitrantag, pbifs, verbosity);
  xml_read_from_stream(is_xml, jpltags, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/IsotopologueRecord");

  irecord = IsotopologueRecord(name, abundance, mass, mytrantag, hitrantag,
                          jpltags);
}


//! Writes IsotopologueRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param irecord  IsotopologueRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const IsotopologueRecord& irecord,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("IsotopologueRecord");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, irecord.Name(), pbofs, "Name", verbosity);
  xml_write_to_stream(os_xml, irecord.Abundance(), pbofs, "Abundance", verbosity);
  xml_write_to_stream(os_xml, irecord.Mass(), pbofs, "Mass", verbosity);
  xml_write_to_stream(os_xml, irecord.MytranTag(), pbofs, "MytranTag", verbosity);
  xml_write_to_stream(os_xml, irecord.HitranTag(), pbofs, "HitranTag", verbosity);
  xml_write_to_stream(os_xml, irecord.JplTags(), pbofs, "JplTags", verbosity);

  close_tag.set_name("/IsotopologueRecord");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== LineMixingRecord ================================================

//! Reads LineMixingRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param lmr      LineMixingRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          LineMixingRecord& lmr,
                          bifstream* pbifs, const Verbosity& verbosity)
{
    ArtsXMLTag    tag(verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("LineMixingRecord");

    SpeciesTag species_tag;

    xml_read_from_stream(is_xml, species_tag, pbifs, verbosity);
    lmr = LineMixingRecord(species_tag.Species(), species_tag.Isotopologue());
    xml_read_from_stream(is_xml, lmr.Quantum(), pbifs, verbosity);
    xml_read_from_stream(is_xml, lmr.Data(), pbifs, verbosity);
    
    tag.read_from_stream(is_xml);
    tag.check_name("/LineMixingRecord");
}


//! Writes LineMixingRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param lmr      LineMixingRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const LineMixingRecord& lmr,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("LineMixingRecord");
    if (name.length())
        open_tag.add_attribute("name", name);
    open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    using global_data::species_data;
    String species_tag = species_name_from_species_index(lmr.Species()) + '-'
        + species_data[lmr.Species()].Isotopologue()[lmr.Isotopologue()].Name();

    xml_write_to_stream(os_xml, SpeciesTag(species_tag), pbofs, "", verbosity);
    xml_write_to_stream(os_xml, lmr.Quantum(), pbofs, "", verbosity);
    xml_write_to_stream(os_xml, lmr.Data(), pbofs, "", verbosity);

    close_tag.set_name("/LineMixingRecord");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';
}



//=== Ppath =====================================================

//! Reads Ppath from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ppath   Ppath return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Ppath& ppath,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Ppath");

  xml_read_from_stream(is_xml, ppath.dim, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.np, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.constant, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.background, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.start_pos, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.start_los, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.start_lstep, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.pos, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.los, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.r, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.lstep, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.end_pos, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.end_los, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.end_lstep, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.nreal, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.ngroup, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.gp_p, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.gp_lat, pbifs, verbosity);
  xml_read_from_stream(is_xml, ppath.gp_lon, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/Ppath");
}


//! Writes Ppath to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ppath   Ppath
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Ppath& ppath,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Ppath");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, ppath.dim, pbofs, "AtmosphericDimensionality", verbosity);
  xml_write_to_stream(os_xml, ppath.np, pbofs,
                      "NumberOfPositionInPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.constant, pbofs,
                      "PropagationPathConstant", verbosity);
  xml_write_to_stream(os_xml, ppath.background, pbofs, "RadiativeBackground", verbosity);
  xml_write_to_stream(os_xml, ppath.start_pos, pbofs,
                      "StartPositionOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.start_los, pbofs,
                      "StartLOSOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.start_lstep, pbofs,
                      "StartLstepOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.pos, pbofs,
                      "PropagationPathPointPositions", verbosity);
  xml_write_to_stream(os_xml, ppath.los, pbofs, "LineOfSight", verbosity);
  xml_write_to_stream(os_xml, ppath.r, pbofs, "PropagationPathPointRadii", verbosity);
  xml_write_to_stream(os_xml, ppath.lstep, pbofs,
                      "PropagationPathPositionLength", verbosity);
  xml_write_to_stream(os_xml, ppath.end_pos, pbofs,
                      "EndPositionOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.end_los, pbofs,
                      "EndLOSOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.end_lstep, pbofs,
                      "EndLstepPropagationPath", verbosity);
  xml_write_to_stream(os_xml, ppath.nreal, pbofs, "RefractiveIndexRealPart", verbosity);
  xml_write_to_stream(os_xml, ppath.ngroup, pbofs, "GroupRefractiveIndex", verbosity);
  xml_write_to_stream(os_xml, ppath.gp_p, pbofs, "PressureGridIndexPosition", verbosity);
  xml_write_to_stream(os_xml, ppath.gp_lat, pbofs,
                      "LatitudeGridIndexPosition", verbosity);
  xml_write_to_stream(os_xml, ppath.gp_lon, pbofs,
                      "LongitudeGridIndexPosition", verbosity);

  close_tag.set_name("/Ppath");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== PropagationMatrix ======================================================

//! Reads PropagationMatrix from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param pm         PropagationMatrix return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream&              is_xml,
                          PropagationMatrix&    pm,
                          bifstream*            pbifs,
                          const Verbosity&      verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("PropagationMatrix");

  try
  {
    Tensor4 d;
    xml_read_from_stream(is_xml, d, pbifs, verbosity);
    pm = PropagationMatrix(d);
  }
  catch (const std::runtime_error &e)
  {
    ostringstream os;
    os << "Error reading PropagationMatrix: "
    << "\n" << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/PropagationMatrix");
}


//! Writes PropagationMatrix to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param pm         PropagationMatrix
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream&                    os_xml,
                         const PropagationMatrix&    pm,
                         bofstream*                  pbofs,
                         const String&               name,
                         const Verbosity&            verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("PropagationMatrix");
  if (name.length())
    open_tag.add_attribute("name", name);
  
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, pm.GetData(), pbofs, "", verbosity);

  close_tag.set_name("/PropagationMatrix");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}


//=== QuantumIdentifier =========================================

//! Reads QuantumIdentifier from XML input stream
/*!
  \param is_xml  XML Input stream
  \param qi      QuantumIdentifier return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          QuantumIdentifier& qi,
                          bifstream* pbifs _U_,
                          const Verbosity& verbosity)
{
    ArtsXMLTag tag(verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("QuantumIdentifier");

    try
    {
        String qi_str;
        parse_xml_tag_content_as_string(is_xml, qi_str);
        qi.SetFromString(qi_str);
    }
    catch (const std::runtime_error &e)
    {
        ostringstream os;
        os << "Error reading QuantumIdentifier: "
        << "\n" << e.what();
        throw runtime_error(os.str());
    }

    tag.read_from_stream(is_xml);
    tag.check_name("/QuantumIdentifier");
}


//! Writes QuantumIdentifier to XML output stream
/*!
  \param os_xml  XML Output stream
  \param qi      QuantumIdentifier
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const QuantumIdentifier& qi,
                         bofstream* pbofs _U_,
                         const String& name, const Verbosity& verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("QuantumIdentifier");
    if (name.length())
        open_tag.add_attribute("name", name);
    open_tag.write_to_stream(os_xml);

    os_xml << qi;

    close_tag.set_name("/QuantumIdentifier");
    close_tag.write_to_stream(os_xml);
    os_xml << endl;
}


//=== QuantumNumberRecord =========================================

//! Reads QuantumNumberRecord from XML input stream
/*!
  \param is_xml  XML Input stream
  \param qnr      QuantumNumberRecord return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          QuantumNumberRecord& qnr,
                          bifstream* pbifs _U_,
                          const Verbosity& verbosity)
{
    ArtsXMLTag tag(verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("QuantumNumberRecord");

    try
    {
        tag.read_from_stream(is_xml);
        tag.check_name("Upper");
        xml_read_from_stream(is_xml, qnr.Upper(), pbifs, verbosity);
        tag.read_from_stream(is_xml);
        tag.check_name("/Upper");
    }
    catch (const std::runtime_error &e)
    {
        ostringstream os;
        os << "Error in upper quantum numbers while reading QuantumNumberRecord: "
        << "\n" << e.what();
        throw runtime_error(os.str());
    }

    try
    {
        tag.read_from_stream(is_xml);
        tag.check_name("Lower");
        xml_read_from_stream(is_xml, qnr.Lower(), pbifs, verbosity);
        tag.read_from_stream(is_xml);
        tag.check_name("/Lower");
    }
    catch (const std::runtime_error &e)
    {
        ostringstream os;
        os << "Error in lower quantum numbers while reading QuantumNumberRecord: "
        << "\n" << e.what();
        throw runtime_error(os.str());
    }

    tag.read_from_stream(is_xml);
    tag.check_name("/QuantumNumberRecord");
}


//! Writes QuantumNumberRecord to XML output stream
/*!
  \param os_xml  XML Output stream
  \param qnr      QuantumNumberRecord
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const QuantumNumberRecord& qnr,
                         bofstream* pbofs _U_,
                         const String& name, const Verbosity& verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("QuantumNumberRecord");
    if (name.length())
        open_tag.add_attribute("name", name);
    open_tag.write_to_stream(os_xml);

    os_xml << '\n';
    open_tag.set_name("Upper");
    open_tag.write_to_stream(os_xml);
    xml_write_to_stream(os_xml, qnr.Upper(), pbofs, "", verbosity);
    close_tag.set_name("/Upper");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';

    open_tag.set_name("Lower");
    open_tag.write_to_stream(os_xml);
    xml_write_to_stream(os_xml, qnr.Lower(), pbofs, "", verbosity);
    close_tag.set_name("/Lower");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';

    close_tag.set_name("/QuantumNumberRecord");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';
}


//=== QuantumNumbers =========================================

//! Reads QuantumNumbers from XML input stream
/*!
  \param is_xml  XML Input stream
  \param qn      QuantumNumbers return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          QuantumNumbers& qn,
                          bifstream* pbifs _U_,
                          const Verbosity& verbosity)
{
    ArtsXMLTag tag(verbosity);
    Index nelem;

    tag.read_from_stream(is_xml);
    tag.check_name("QuantumNumbers");

    tag.get_attribute_value("nelem", nelem);

    Index n;
    try
    {
        for (n = 0; n < nelem; n++)
            is_xml >> qn;
    }
    catch (const std::runtime_error &e)
    {
        ostringstream os;
        os << "Error reading QuantumNumbers: "
        << "\n Element: " << n
        << "\n" << e.what();
        throw runtime_error(os.str());
    }

    tag.read_from_stream(is_xml);
    tag.check_name("/QuantumNumbers");
}


//! Writes QuantumNumbers to XML output stream
/*!
  \param os_xml  XML Output stream
  \param qn      QuantumNumbers
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const QuantumNumbers& qn,
                         bofstream* pbofs _U_,
                         const String& name, const Verbosity& verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("QuantumNumbers");
    if (name.length())
        open_tag.add_attribute("name", name);
    open_tag.add_attribute("nelem", qn.GetNumbers().size());
    open_tag.write_to_stream(os_xml);

    os_xml << " " << qn << " ";

    close_tag.set_name("/QuantumNumbers");
    close_tag.write_to_stream(os_xml);
}


//=== RetrievalQuantity =========================================

//! Reads RetrievalQuantity from XML input stream
/*!
  \param is_xml  XML Input stream
  \param rq      RetrievalQuantity return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          RetrievalQuantity& rq,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag     tag(verbosity);
  String         maintag;
  String         subtag;
  String         subsubtag;
  String         mode;
  Index          analytical;
  Numeric        perturbation;
  ArrayOfVector  grids;

  tag.read_from_stream(is_xml);
  tag.check_name("RetrievalQuantity");

  xml_read_from_stream(is_xml, maintag, pbifs, verbosity);
  xml_read_from_stream(is_xml, subtag, pbifs, verbosity);
  xml_read_from_stream(is_xml, subsubtag, pbifs, verbosity);
  xml_read_from_stream(is_xml, mode, pbifs, verbosity);
  xml_read_from_stream(is_xml, analytical, pbifs, verbosity);
  xml_read_from_stream(is_xml, perturbation, pbifs, verbosity);
  xml_read_from_stream(is_xml, grids, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/RetrievalQuantity");

  rq = RetrievalQuantity(maintag, subtag, subsubtag, mode, analytical,
                         perturbation, grids);
}


//! Writes RetrievalQuantity to XML output stream
/*!
  \param os_xml  XML Output stream
  \param rq      RetrievalQuantity
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const RetrievalQuantity& rq,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("RetrievalQuantity");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, rq.MainTag(), pbofs, "MainTag", verbosity);
  xml_write_to_stream(os_xml, rq.Subtag(), pbofs, "Subtag", verbosity);
  xml_write_to_stream(os_xml, rq.SubSubtag(), pbofs, "SubSubtag", verbosity);
  xml_write_to_stream(os_xml, rq.Mode(), pbofs, "Mode", verbosity);
  xml_write_to_stream(os_xml, rq.Analytical(), pbofs, "Analytical", verbosity);
  xml_write_to_stream(os_xml, rq.Perturbation(), pbofs, "Perturbation", verbosity);
  xml_write_to_stream(os_xml, rq.Grids(), pbofs, "Grids", verbosity);

  close_tag.set_name("/RetrievalQuantity");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== SingleScatteringData ======================================

//! Reads SingleScatteringData from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ssdata  SingleScatteringData return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          SingleScatteringData& ssdata,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("SingleScatteringData");
  tag.get_attribute_value("version", version);

  if (version == "3")
    {
      String ptype_string;
      xml_read_from_stream(is_xml, ptype_string, pbifs, verbosity);
      ssdata.ptype = PTypeFromString(ptype_string);
    }
  else if (version == "2")
    {
      String ptype_string;
      xml_read_from_stream(is_xml, ptype_string, pbifs, verbosity);
      ssdata.ptype = PType2FromString(ptype_string);
    }
  else
    {
      Index ptype;
      xml_read_from_stream(is_xml, ptype, pbifs, verbosity);
      if (ptype != PTYPE_GENERAL &&
          ptype != PTYPE_TOTAL_RND &&
          ptype != PTYPE_AZIMUTH_RND)
      {
            ostringstream os;
            os << "Ptype value (" << ptype << ") is wrong."
            << "It must be \n"
            << PTYPE_TOTAL_RND << " - totally randomly oriented particles,\n"
            << PTYPE_AZIMUTH_RND << " - azimuthally randomly oriented particles, or\n"
            << PTYPE_GENERAL << " - arbitrary oriented particles.\n";
            throw runtime_error( os.str() );
      }
      ssdata.ptype = PType(ptype);
    }
  xml_read_from_stream(is_xml, ssdata.description, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.f_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.T_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.za_grid, pbifs, verbosity);
  /* Verify that we have a good coverage for the za grid */
  if ((ssdata.za_grid[0] > 1) || ssdata.za_grid[ssdata.za_grid.nelem()-1] < 179)
    {
      ostringstream os;
      os << "Missing data in xml-stream. Expected za_grid: [0, 180]. "
         << "Found za_grid: [" << ssdata.za_grid[0]
         << ", "
         << ssdata.za_grid[ssdata.za_grid.nelem()-1]
         << "]";
      throw runtime_error(os.str());
    }
  xml_read_from_stream(is_xml, ssdata.aa_grid, pbifs, verbosity);

  xml_read_from_stream(is_xml, ssdata.pha_mat_data, pbifs, verbosity);
  if (ssdata.pha_mat_data.nlibraries() != ssdata.f_grid.nelem())
    {
      throw runtime_error("Number of frequencies in f_grid and pha_mat_data "
                          "not matching!!!");
    }

  xml_read_from_stream(is_xml, ssdata.ext_mat_data, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.abs_vec_data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/SingleScatteringData");

  if (version != "3" && ssdata.ptype == PTYPE_AZIMUTH_RND)
  {
      ConvertAzimuthallyRandomSingleScatteringData(ssdata);
  }

  chk_scat_data(ssdata, verbosity);
}


//! Writes SingleScatteringData to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ssdata  SingleScatteringData
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const SingleScatteringData& ssdata,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SingleScatteringData");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  os_xml << '\n';
  xml_write_to_stream(os_xml, PTypeToString(ssdata.ptype),
                      pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.description, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.f_grid, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.T_grid, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.za_grid, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.aa_grid, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.pha_mat_data, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.ext_mat_data, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, ssdata.abs_vec_data, pbofs, "", verbosity);

  close_tag.set_name("/SingleScatteringData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== ScatteringMetaData ======================================

//! Reads ScatteringMetaData from XML input stream
/*!
  \param is_xml  XML Input stream
  \param smdata  ScatteringMetaData return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ScatteringMetaData& smdata,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("ScatteringMetaData");
  tag.get_attribute_value("version", version);

  if (version != "3")
    {
        ostringstream os;
        os << "Only ScatteringMetaData version 3 can be handled. "
           << "Versions 1 and 2 are obsolete.";
        throw runtime_error(os.str());
    }
      
  xml_read_from_stream(is_xml, smdata.description, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.source, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.refr_index, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.mass, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.diameter_max, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.diameter_volume_equ, pbifs, verbosity);
  xml_read_from_stream(is_xml, smdata.diameter_area_equ_aerodynamical, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/ScatteringMetaData");
}


//! Writes ScatteringMetaData to XML output stream
/*!
  \param os_xml  XML Output stream
  \param smdata  ScatteringMetaData
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ScatteringMetaData& smdata,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("ScatteringMetaData");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, smdata.description, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.source, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.refr_index, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.mass, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.diameter_max, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.diameter_volume_equ, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.diameter_area_equ_aerodynamical, pbofs, "", verbosity);
  
  close_tag.set_name("/ScatteringMetaData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}



//=== SLIData2 =====================================================
//! Reads SLIData2 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param slidata  SLIData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/

void xml_read_from_stream(istream& is_xml,
                          SLIData2& slidata,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("SLIData2");

  xml_read_from_stream(is_xml, slidata.x1a, pbifs, verbosity);
  xml_read_from_stream(is_xml, slidata.x2a, pbifs, verbosity);
  xml_read_from_stream(is_xml, slidata.ya, pbifs, verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("/SLIData2");
}


void xml_write_to_stream(ostream& os_xml,
                         const SLIData2& slidata,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SLIData2");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, slidata.x1a, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, slidata.x2a, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, slidata.ya, pbofs, "", verbosity);
  
  close_tag.set_name("/SLIData2");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== SpeciesAuxData ===========================================

//! Reads SpeciesAuxData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param sap      SpeciesAuxData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream&           is_xml,
                          SpeciesAuxData&    sad,
                          bifstream*         pbifs _U_,
                          const Verbosity&   verbosity)
{
    CREATE_OUT2;

    ArtsXMLTag tag(verbosity);
    Index nelem;
    Index nparam;

    tag.read_from_stream(is_xml);
    tag.check_name("SpeciesAuxData");

    Index version;
    tag.get_attribute_value("version", version);

    if (version == 1)
    {
        tag.get_attribute_value("nelem", nelem);
        tag.get_attribute_value("nparam", nparam);

        Index n = 0;
        try
        {
            ArrayOfString artstags;
            sad.InitFromSpeciesData();
            for (n = 0; n < nelem; n++)
            {
                String artstag;
                sad.ReadFromStream(artstag, is_xml, nparam, verbosity);

                if (find_first(artstags, artstag) == -1)
                    artstags.push_back(artstag);
                else
                {
                    ostringstream os;
                    os << "SpeciesAuxData for " << artstag << " already defined.\n"
                    << "Duplicates are not allowed in input file.";
                    throw runtime_error(os.str());
                }
            }
        }
        catch (const std::runtime_error &e)
        {
            ostringstream os;
            os << "Error reading SpeciesAuxData: "
            << "\n Element: " << n
            << "\n" << e.what();
            throw runtime_error(os.str());
        }
    }
    else if (version == 2)
    {
        tag.get_attribute_value("nelem", nelem);
        Index n = 0;
        try
        {
            ArrayOfString artstags;
            sad.InitFromSpeciesData();
            for (n = 0; n < nelem; n++)
            {
                String artstag;
                String auxtype;
                ArrayOfGriddedField1 auxdata;

                xml_read_from_stream(is_xml, artstag, pbifs, verbosity);
                xml_read_from_stream(is_xml, auxtype, pbifs, verbosity);
                xml_read_from_stream(is_xml, auxdata, pbifs, verbosity);

                sad.setParam(artstag, auxtype, auxdata);

                if (find_first(artstags, artstag) == -1)
                    artstags.push_back(artstag);
                else
                {
                    ostringstream os;
                    os << "SpeciesAuxData for " << artstag << " already defined.\n"
                    << "Duplicates are not allowed in input file.";
                    throw runtime_error(os.str());
                }
            }
        }
        catch (const std::runtime_error &e)
        {
            ostringstream os;
            os << "Error reading SpeciesAuxData: "
            << "\n Element: " << n
            << "\n" << e.what();
            throw runtime_error(os.str());
        }
    }
    else
    {
        ostringstream os;
        os << "Unsupported SpeciesAuxData version number: " << version << ", expected 1 or 2.";
        throw std::runtime_error(os.str());
    }

    tag.read_from_stream(is_xml);
    tag.check_name("/SpeciesAuxData");
}


//! Writes SpeciesAuxData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param sap      SpeciesAuxData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream&                 os_xml,
                         const SpeciesAuxData&    sad,
                         bofstream* pbofs         _U_,
                         const String&            name,
                         const Verbosity&         verbosity)

{
    using global_data::species_data;

    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    Index nelem = 0;
    for (Index isp = 0; isp < sad.nspecies(); isp++)
        for (Index iso = 0; iso < sad.nisotopologues(isp); iso++)
            if (sad.getParam(isp, iso).nelem())
                nelem++;

    open_tag.set_name("SpeciesAuxData");
    if (name.length())
        open_tag.add_attribute("name", name);

    open_tag.add_attribute("version", 2);
    open_tag.add_attribute("nelem", nelem);

    open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    xml_set_stream_precision(os_xml);

    for (Index isp = 0; isp < sad.nspecies(); isp++)
    {
        const String spname = species_data[isp].Name();
        for (Index iso = 0; iso < sad.nisotopologues(isp); iso++)
        {
            if (sad.getParam(isp, iso).nelem())
            {
                const String artstag =
                spname + "-" + species_data[isp].Isotopologue()[iso].Name();

                xml_write_to_stream(os_xml, artstag, pbofs, "", verbosity);
                xml_write_to_stream(os_xml, sad.getTypeString(isp, iso), pbofs, "", verbosity);
                xml_write_to_stream(os_xml, sad.getParam(isp, iso), pbofs, "", verbosity);
            }
        }
    }

    close_tag.set_name("/SpeciesAuxData");
    close_tag.write_to_stream(os_xml);
    
    os_xml << '\n';
}


//=== SpeciesRecord ================================================

//! Reads SpeciesRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param srecord  SpeciesRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          SpeciesRecord& srecord,
                          bifstream* pbifs, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);
  String sname;
  Index degfr;
  Array<IsotopologueRecord> airecord;

  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesRecord");

  xml_read_from_stream(is_xml, sname, pbifs, verbosity);
  xml_read_from_stream(is_xml, degfr, pbifs, verbosity);
  xml_read_from_stream(is_xml, airecord, pbifs, verbosity);

  srecord = SpeciesRecord(sname.c_str(), degfr, airecord);

  tag.read_from_stream(is_xml);
  tag.check_name("/SpeciesRecord");
}


//! Writes SpeciesRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param srecord  SpeciesRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const SpeciesRecord& srecord,
                         bofstream* pbofs,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SpeciesRecord");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, srecord.Name(), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, srecord.Degfr(), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, srecord.Isotopologue(), pbofs, "", verbosity);

  close_tag.set_name("/SpeciesRecord");
  close_tag.write_to_stream(os_xml);
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
void xml_read_from_stream(istream& is_xml,
                          SpeciesTag& stag,
                          bifstream* /* pbifs */, const Verbosity& verbosity)
{
  ArtsXMLTag tag(verbosity);
  stringbuf  strbuf;
  char dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesTag");

  // Skip whitespaces
  bool string_starts_with_quotes = true;
  do {
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
    } while (is_xml.good() && dummy != '"' && string_starts_with_quotes);

  // Throw exception if first char after whitespaces is not a quote
  if (!string_starts_with_quotes)
    {
      xml_parse_error("SpeciesTag must begin with \"");
    }

  is_xml.get(strbuf, '"');
  if (is_xml.fail())
    {
      xml_parse_error("SpeciesTag must end with \"");
    }

  stag = SpeciesTag(strbuf.str());

  // Ignore quote
  is_xml >> dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("/SpeciesTag");
}


//! Writes SpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param stag    SpeciesTag
  \param name    Optional name attribute
*/
/*  param pbofs   Pointer to binary file stream. NULL for ASCII output.
                 Ignore because SpeciesTag is always stored in ASCII format. */
void xml_write_to_stream(ostream& os_xml,
                         const SpeciesTag& stag,
                         bofstream* /* pbofs */,
                         const String& name, const Verbosity& verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SpeciesTag");
  if (name.length())
    open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  os_xml << '\"' << stag.Name() << '\"';

  close_tag.set_name("/SpeciesTag");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== StokesVector ======================================================

//! Reads StokesVector from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param pm         StokesVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream&              is_xml,
                          StokesVector&         sv,
                          bifstream*            pbifs,
                          const Verbosity&      verbosity)
{
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("StokesVector");

  try
  {
    Tensor4 d;
    xml_read_from_stream(is_xml, d, pbifs, verbosity);
    sv = StokesVector(d);
  }
  catch (const std::runtime_error &e)
  {
    ostringstream os;
    os << "Error reading StokesVector: "
    << "\n" << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/StokesVector");
}


//! Writes StokesVector to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param pm         StokesVector
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream&                    os_xml,
                         const StokesVector&         sv,
                         bofstream*                  pbofs,
                         const String&               name,
                         const Verbosity&            verbosity)
{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("StokesVector");
  if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, sv.GetData(), pbofs, "", verbosity);

  close_tag.set_name("/StokesVector");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}


//=== TelsemAtlas ======================================================

//! Reads TelsemAtlas from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param pm         TelsemAtlas return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream&              is_xml,
                          TelsemAtlas&          ta,
                          bifstream*            pbifs,
                          const Verbosity&      verbosity)
{
    ArtsXMLTag tag(verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("TelsemAtlas");

    xml_read_from_stream(is_xml, ta.ndat, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.nchan, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.name, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.month, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.dlat, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.emis, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.correl, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.emis_err, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.classes1, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.classes2, pbifs, verbosity);
    xml_read_from_stream(is_xml, ta.cellnums, pbifs, verbosity);
    ta.telsem_calc_correspondence();
    tag.read_from_stream(is_xml);
    tag.check_name("/TelsemAtlas");
}


//! Writes TelsemAtlas to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param pm         TelsemAtlas
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream&                    os_xml,
                         const TelsemAtlas&          ta,
                         bofstream*                  pbofs,
                         const String&               name,
                         const Verbosity&            verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("TelsemAtlas");
    if (name.length())
        open_tag.add_attribute("name", name);

    open_tag.write_to_stream(os_xml);
    os_xml << '\n';
    xml_write_to_stream(os_xml, ta.ndat, pbofs, "ndat", verbosity);
    xml_write_to_stream(os_xml, ta.nchan, pbofs, "nchan", verbosity);
    xml_write_to_stream(os_xml, ta.name, pbofs, "name", verbosity);
    xml_write_to_stream(os_xml, ta.month, pbofs, "month", verbosity);
    xml_write_to_stream(os_xml, ta.dlat, pbofs, "dlat", verbosity);
    xml_write_to_stream(os_xml, ta.emis, pbofs, "emis", verbosity);
    xml_write_to_stream(os_xml, ta.correl, pbofs, "correl", verbosity);
    xml_write_to_stream(os_xml, ta.emis_err, pbofs, "emis_err", verbosity);
    xml_write_to_stream(os_xml, ta.classes1, pbofs, "class1", verbosity);
    xml_write_to_stream(os_xml, ta.classes2, pbofs, "class2", verbosity);
    xml_write_to_stream(os_xml, ta.cellnums, pbofs, "cellnum", verbosity);
    close_tag.set_name("/TelsemAtlas");
    close_tag.write_to_stream(os_xml);

    os_xml << '\n';
}


//=== XsecRecord ======================================================

//! Reads XsecData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param xd         XsecData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          XsecRecord& xd,
                          bifstream *pbifs,
                          const Verbosity& verbosity)
{
    CREATE_OUT2;
    ArtsXMLTag tag(verbosity);

    tag.read_from_stream(is_xml);
    tag.check_name("XsecRecord");

    String species_name;
    xml_read_from_stream(is_xml, species_name, pbifs, verbosity);

    const Index species = species_index_from_species_name(species_name);
    if (species == -1)
    {
        ostringstream os;
        os << "  Unknown species in XsecRecord: " << species_name;
        throw std::runtime_error(os.str());
    }
    xd.mspecies = species;

    Vector fmin, fmax;
    xml_read_from_stream(is_xml, xd.mcoeffs, pbifs, verbosity);
    xml_read_from_stream(is_xml, fmin, pbifs, verbosity);
    xml_read_from_stream(is_xml, fmax, pbifs, verbosity);
    xml_read_from_stream(is_xml, xd.mrefpressure, pbifs, verbosity);
    xml_read_from_stream(is_xml, xd.mreftemperature, pbifs, verbosity);
    xml_read_from_stream(is_xml, xd.mxsecs, pbifs, verbosity);

    xml_read_from_stream(is_xml, xd.mtslope, pbifs, verbosity);
    xml_read_from_stream(is_xml, xd.mtintersect, pbifs, verbosity);
    if (xd.mtslope.nelem() != xd.mxsecs.nelem()
        || xd.mtintersect.nelem() != xd.mxsecs.nelem())
    {
        std::ostringstream os;
        os << "  Bugged input data for " << species_name << ". "
           << "Number of bands mismatch.\n"
           << "Cross sections:      " << xd.mxsecs.nelem() << "\n"
           << "Tfit slope data:     " << xd.mtslope.nelem() << "\n"
           << "Tfit intersect data: " << xd.mtintersect.nelem() << "\n";
        throw std::runtime_error(os.str());
    }
    for (Index i=0; i < xd.mxsecs.nelem(); i++)
    {
        if (xd.mtslope[i].nelem() <= 1)
        {
            out2 << "  Warning: No temperature fit available for "
                 << species_name << " band " << i << "\n";
        }
        else if (xd.mtslope[i].nelem() != xd.mxsecs[i].nelem())
        {
            std::ostringstream os;
            os << "  Bugged input data for " << species_name << ". "
               << "Length of cross sections (" << xd.mxsecs[i].nelem()
               << ") does not match length of temperature fit data ("
               << xd.mtslope[i].nelem() << ").";
            throw std::runtime_error(os.str());
        }

        if (xd.mtslope[i].nelem() != xd.mtintersect[i].nelem())
        {
            std::ostringstream os;
            os << "  Bugged input data for " << species_name << ". "
               << "Length of temperature fit slope data (" << xd.mtslope[i].nelem()
               << ") does not match length of temperature fit intersect data ("
               << xd.mtintersect[i].nelem() << ") for band " << i << ".";
            throw std::runtime_error(os.str());
        }
    }
    tag.read_from_stream(is_xml);

    const Index ndatasets = xd.mxsecs.nelem();
    xd.mfgrids.resize(ndatasets);
    for (Index i = 0; i < ndatasets; i++)
    {
        nlinspace(xd.mfgrids[i], fmin[i], fmax[i], xd.mxsecs[i].nelem());
    }
    tag.check_name("/XsecRecord");
}


//! Writes XsecData to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param xd         XsecData
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const XsecRecord& xd,
                         bofstream *pbofs,
                         const String& name,
                         const Verbosity& verbosity)
{
    ArtsXMLTag open_tag(verbosity);
    ArtsXMLTag close_tag(verbosity);

    open_tag.set_name("XsecRecord");
    if (name.length())
        open_tag.add_attribute("name", name);

    open_tag.write_to_stream(os_xml);
    os_xml << '\n';
    const Index ndatasets = xd.Coeffs().nelem();
    Vector fmin(ndatasets), fmax(ndatasets);
    for (Index i = 0; i < ndatasets; i++)
    {
        fmin[i] = xd.Fgrids()[i][0];
        fmax[i] = xd.Fgrids()[i][xd.Fgrids()[i].nelem() - 1];
    }

    xml_write_to_stream(os_xml, xd.SpeciesName(), pbofs, "species", verbosity);
    xml_write_to_stream(os_xml, xd.Coeffs(), pbofs, "coefs", verbosity);
    xml_write_to_stream(os_xml, fmin, pbofs, "fmin", verbosity);
    xml_write_to_stream(os_xml, fmax, pbofs, "fmax", verbosity);
    xml_write_to_stream(os_xml, xd.RefPressure(), pbofs, "refpressure", verbosity);
    xml_write_to_stream(os_xml, xd.RefTemperature(), pbofs, "reftemperature", verbosity);
    xml_write_to_stream(os_xml, xd.Xsecs(), pbofs, "xsec", verbosity);
    xml_write_to_stream(os_xml, xd.TemperatureSlope(), pbofs, "xsec", verbosity);
    xml_write_to_stream(os_xml, xd.TemperatureIntersect(), pbofs, "xsec", verbosity);
    close_tag.set_name("/XsecRecord");
    close_tag.write_to_stream(os_xml);

    os_xml << '\n';
}


////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// FIXME: These should be implemented, sooner or later...

void xml_read_from_stream(istream&,
                          Agenda&,
                          bifstream* /* pbifs */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


void xml_write_to_stream(ostream&,
                         const Agenda&,
                         bofstream* /* pbofs */,
                         const String& /* name */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


//=== MCAntenna ================================================

void xml_read_from_stream(istream&,
                          MCAntenna&,
                          bifstream* /* pbifs */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


void xml_write_to_stream(ostream&,
                         const MCAntenna&,
                         bofstream* /* pbofs */,
                         const String& /* name */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


//=== TessemNN ================================================

void xml_read_from_stream(istream&,
                          TessemNN&,
                          bifstream* /* pbifs */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


void xml_write_to_stream(ostream&,
                         const TessemNN&,
                         bofstream* /* pbofs */,
                         const String& /* name */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


//=== Verbosity ================================================

void xml_read_from_stream(istream&,
                          Verbosity&,
                          bifstream* /* pbifs */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}


void xml_write_to_stream(ostream&,
                         const Verbosity&,
                         bofstream* /* pbofs */,
                         const String& /* name */, const Verbosity&)
{
  throw runtime_error("Method not implemented!");
}

