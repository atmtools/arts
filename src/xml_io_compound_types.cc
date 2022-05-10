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
#include "cloudbox.h"
#include "debug.h"
#include "global_data.h"
#include "predefined/predef_data.h"
#include "xml_io.h"
#include <sstream>

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== CIARecord ================================================

//! Reads CIARecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param irecord  CIARecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          CIARecord& cr,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  String name;
  String molecule1;
  String molecule2;
  Species::Species species1;
  Species::Species species2;

  tag.read_from_stream(is_xml);
  tag.check_name("CIARecord");
  tag.get_attribute_value("molecule1", molecule1);
  tag.get_attribute_value("molecule2", molecule2);

  species1 = Species::fromShortName(molecule1);
  species2 = Species::fromShortName(molecule2);

  if (not good_enum(species1)) {
    ostringstream os;
    os << "Unknown species (1st molecule) in CIARecord: " << molecule1;
    throw runtime_error(os.str());
  }
  if (not good_enum(species2)) {
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
                         const String& name _U_,
                         const Verbosity& verbosity) {
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
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  String name, type;
  Index n_blocks, row_start, row_extent, column_start, column_extent, row_index,
      column_index, is_inverse;

  tag.read_from_stream(is_xml);
  tag.check_name("CovarianceMatrix");
  tag.get_attribute_value("n_blocks", n_blocks);

  covmat = CovarianceMatrix();
  for (Index i = 0; i < n_blocks; i++) {
    tag.read_from_stream(is_xml);
    tag.check_name("Block");

    tag.get_attribute_value("row_index", row_index);
    tag.get_attribute_value("column_index", column_index);
    tag.get_attribute_value("row_start", row_start);
    tag.get_attribute_value("row_extent", row_extent);
    tag.get_attribute_value("column_start", column_start);
    tag.get_attribute_value("column_extent", column_extent);
    tag.get_attribute_value("type", type);
    tag.get_attribute_value("is_inverse", is_inverse);

    Range row_range(row_start, row_extent);
    Range column_range(column_start, column_extent);
    if (type == "Matrix") {
      std::shared_ptr<Matrix> M =
          std::make_shared<Matrix>(row_extent, column_extent);
      xml_read_from_stream(is_xml, *M, pbifs, verbosity);
      if (!is_inverse) {
        covmat.correlations_.emplace_back(
            row_range,
            column_range,
            std::make_pair(row_index, column_index),
            M);
      } else {
        covmat.inverses_.emplace_back(row_range,
                                      column_range,
                                      std::make_pair(row_index, column_index),
                                      M);
      }
    } else if (type == "Sparse") {
      std::shared_ptr<Sparse> M =
          std::make_shared<Sparse>(row_extent, column_extent);
      xml_read_from_stream(is_xml, *M, pbifs, verbosity);
      if (!is_inverse) {
        covmat.correlations_.emplace_back(
            row_range,
            column_range,
            std::make_pair(row_index, column_index),
            M);
      } else {
        covmat.inverses_.emplace_back(row_range,
                                      column_range,
                                      std::make_pair(row_index, column_index),
                                      M);
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
                         const String& name _U_,
                         const Verbosity& verbosity) {
  ArtsXMLTag covmat_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  covmat_tag.set_name("CovarianceMatrix");
  covmat_tag.add_attribute(
      "n_blocks", Index(covmat.correlations_.size() + covmat.inverses_.size()));
  covmat_tag.write_to_stream(os_xml);
  os_xml << '\n';
  for (const Block& c : covmat.correlations_) {
    ArtsXMLTag block_tag(verbosity);
    block_tag.set_name("Block");

    Index i, j;
    std::tie(i, j) = c.get_indices();
    block_tag.add_attribute("row_index", i);
    block_tag.add_attribute("column_index", j);

    Range row_range = c.get_row_range();
    Range column_range = c.get_column_range();
    block_tag.add_attribute("row_start", row_range.get_start());
    block_tag.add_attribute("row_extent", row_range.get_extent());
    block_tag.add_attribute("column_start", column_range.get_start());
    block_tag.add_attribute("column_extent", column_range.get_extent());
    block_tag.add_attribute("is_inverse", Index(0));
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
  for (const Block& c : covmat.inverses_) {
    ArtsXMLTag block_tag(verbosity);
    block_tag.set_name("Block");

    Index i, j;
    std::tie(i, j) = c.get_indices();
    block_tag.add_attribute("row_index", i);
    block_tag.add_attribute("column_index", j);

    Range row_range = c.get_row_range();
    Range column_range = c.get_column_range();
    block_tag.add_attribute("row_start", row_range.get_start());
    block_tag.add_attribute("row_extent", row_range.get_extent());
    block_tag.add_attribute("column_start", column_range.get_start());
    block_tag.add_attribute("column_extent", column_range.get_extent());
    block_tag.add_attribute("is_inverse", Index(1));
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

//=== EnergyLevelMap ===========================================================

//! Reads EnergyLevelMap from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gal     EnergyLevelMap return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          EnergyLevelMap& elm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  
  tag.check_name("EnergyLevelMap");
  String type;
  tag.get_attribute_value("type", type);
  elm.type = toEnergyLevelMapTypeOrThrow(type);

  xml_read_from_stream(is_xml, elm.levels, pbifs, verbosity);
  xml_read_from_stream(is_xml, elm.value, pbifs, verbosity);
  xml_read_from_stream(is_xml, elm.vib_energy, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/EnergyLevelMap");
  
  elm.ThrowIfNotOK();
}

//! Writes EnergyLevelMap to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gal     EnergyLevelMap
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const EnergyLevelMap& elm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("EnergyLevelMap");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("type", toString(elm.type));
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, elm.levels, pbofs, "Energy Levels", verbosity);
  xml_write_to_stream(os_xml, elm.value, pbofs, "Level Data", verbosity);
  xml_write_to_stream(os_xml, elm.vib_energy, pbofs, "Level Energy", verbosity);

  close_tag.set_name("/EnergyLevelMap");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GasAbsLookup");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, gal.species, pbofs, "", verbosity);
  xml_write_to_stream(
      os_xml, gal.nonlinear_species, pbofs, "NonlinearSpecies", verbosity);
  xml_write_to_stream(os_xml, gal.f_grid, pbofs, "FrequencyGrid", verbosity);
  xml_write_to_stream(os_xml, gal.p_grid, pbofs, "PressureGrid", verbosity);
  xml_write_to_stream(
      os_xml, gal.vmrs_ref, pbofs, "ReferenceVmrProfiles", verbosity);
  xml_write_to_stream(
      os_xml, gal.t_ref, pbofs, "ReferenceTemperatureProfile", verbosity);
  xml_write_to_stream(
      os_xml, gal.t_pert, pbofs, "TemperaturePerturbations", verbosity);
  xml_write_to_stream(os_xml,
                      gal.nls_pert,
                      pbofs,
                      "NonlinearSpeciesVmrPerturbations",
                      verbosity);
  xml_write_to_stream(
      os_xml, gal.xsec, pbofs, "AbsorptionCrossSections", verbosity);

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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  for (Index i = 0; i < gfield.get_dim(); i++) {
    tag.read_from_stream(is_xml);
    if (tag.get_name() == "Vector") {
      String s;
      tag.get_attribute_value("name", s);
      if (s.length()) gfield.set_grid_name(i, s);

      Vector v;
      xml_parse_from_stream(is_xml, v, pbifs, tag, verbosity);
      gfield.set_grid(i, v);
      tag.read_from_stream(is_xml);
      tag.check_name("/Vector");
    } else if (tag.get_name() == "Array") {
      String s;
      tag.get_attribute_value("name", s);
      if (s.length()) gfield.set_grid_name(i, s);

      tag.get_attribute_value("type", s);
      if (s == "String") {
        ArrayOfString as;
        xml_parse_from_stream(is_xml, as, pbifs, tag, verbosity);
        gfield.set_grid(i, as);
        tag.read_from_stream(is_xml);
        tag.check_name("/Array");

      } else {
        xml_parse_error(
            "Grids must be of type *Vector* or *ArrayOfString*\n"
            "but *ArrayOf" +
            s + "* found.");
      }
    } else {
      ostringstream os;
      os << "Grids must be of type *Vector* or *ArrayOfString*\n"
         << "but tag <" + tag.get_name() + "> found.";
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
void xml_write_to_stream(ostream& os_xml,
                         const GriddedField& gfield,
                         bofstream* pbofs,
                         const String& /* name */,
                         const Verbosity& verbosity) {
  for (Index i = 0; i < gfield.get_dim(); i++) {
    switch (gfield.get_grid_type(i)) {
      case GRID_TYPE_NUMERIC:
        xml_write_to_stream(os_xml,
                            gfield.get_numeric_grid(i),
                            pbofs,
                            gfield.get_grid_name(i),
                            verbosity);
        break;
      case GRID_TYPE_STRING:
        xml_write_to_stream(os_xml,
                            gfield.get_string_grid(i),
                            pbofs,
                            gfield.get_grid_name(i),
                            verbosity);
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField1");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField2");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField3");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField4");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField5");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("GriddedField6");

  String s;
  tag.get_attribute_value("name", s);
  if (s.length()) gfield.set_name(s);

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
                         const String& name,
                         const Verbosity& verbosity) {
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("GridPos");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml,
                      gpos.idx,
                      pbofs,
                      "OriginalGridIndexBelowInterpolationPoint",
                      verbosity);
  xml_write_to_stream(
      os_xml, gpos.fd[0], pbofs, "FractionalDistanceToNextPoint_1", verbosity);
  xml_write_to_stream(
      os_xml, gpos.fd[1], pbofs, "FractionalDistanceToNextPoint_2", verbosity);

  close_tag.set_name("/GridPos");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== HitranRelaxationMatrixData ================================================

//! Reads HitranRelaxationMatrixData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param hitran   HitranRelaxationMatrixData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          HitranRelaxationMatrixData& hitran,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("HitranRelaxationMatrixData");

  xml_read_from_stream(is_xml, hitran.W0pp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0pp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0rp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0rp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0qp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0qp, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0pr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0pr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0rr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0rr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0qr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0qr, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0pq, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0pq, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0rq, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0rq, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.W0qq, pbifs, verbosity);
  xml_read_from_stream(is_xml, hitran.B0qq, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/HitranRelaxationMatrixData");
}

//! Writes HitranRelaxationMatrixData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param hitran   HitranRelaxationMatrixData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const HitranRelaxationMatrixData& hitran,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("HitranRelaxationMatrixData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  
  xml_write_to_stream(os_xml, hitran.W0pp, pbofs, "W0pp", verbosity);
  xml_write_to_stream(os_xml, hitran.B0pp, pbofs, "B0pp", verbosity);
  xml_write_to_stream(os_xml, hitran.W0rp, pbofs, "W0rp", verbosity);
  xml_write_to_stream(os_xml, hitran.B0rp, pbofs, "B0rp", verbosity);
  xml_write_to_stream(os_xml, hitran.W0qp, pbofs, "W0qp", verbosity);
  xml_write_to_stream(os_xml, hitran.B0qp, pbofs, "B0qp", verbosity);
  xml_write_to_stream(os_xml, hitran.W0pr, pbofs, "W0pr", verbosity);
  xml_write_to_stream(os_xml, hitran.B0pr, pbofs, "B0pr", verbosity);
  xml_write_to_stream(os_xml, hitran.W0rr, pbofs, "W0rr", verbosity);
  xml_write_to_stream(os_xml, hitran.B0rr, pbofs, "B0rr", verbosity);
  xml_write_to_stream(os_xml, hitran.W0qr, pbofs, "W0qr", verbosity);
  xml_write_to_stream(os_xml, hitran.B0qr, pbofs, "B0qr", verbosity);
  xml_write_to_stream(os_xml, hitran.W0pq, pbofs, "W0pq", verbosity);
  xml_write_to_stream(os_xml, hitran.B0pq, pbofs, "B0pq", verbosity);
  xml_write_to_stream(os_xml, hitran.W0rq, pbofs, "W0rq", verbosity);
  xml_write_to_stream(os_xml, hitran.B0rq, pbofs, "B0rq", verbosity);
  xml_write_to_stream(os_xml, hitran.W0qq, pbofs, "W0qq", verbosity);
  xml_write_to_stream(os_xml, hitran.B0qq, pbofs, "B0qq", verbosity);

  close_tag.set_name("/HitranRelaxationMatrixData");
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Ppath");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(
      os_xml, ppath.dim, pbofs, "AtmosphericDimensionality", verbosity);
  xml_write_to_stream(
      os_xml, ppath.np, pbofs, "NumberOfPositionInPropagationPath", verbosity);
  xml_write_to_stream(
      os_xml, ppath.constant, pbofs, "PropagationPathConstant", verbosity);
  xml_write_to_stream(
      os_xml, ppath.background, pbofs, "RadiativeBackground", verbosity);
  xml_write_to_stream(os_xml,
                      ppath.start_pos,
                      pbofs,
                      "StartPositionOfPropagationPath",
                      verbosity);
  xml_write_to_stream(
      os_xml, ppath.start_los, pbofs, "StartLOSOfPropagationPath", verbosity);
  xml_write_to_stream(os_xml,
                      ppath.start_lstep,
                      pbofs,
                      "StartLstepOfPropagationPath",
                      verbosity);
  xml_write_to_stream(
      os_xml, ppath.pos, pbofs, "PropagationPathPointPositions", verbosity);
  xml_write_to_stream(os_xml, ppath.los, pbofs, "LineOfSight", verbosity);
  xml_write_to_stream(
      os_xml, ppath.r, pbofs, "PropagationPathPointRadii", verbosity);
  xml_write_to_stream(
      os_xml, ppath.lstep, pbofs, "PropagationPathPositionLength", verbosity);
  xml_write_to_stream(
      os_xml, ppath.end_pos, pbofs, "EndPositionOfPropagationPath", verbosity);
  xml_write_to_stream(
      os_xml, ppath.end_los, pbofs, "EndLOSOfPropagationPath", verbosity);
  xml_write_to_stream(
      os_xml, ppath.end_lstep, pbofs, "EndLstepPropagationPath", verbosity);
  xml_write_to_stream(
      os_xml, ppath.nreal, pbofs, "RefractiveIndexRealPart", verbosity);
  xml_write_to_stream(
      os_xml, ppath.ngroup, pbofs, "GroupRefractiveIndex", verbosity);
  xml_write_to_stream(
      os_xml, ppath.gp_p, pbofs, "PressureGridIndexPosition", verbosity);
  xml_write_to_stream(
      os_xml, ppath.gp_lat, pbofs, "LatitudeGridIndexPosition", verbosity);
  xml_write_to_stream(
      os_xml, ppath.gp_lon, pbofs, "LongitudeGridIndexPosition", verbosity);

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
void xml_read_from_stream(istream& is_xml,
                          PropagationMatrix& pm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("PropagationMatrix");

  try {
    Tensor4 d;
    xml_read_from_stream(is_xml, d, pbifs, verbosity);
    Index naa = d.nbooks();
    Index nza = d.npages();
    Index nf = d.nrows();
    Index nstokes_needed = d.ncols();
    pm = PropagationMatrix(nf, need2stokes<true>(nstokes_needed), nza, naa);
    pm.Data() = std::move(d); // destructive takeover
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading PropagationMatrix: "
       << "\n"
       << e.what();
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
void xml_write_to_stream(ostream& os_xml,
                         const PropagationMatrix& pm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("PropagationMatrix");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, pm.Data(), pbofs, "", verbosity);

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
                          const Verbosity& verbosity) {
  static_assert(QuantumIdentifier::version == 1);

  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("QuantumIdentifier");
  
  Index version;
  if (tag.has_attribute("version")) {
    tag.get_attribute_value("version", version);
  } else {
    version = 0;
  }
  ARTS_USER_ERROR_IF(
      QuantumIdentifier::version < version,
      "The version of this quantum identifier is too new.  You need to upgrade ARTS to use it.")

  try {
    String qi_str;
    parse_xml_tag_content_as_string(is_xml, qi_str);
    qi = QuantumIdentifier(qi_str, version);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading QuantumIdentifier: "
       << "\n"
       << e.what();
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  static_assert(QuantumIdentifier::version == 1);

  open_tag.set_name("QuantumIdentifier");
  open_tag.add_attribute("version", QuantumIdentifier::version);
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  os_xml << qi;

  close_tag.set_name("/QuantumIdentifier");
  close_tag.write_to_stream(os_xml);
  os_xml << endl;
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Jacobian::Target target;
  String subtag;
  String subsubtag;
  String mode;
  ArrayOfVector grids;

  tag.read_from_stream(is_xml);
  tag.check_name("RetrievalQuantity");

  xml_read_from_stream(is_xml, target, pbifs, verbosity);
  xml_read_from_stream(is_xml, subtag, pbifs, verbosity);
  xml_read_from_stream(is_xml, subsubtag, pbifs, verbosity);
  xml_read_from_stream(is_xml, mode, pbifs, verbosity);
  xml_read_from_stream(is_xml, grids, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/RetrievalQuantity");

  rq = RetrievalQuantity(
    target, subtag, subsubtag, mode, target.perturbation, grids);
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("RetrievalQuantity");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, rq.Target(), pbofs, "", verbosity);
  xml_write_to_stream(os_xml, rq.Subtag(), pbofs, "Subtag", verbosity);
  xml_write_to_stream(os_xml, rq.SubSubtag(), pbofs, "SubSubtag", verbosity);
  xml_write_to_stream(os_xml, rq.Mode(), pbofs, "Mode", verbosity);
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("SingleScatteringData");
  tag.get_attribute_value("version", version);

  if (version == "3") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs, verbosity);
    ssdata.ptype = PTypeFromString(ptype_string);
  } else if (version == "2") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs, verbosity);
    ssdata.ptype = PType2FromString(ptype_string);
  } else {
    Index ptype;
    xml_read_from_stream(is_xml, ptype, pbifs, verbosity);
    if (ptype != PTYPE_GENERAL && ptype != PTYPE_TOTAL_RND &&
        ptype != PTYPE_AZIMUTH_RND) {
      ostringstream os;
      os << "Ptype value (" << ptype << ") is wrong."
         << "It must be \n"
         << PTYPE_TOTAL_RND << " - totally randomly oriented particles,\n"
         << PTYPE_AZIMUTH_RND
         << " - azimuthally randomly oriented particles, or\n"
         << PTYPE_GENERAL << " - arbitrary oriented particles.\n";
      throw runtime_error(os.str());
    }
    ssdata.ptype = PType(ptype);
  }
  xml_read_from_stream(is_xml, ssdata.description, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.f_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.T_grid, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.za_grid, pbifs, verbosity);
  /* Verify that we have a good coverage for the za grid */
  if ((ssdata.za_grid[0] > 1) ||
      ssdata.za_grid[ssdata.za_grid.nelem() - 1] < 179) {
    ostringstream os;
    os << "Missing data in xml-stream. Expected za_grid: [0, 180]. "
       << "Found za_grid: [" << ssdata.za_grid[0] << ", "
       << ssdata.za_grid[ssdata.za_grid.nelem() - 1] << "]";
    throw runtime_error(os.str());
  }
  xml_read_from_stream(is_xml, ssdata.aa_grid, pbifs, verbosity);

  xml_read_from_stream(is_xml, ssdata.pha_mat_data, pbifs, verbosity);
  if (ssdata.pha_mat_data.nlibraries() != ssdata.f_grid.nelem()) {
    throw runtime_error(
        "Number of frequencies in f_grid and pha_mat_data "
        "not matching!!!");
  }

  xml_read_from_stream(is_xml, ssdata.ext_mat_data, pbifs, verbosity);
  xml_read_from_stream(is_xml, ssdata.abs_vec_data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/SingleScatteringData");

  if (version != "3" && ssdata.ptype == PTYPE_AZIMUTH_RND) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SingleScatteringData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  os_xml << '\n';
  xml_write_to_stream(
      os_xml, PTypeToString(ssdata.ptype), pbofs, "", verbosity);
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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("ScatteringMetaData");
  tag.get_attribute_value("version", version);

  if (version != "3") {
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
  xml_read_from_stream(
      is_xml, smdata.diameter_area_equ_aerodynamical, pbifs, verbosity);

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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("ScatteringMetaData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, smdata.description, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.source, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.refr_index, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.mass, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.diameter_max, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, smdata.diameter_volume_equ, pbofs, "", verbosity);
  xml_write_to_stream(
      os_xml, smdata.diameter_area_equ_aerodynamical, pbofs, "", verbosity);

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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SLIData2");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, slidata.x1a, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, slidata.x2a, pbofs, "", verbosity);
  xml_write_to_stream(os_xml, slidata.ya, pbofs, "", verbosity);

  close_tag.set_name("/SLIData2");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SpeciesIsotopologueRatios ===========================================

//! Reads SpeciesIsotopologueRatios from XML input stream
/*!
  \param is_xml   XML Input stream
  \param sap      SpeciesIsotopologueRatios return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          SpeciesIsotopologueRatios& iso_rat,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(pbifs, "No support for binary IO for (SpeciesIsotopologueRatios)")
  
  CREATE_OUT2;
  
  iso_rat = SpeciesIsotopologueRatios{};

  ArtsXMLTag tag(verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesIsotopologueRatios");
  
  Index nelem;
  tag.get_attribute_value("nelem", nelem);
  
  String name;
  Numeric val;
  for (Index n = 0; n < nelem; n++) {
    is_xml >> name >> double_imanip() >> val;
    
    if (is_xml.fail()) {
      ostringstream os;
      os << " near "
      << "\n  Element: " << n;
      xml_data_parse_error(tag, os.str());
    }
    
    const Index i = Species::find_species_index(Species::Tag(name).Isotopologue());
    ARTS_USER_ERROR_IF(i < 0 or i >= iso_rat.maxsize,
                       "Species: ", name, " cannot be understood as a species by your compiled version of ARTS")
    
    iso_rat.data[i] = val;
  }
  
  tag.read_from_stream(is_xml);
  tag.check_name("/SpeciesIsotopologueRatios");
}

//! Writes SpeciesIsotopologueRatios to XML output stream
/*!
  \param os_xml   XML Output stream
  \param sap      SpeciesIsotopologueRatios
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const SpeciesIsotopologueRatios& iso_rat,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity)

{
  ARTS_USER_ERROR_IF(pbofs, "No support for binary IO for (SpeciesIsotopologueRatios)")
  
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);
  
  open_tag.set_name("SpeciesIsotopologueRatios");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", iso_rat.maxsize);
  
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  
  xml_set_stream_precision(os_xml);
  os_xml << iso_rat << '\n';
  
  close_tag.set_name("/SpeciesIsotopologueRatios");
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
                          bifstream* /* pbifs */,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  stringbuf strbuf;
  char dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesTag");

  // Skip whitespaces
  bool string_starts_with_quotes = true;
  do {
    is_xml >> dummy;
    switch (dummy) {
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
  if (!string_starts_with_quotes) {
    xml_parse_error("SpeciesTag must begin with \"");
  }

  is_xml.get(strbuf, '"');
  if (is_xml.fail()) {
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
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("SpeciesTag");
  if (name.length()) open_tag.add_attribute("name", name);
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
 * \param sv         StokesVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          StokesVector& sv,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("StokesVector");

  try {
    Tensor4 d;
    xml_read_from_stream(is_xml, d, pbifs, verbosity);
    Index naa = d.nbooks();
    Index nza = d.npages();
    Index nf = d.nrows();
    Index nstokes_needed = d.ncols();
    sv = StokesVector(nf, need2stokes<false>(nstokes_needed), nza, naa);
    sv.Data() = std::move(d); // destructive takeover
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading StokesVector: "
       << "\n"
       << e.what();
    throw runtime_error(os.str());
  }
  
  tag.read_from_stream(is_xml);
  tag.check_name("/StokesVector");
}

//! Writes StokesVector to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param sv         StokesVector
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const StokesVector& sv,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("StokesVector");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, sv.Data(), pbofs, "", verbosity);

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
void xml_read_from_stream(istream& is_xml,
                          TelsemAtlas& ta,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
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
void xml_write_to_stream(ostream& os_xml,
                         const TelsemAtlas& ta,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("TelsemAtlas");
  if (name.length()) open_tag.add_attribute("name", name);

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
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  CREATE_OUT2;
  ArtsXMLTag tag(verbosity);
  Index version;

  tag.read_from_stream(is_xml);
  tag.check_name("XsecRecord");
  tag.get_attribute_value("version", version);
  xd.SetVersion(version);

  String species_name;
  xml_read_from_stream(is_xml, species_name, pbifs, verbosity);

  const Species::Species species = Species::fromShortName(species_name);
  if (not good_enum(species)) {
    ostringstream os;
    os << "  Unknown species in XsecRecord: " << species_name;
    throw std::runtime_error(os.str());
  }
  xd.SetSpecies(species);

  ARTS_USER_ERROR_IF(version != 2, "Only XsecRecord version 2 is supported")

  xml_read_from_stream(is_xml, xd.FitMinPressures(), pbifs, verbosity);
  xml_read_from_stream(is_xml, xd.FitMaxPressures(), pbifs, verbosity);
  xml_read_from_stream(is_xml, xd.FitMinTemperatures(), pbifs, verbosity);
  xml_read_from_stream(is_xml, xd.FitMaxTemperatures(), pbifs, verbosity);
  xml_read_from_stream(is_xml, xd.FitCoeffs(), pbifs, verbosity);

  for (const auto& fitcoeffs : xd.FitCoeffs()) {
    const Index ncoeff = fitcoeffs.data.ncols();
    ARTS_USER_ERROR_IF(
        ncoeff != 4, "Wrong number of coefficients, expected 4, found ", ncoeff)
  }

  tag.read_from_stream(is_xml);
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
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("XsecRecord");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", xd.Version());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  xml_write_to_stream(os_xml, xd.SpeciesName(), pbofs, "species", verbosity);

  xml_write_to_stream(os_xml,
                      xd.FitMinPressures(),
                      pbofs,
                      "Mininum pressures from fit",
                      verbosity);
  xml_write_to_stream(os_xml,
                      xd.FitMaxPressures(),
                      pbofs,
                      "Maximum pressures from fit",
                      verbosity);
  xml_write_to_stream(os_xml,
                      xd.FitMinTemperatures(),
                      pbofs,
                      "Mininum temperatures from fit",
                      verbosity);
  xml_write_to_stream(os_xml,
                      xd.FitMaxTemperatures(),
                      pbofs,
                      "Maximum temperatures from fit",
                      verbosity);
  xml_write_to_stream(
      os_xml, xd.FitCoeffs(), pbofs, "Fit coefficients", verbosity);

  close_tag.set_name("/XsecRecord");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== MapOfErrorCorrectedSuddenData ======================================================

//! Reads MapOfErrorCorrectedSuddenData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param rvb        MapOfErrorCorrectedSuddenData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          MapOfErrorCorrectedSuddenData& rvb,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(pbifs not_eq nullptr, "No binary data")
  
  CREATE_OUT2;
  ArtsXMLTag open_tag(verbosity);
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("MapOfErrorCorrectedSuddenData");
  
  Index nelem;
  open_tag.get_attribute_value("nelem", nelem);
  
  for (Index i=0; i<nelem; i++) {
    ArtsXMLTag internal_open_tag(verbosity);
    internal_open_tag.read_from_stream(is_xml);
    internal_open_tag.check_name("ErrorCorrectedSuddenData");
    
    // Get key
    String val;
    internal_open_tag.get_attribute_value("key", val);
    auto& data = rvb[QuantumIdentifier(val)];
    
    // Get size
    Index nelem_specs;
    internal_open_tag.get_attribute_value("nelem", nelem_specs);
    
    // Get values
    for (Index j=0; j<nelem_specs; j++) {
      SpeciesErrorCorrectedSuddenData secds;
      is_xml >> secds;
      data[secds.spec] = secds;
    }
    
    ArtsXMLTag internal_close_tag(verbosity);
    internal_close_tag.read_from_stream(is_xml);
    internal_close_tag.check_name("/ErrorCorrectedSuddenData");
    
  }
  
  ArtsXMLTag close_tag(verbosity);
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/MapOfErrorCorrectedSuddenData");
  
  // Sanity check, it is not OK to not have AIR as catch-all broadener
  for (auto& x: rvb) {
    bool found_air=false;
    for (auto& y: x.data) {
      found_air = found_air or (y.spec == Species::Species::Bath);
    }
    ARTS_USER_ERROR_IF(not found_air,
      "Incomplete ErrorCorrectedSuddenData, must contain air, contains:\n",
      x)
  }
}

//! Writes MapOfErrorCorrectedSuddenData to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param rvb        MapOfErrorCorrectedSuddenData
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const MapOfErrorCorrectedSuddenData& rvb,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(pbofs not_eq nullptr, "No binary data")
  
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("MapOfErrorCorrectedSuddenData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", rvb.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  
  xml_set_stream_precision(os_xml);
  
  for (auto& r: rvb) {
    ArtsXMLTag internal_open_tag(verbosity);
    internal_open_tag.set_name("ErrorCorrectedSuddenData");
    internal_open_tag.add_attribute("key", var_string(r.id));
    internal_open_tag.add_attribute("nelem", r.data.nelem());
    internal_open_tag.write_to_stream(os_xml);
    os_xml << '\n';
    
    // Set values
    os_xml << r << '\n';
    
    ArtsXMLTag internal_close_tag(verbosity);
    internal_close_tag.set_name("/ErrorCorrectedSuddenData");
    internal_close_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }
  
  close_tag.set_name("/MapOfErrorCorrectedSuddenData");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}



//=== PredefinedModelData =========================================
/*!
 * \param is_xml     XML Input stream
 * \param pmd        PredefinedModelData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          PredefinedModelData& pmd,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(pbifs, "No binary data")

  pmd = PredefinedModelData{};  // overwrite

  CREATE_OUT2;
  ArtsXMLTag open_tag(verbosity);
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("PredefinedModelData");

  Index nelem;
  open_tag.get_attribute_value("nelem", nelem);

  for (Index i = 0; i < nelem; i++) {
    ArtsXMLTag internal_open_tag(verbosity);
    internal_open_tag.read_from_stream(is_xml);
    internal_open_tag.check_name("Data");

    // Get key
    String key_str;
    internal_open_tag.get_attribute_value("key", key_str);
    auto key =
        Absorption::PredefinedModel::toPredefinedModelDataKeyOrThrow(key_str);

    String sizes_str;
    internal_open_tag.get_attribute_value("sizes", sizes_str);

    Index sizes_len;
    internal_open_tag.get_attribute_value("sizes_nelem", sizes_len);

    std::vector<std::size_t> sizes(sizes_len);
    std::istringstream values(sizes_str);
    for (auto& sz : sizes) values >> sz;

    pmd.resize(sizes, key);
    pmd.set_data_from_stream(is_xml, key);

    ArtsXMLTag internal_close_tag(verbosity);
    internal_close_tag.read_from_stream(is_xml);
    internal_close_tag.check_name("/Data");
  }

  ArtsXMLTag close_tag(verbosity);
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/PredefinedModelData");
}

/*!
 * \param os_xml     XML Output stream
 * \param rvb        PredefinedModelData
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const PredefinedModelData& pmd,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(pbofs, "No binary data")

  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("PredefinedModelData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", Index(pmd.size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  const auto keys = pmd.keys();

  for (auto& key : keys) {
    ArtsXMLTag internal_open_tag(verbosity);
    internal_open_tag.set_name("Data");
    internal_open_tag.add_attribute("key", toString(key));

    auto sizes = pmd.data_size(key);
    internal_open_tag.add_attribute("sizes_nelem", Index(sizes.size()));

    String sizes_str = "";
    for (std::size_t i = 0; i < sizes.size(); i++) {
      if (i > 0) sizes_str += " ";
      sizes_str += var_string(sizes[i]);
    }
    internal_open_tag.add_attribute("sizes", sizes_str);

    internal_open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    // Set values
    pmd.output_data_to_stream(os_xml, key);
    os_xml << '\n';

    ArtsXMLTag internal_close_tag(verbosity);
    internal_close_tag.set_name("/Data");
    internal_close_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }

  close_tag.set_name("/PredefinedModelData");
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
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const Agenda&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== MCAntenna ================================================

void xml_read_from_stream(istream&,
                          MCAntenna&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const MCAntenna&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== TessemNN ================================================

void xml_read_from_stream(istream&,
                          TessemNN&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const TessemNN&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== Verbosity ================================================

void xml_read_from_stream(istream&,
                          Verbosity&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const Verbosity&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== CallbackFunction =========================================

void xml_read_from_stream(istream&,
                          CallbackFunction&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const CallbackFunction&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR("Method not implemented!");
}
