#include <rtepack.h>

#include "xml_io.h"

#include "xml_io_array_macro.h"

//=== Propmat ================================================================

//! Reads Propmat from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Propmat return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          Propmat& pm,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;
  
  tag.read_from_stream(is_xml);
  tag.check_name("Propmat");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 7);
  else
    is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >> pm.W();

  tag.read_from_stream(is_xml);
  tag.check_name("/Propmat");
}

//! Writes Propmat to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Propmat
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const Propmat& pm,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Propmat");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  
  if (pbofs)
    *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V() << pm.W();
  else
    os_xml << ' ' << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D() << ' ' << pm.U() << ' ' << pm.V() << ' ' << pm.W() << ' ';

  close_tag.set_name("/Propmat");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Stokvec ================================================================

//! Reads Stokvec from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Stokvec return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          Stokvec& pm,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;
  
  tag.read_from_stream(is_xml);
  tag.check_name("Stokvec");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 4);
  else
    is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();

  tag.read_from_stream(is_xml);
  tag.check_name("/Stokvec");
}

//! Writes Stokvec to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Stokvec
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const Stokvec& pm,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Stokvec");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  
  if (pbofs)
    *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
  else
    os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' ' << pm.V() << ' ';

  close_tag.set_name("/Stokvec");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}


//=== Muelmat ================================================================

//! Reads Muelmat from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Muelmat return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          Muelmat& pm,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;
  
  tag.read_from_stream(is_xml);
  tag.check_name("Muelmat");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 16);
  else
    is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >> 
              pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >> 
              pm.data[8] >> pm.data[9] >> pm.data[10] >> pm.data[11] >> 
              pm.data[12] >> pm.data[13] >> pm.data[14] >> pm.data[15];

  tag.read_from_stream(is_xml);
  tag.check_name("/Muelmat");
}

//! Writes Muelmat to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Muelmat
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml, const Muelmat &pm,
                         bofstream *pbofs [[maybe_unused]], const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Muelmat");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  if (pbofs)
    *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3] << pm.data[4]
           << pm.data[5] << pm.data[6] << pm.data[7] << pm.data[8] << pm.data[9]
           << pm.data[10] << pm.data[11] << pm.data[12] << pm.data[13]
           << pm.data[14] << pm.data[15];
  else
    os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2] << ' '
           << pm.data[3] << '\n'
           << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6] << ' '
           << pm.data[7] << '\n'
           << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
           << ' ' << pm.data[11] << '\n'
           << ' ' << pm.data[12] << ' ' << pm.data[13] << ' ' << pm.data[14]
           << ' ' << pm.data[15] << ' ';

  close_tag.set_name("/Muelmat");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== PropmatVector ================================================================

//! Reads PropmatVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      PropmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, PropmatVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("PropmatVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 7);
    else
      is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >>
          pm.W();
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/PropmatVector");
}

//! Writes PropmatVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      PropmatVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const PropmatVector& pmv,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("PropmatVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
    if (pbofs)
      *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V()
             << pm.W();
    else
      os_xml << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D() << ' '
             << pm.U() << ' ' << pm.V() << ' ' << pm.W() << '\n';
  }

  close_tag.set_name("/PropmatVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecVector ================================================================

//! Reads StokvecVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, StokvecVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 4);
    else
      is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecVector");
}

//! Writes StokvecVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const StokvecVector& pmv,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
  if (pbofs)
    *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
  else
    os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' ' << pm.V() << '\n';
  }

  close_tag.set_name("/StokvecVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== MuelmatVector ================================================================

//! Reads MuelmatVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      MuelmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, MuelmatVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("MuelmatVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 16);
    else
      is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >> 
                pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >> 
                pm.data[8] >> pm.data[9] >> pm.data[10] >> pm.data[11] >> 
                pm.data[12] >> pm.data[13] >> pm.data[14] >> pm.data[15];
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/MuelmatVector");
}

//! Writes MuelmatVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      MuelmatVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const MuelmatVector& pmv,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MuelmatVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
    if (pbofs)
      *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3] << pm.data[4]
            << pm.data[5] << pm.data[6] << pm.data[7] << pm.data[8] << pm.data[9]
            << pm.data[10] << pm.data[11] << pm.data[12] << pm.data[13]
            << pm.data[14] << pm.data[15];
    else
      os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2] << ' '
            << pm.data[3] << '\n'
            << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6] << ' '
            << pm.data[7] << '\n'
            << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
            << ' ' << pm.data[11] << '\n'
            << ' ' << pm.data[12] << ' ' << pm.data[13] << ' ' << pm.data[14]
            << ' ' << pm.data[15] << '\n';
  }

  close_tag.set_name("/MuelmatVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== PropmatMatrix ================================================================

//! Reads PropmatMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      PropmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, PropmatMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("PropmatMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 7);
      else
        is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >>
            pm.W();
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/PropmatMatrix");
}

//! Writes PropmatMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      PropmatMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const PropmatMatrix& pmm,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("PropmatMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmv : pmm) {
    for (const auto &pm : pmv) {
      if (pbofs)
        *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V()
               << pm.W();
      else
        os_xml << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D()
               << ' ' << pm.U() << ' ' << pm.V() << ' ' << pm.W() << '\n';
    }
  }

  close_tag.set_name("/PropmatMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecMatrix ================================================================

//! Reads StokvecMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, StokvecMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 4);
      else
        is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecMatrix");
}

//! Writes StokvecMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const StokvecMatrix& pmm,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmv : pmm) {
    for (const auto &pm : pmv) {
      if (pbofs)
        *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
      else
        os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
               << pm.V() << '\n';
    }
  }

  close_tag.set_name("/StokvecMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== MuelmatMatrix ================================================================

//! Reads MuelmatMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      MuelmatMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml, MuelmatMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("MuelmatMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 16);
      else
        is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >>
            pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >>
            pm.data[8] >> pm.data[9] >> pm.data[10] >> pm.data[11] >>
            pm.data[12] >> pm.data[13] >> pm.data[14] >> pm.data[15];
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/MuelmatMatrix");
}

//! Writes MuelmatMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      MuelmatMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const MuelmatMatrix& pmm,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MuelmatMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto& pmv : pmm) {
      for (const auto &pm : pmv) {
        if (pbofs)
          *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3]
                 << pm.data[4] << pm.data[5] << pm.data[6] << pm.data[7]
                 << pm.data[8] << pm.data[9] << pm.data[10] << pm.data[11]
                 << pm.data[12] << pm.data[13] << pm.data[14] << pm.data[15];
        else
          os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2]
                 << ' ' << pm.data[3] << '\n'
                 << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6]
                 << ' ' << pm.data[7] << '\n'
                 << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
                 << ' ' << pm.data[11] << '\n'
                 << ' ' << pm.data[12] << ' ' << pm.data[13] << ' '
                 << pm.data[14] << ' ' << pm.data[15] << '\n';
      }
    }

  close_tag.set_name("/MuelmatMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

// arrays
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfPropmatMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfStokvecMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMuelmatMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfPropmatMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfStokvecMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfMuelmatMatrix)
