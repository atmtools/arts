////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#include "xml_io.h"
#include "bifstream.h"
#include "bofstream.h"
#include "debug.h"
#include "file.h"
#include "parameters.h"


////////////////////////////////////////////////////////////////////////////
//   ArtsXMLTag implementation
////////////////////////////////////////////////////////////////////////////

void ArtsXMLTag::add_attribute(const String& aname, const std::vector<QuantumNumberType>& value) {
  std::ostringstream v;
  
  if(value.size() == 0)
    v << "";
  else {
    for(size_t i=0; i<value.size()-1; i++)
      v << value[i] << ' ';
    v << value.back();
  }
  
  add_attribute(aname, v.str());
}

void ArtsXMLTag::add_attribute(const String& aname, const ArrayOfSpecies& value, const bool self, const bool bath) {
  std::ostringstream v;
  
  if(self)
    v << LineShape::self_broadening;
  for(Index i=Index(self); i<value.nelem()-Index(bath); i++)
    v << ' ' << Species::toShortName(value[i]);
  if(bath) {
    v << ' ' << LineShape::bath_broadening;
  }
  
  add_attribute(aname, v.str());
}

void ArtsXMLTag::get_attribute_value(const String& aname, SpeciesTag& value) {
  String attribute_value;
  
  get_attribute_value(aname, attribute_value);
  value = SpeciesTag(attribute_value);
}

void ArtsXMLTag::get_attribute_value(const String& aname, ArrayOfSpecies& value, bool& self, bool& bath) {
  value.resize(0);
  self=false;
  bath=false;
  
  String attribute_value;
  std::istringstream strstr("");
  
  get_attribute_value(aname, attribute_value);
  if (attribute_value.nelem() == 0) return;
  
  strstr.str(attribute_value);
  String val;
  
  while(not strstr.eof()) {
    strstr >> val;
    if (strstr.fail()) {
      xml_parse_error("Error while parsing value of " + aname + " from <" + name +
      ">");
    }
    
    if(val == LineShape::self_broadening) {
      value.push_back(Species::Species::FINAL);
      self = true;
    }
    else if(val == LineShape::bath_broadening) {
      value.push_back(Species::Species::Bath);
      bath = true;
    }
    else {
      Species::Species x = Species::fromShortName(val);
      ARTS_USER_ERROR_IF(not good_enum(x), "Species: ", val, " cannot be understood")
      value.push_back(x);
    }
  }
}

void ArtsXMLTag::get_attribute_value(const String& aname, std::vector<QuantumNumberType>& value) {
  value.resize(0);
  
  String attribute_value;
  std::istringstream strstr("");
  
  get_attribute_value(aname, attribute_value);
  if (attribute_value.nelem() == 0) return;
  
  strstr.str(attribute_value);
  String val;
  
  while(not strstr.eof()) {
    strstr >> val;
    if (strstr.fail()) {
      xml_parse_error("Error while parsing value of " + aname + " from <" + name +
      ">");
    }
    value.push_back(Quantum::Number::toType(val));
  }
}

void xml_find_and_open_input_file(std::shared_ptr<std::istream>& ifs,
                                  const String& filename) {
  String xml_file = filename;
  find_xml_file(xml_file);

  // Open input stream:
  if (xml_file.substr(xml_file.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::shared_ptr<std::istream>(new igzstream());
    xml_open_input_file(
        *(static_cast<igzstream*>(ifs.get())), xml_file);
  }
#else
  {
    throw std::runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be read.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ifs = std::shared_ptr<std::istream>(new std::ifstream());
    xml_open_input_file(
        *(static_cast<std::ifstream*>(ifs.get())), xml_file);
  }
}

////////////////////////////////////////////////////////////////////////////
//   General XML functions (file header, start root tag, end root tag)
////////////////////////////////////////////////////////////////////////////

//! Throws XML parser runtime error
/*!
  This is used quite often inside the data parsing routines so it's a
  function for itself.

  \param tag        ArtsXMLTag
  \param str_error  Error description
*/
void xml_data_parse_error(ArtsXMLTag& tag, String str_error) {
  std::ostringstream os;
  os << "XML data parse error: Error reading ";
  tag.write_to_stream(os);
  os << str_error << "\n"
     << "Check syntax of XML file. A possible cause is that the file "
     << "contains NaN or Inf values.\n";
  throw std::runtime_error(os.str());
}

////////////////////////////////////////////////////////////////////////////
//   Default file name
////////////////////////////////////////////////////////////////////////////

//! Gives the default filename for the XML formats.
/*!
  The default name is only used if the filename is empty.

  \param filename filename
  \param varname variable name
*/
void filename_xml(const String& filename) {
  ARTS_USER_ERROR_IF(filename == "", "Must have filename")
}

//! Gives the default filename, with file index, for the XML formats.
/*!
  The default name is only used if the filename is empty.

  \param[out] filename   filename
  \param[in]  file_index Index appended to the filename
  \param[in]  digits     Width for padding with zeros
*/
void filename_xml_with_index(String& filename,
                             const Index& file_index,
                             const Index& digits) {
  ARTS_USER_ERROR_IF ("" == filename, "Must have filename")
  var_string(filename, ".", std::setw((int)digits), std::setfill('0'), file_index, ".xml");
}
