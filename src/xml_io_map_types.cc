#include <xml_io_map_types.h>

//! Helper macro for Map types
#define TMPL_XML_READ_WRITE_STREAM_MAP(T)                                     \
  void xml_read_from_stream(std::istream& is_xml, T& map, bifstream* pbifs) { \
    xml_read(is_xml, map, pbifs);                                             \
  }                                                                           \
                                                                              \
  void xml_write_to_stream(std::ostream& os_xml,                              \
                           const T& map,                                      \
                           bofstream* pbofs,                                  \
                           const String& name) {                              \
    xml_write(os_xml, map, pbofs, name);                                      \
  }

TMPL_XML_READ_WRITE_STREAM_MAP(AbsorptionBands)
TMPL_XML_READ_WRITE_STREAM_MAP(AbsorptionLookupTables)
TMPL_XML_READ_WRITE_STREAM_MAP(SpeciesEnumVectors)
TMPL_XML_READ_WRITE_STREAM_MAP(AtmPoint::SpeciesMap)
TMPL_XML_READ_WRITE_STREAM_MAP(AtmPoint::SpeciesIsotopeMap)
TMPL_XML_READ_WRITE_STREAM_MAP(AtmPoint::NlteMap)
TMPL_XML_READ_WRITE_STREAM_MAP(AtmPoint::ScatteringSpeciesMap)
// TMPL_XML_READ_WRITE_STREAM_MAP(AtmField::map_type<AtmKey>)
// TMPL_XML_READ_WRITE_STREAM_MAP(AtmField::map_type<SpeciesEnum>)
// TMPL_XML_READ_WRITE_STREAM_MAP(AtmField::map_type<SpeciesIsotope>)
// TMPL_XML_READ_WRITE_STREAM_MAP(AtmField::map_type<QuantumIdentifier>)
// TMPL_XML_READ_WRITE_STREAM_MAP(AtmField::map_type<ScatteringSpeciesProperty>)
TMPL_XML_READ_WRITE_STREAM_MAP(NonlteLineFluxProfile)
