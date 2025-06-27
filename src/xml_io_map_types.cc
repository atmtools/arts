#include <xml_io_map_types.h>

//! Helper macro for Map types
#define TMPL_XML_READ_WRITE_STREAM_UMAP(Key, Type)                   \
  void xml_read_from_stream(std::istream& is_xml,                    \
                            std::unordered_map<Key, Type>& map,      \
                            bifstream* pbifs) {                      \
    xml_read(is_xml, map, pbifs);                                    \
  }                                                                  \
                                                                     \
  void xml_write_to_stream(std::ostream& os_xml,                     \
                           const std::unordered_map<Key, Type>& map, \
                           bofstream* pbofs,                         \
                           const String& name) {                     \
    xml_write(os_xml, map, pbofs, name);                             \
  }

TMPL_XML_READ_WRITE_STREAM_UMAP(QuantumIdentifier, AbsorptionBand)
TMPL_XML_READ_WRITE_STREAM_UMAP(QuantumIdentifier, GriddedField1)
TMPL_XML_READ_WRITE_STREAM_UMAP(QuantumIdentifier, Numeric)
TMPL_XML_READ_WRITE_STREAM_UMAP(QuantumIdentifier, Vector)

TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesEnum, AbsorptionLookupTable)
TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesEnum, Numeric)
TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesEnum, Vector)
TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesEnum, LinemixingSingleEcsData)

TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesIsotope, Numeric)
TMPL_XML_READ_WRITE_STREAM_UMAP(SpeciesIsotope, LinemixingSpeciesEcsData)

TMPL_XML_READ_WRITE_STREAM_UMAP(ScatteringSpeciesProperty, Numeric)
