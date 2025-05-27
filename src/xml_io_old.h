#pragma once

#include <lbl_data.h>

struct ArtscatMeta {
  lbl::line data;
  QuantumIdentifier quantumidentity{};
  bool bad{true};
};

using ArrayOfArtscatMeta = std::vector<ArtscatMeta>;

void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfArtscatMeta& meta,
                          bifstream*);