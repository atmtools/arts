#include "scattering_species.h"
#include "m_append.h"
#include "m_select.h"


std::ostream & operator<<(std::ostream &out, const ScatteringSpecies &s) {
    out << "A scattering species.";
    return out;
}


void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ArrayOfScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity) {
    Append(reinterpret_cast<Array<ScatteringSpecies>&>(out),
           out_name,
           reinterpret_cast<const Array<ScatteringSpecies>&>(in),
           direction,
           in_name,
           direction_name,
           verbosity);
}

void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity) {
    Append(reinterpret_cast<Array<ScatteringSpecies>&>(out),
           out_name,
           in,
           direction,
           in_name,
           direction_name,
           verbosity);
}

void Select(  // WS Generic Output:
    ArrayOfScatteringSpecies& needles,
    // WS Generic Input:
    const ArrayOfScatteringSpecies& haystack,
    const ArrayOfIndex& needleind,
    const Verbosity& verbosity) {
    Select(reinterpret_cast<Array<ScatteringSpecies> &>(needles),
           reinterpret_cast<const Array<ScatteringSpecies> &>(haystack),
           needleind,
           verbosity);
}
