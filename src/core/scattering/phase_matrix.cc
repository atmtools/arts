#include "phase_matrix.h"

namespace scattering {

std::ostream &operator<<(std::ostream &out, Format format) {
  switch (format) {
    case Format::TRO:
      out << "TRO";
      break;
    case Format::ARO:
      out << "ARO";
      break;
    case Format::General:
      out << "General";
      break;
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, Representation repr) {
  switch (repr) {
    case Representation::Gridded:
      out << "gridded";
      break;
    case Representation::Spectral:
      out << "spectral";
      break;
    case Representation::DoublySpectral:
      out << "doubly-spectral";
      break;
  }
  return out;
}

}  // namespace scattering
