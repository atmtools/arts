#include "scattering_species.h"

namespace scattering {

  std::ostream& operator<<(std::ostream& os,
                           const Species& /*species*/) {
    os << "A scattering species." << std::endl;
    return os;
}

}  // namespace Scattering
