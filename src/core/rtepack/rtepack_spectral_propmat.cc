#include "rtepack_spectral_propmat.h"

namespace rtepack {
void SpectralPropmat::init(const Size nf, const Size nq) {
  K.resize(nf);
  dK.resize(nq, nf);
  nlte.resize(nf);
  dnlte.resize(nq, nf);

  K       = propmat{0.0};
  dK      = propmat{0.0};
  nlte    = stokvec{0.0};
  dnlte   = stokvec{0.0};
  nlte_on = false;
}

bool SpectralPropmat::ok() const {
  return K.size() == nlte.size() and dK.ncols() == dnlte.ncols() and dK.nrows() == dnlte.nrows() and
         static_cast<Size>(dK.ncols()) == K.size();
}

std::array<Index, 2> SpectralPropmat::shape() const noexcept { return {static_cast<Index>(K.size()), dK.nrows()}; }
}  // namespace rtepack
