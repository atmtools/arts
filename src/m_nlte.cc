#include <workspace.h>

void absorption_bandsSetNonLTE(AbsorptionBands& absorption_bands) {
  for (auto& [_, band] : absorption_bands) {
    band.lineshape = LineByLineLineshape::VP_LINE_NLTE;
  }
}

void atmospheric_fieldInitializeNonLTE(AtmField& atmospheric_field,
                                       const AbsorptionBands& absorption_bands,
                                       const Numeric& normalizing_factor) try {
  atmospheric_field.nlte() = lbl::nlte::from_lte(
      atmospheric_field, absorption_bands, normalizing_factor);
}
ARTS_METHOD_ERROR_CATCH
