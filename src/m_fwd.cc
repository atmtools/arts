#include <fwd.h>
#include <workspace.h>

#include <memory>

void spectral_radiance_operatorGeometricPlanar(
    const Workspace& ws,
    SpectralRadianceOperator& spectral_radiance_operator,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const AscendingGrid& altitude_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& cia_extrapolation,
    const Index& cia_robust) {
  ARTS_USER_ERROR_IF(altitude_grid.size() < 2, "Must have some type of path")

  using lines_t = AbsorptionBands;
  using cia_t = ArrayOfCIARecord;
  using xsec_t = ArrayOfXsecRecord;
  using predef_t = PredefinedModelData;

  const String lines_str = "absorption_bands";
  const String cia_str = "absorption_cia_data";
  const String xsec_str = "absorption_xsec_fit_data";
  const String predef_str = "absorption_predefined_model_data";

  auto lines = ws.wsv_and_contains(lines_str)
                   ? ws.share(lines_str)->share<lines_t>()
                   : std::shared_ptr<lines_t>{};
  auto cia = ws.wsv_and_contains(cia_str) ? ws.share(cia_str)->share<cia_t>()
                                          : std::shared_ptr<cia_t>{};
  auto xsec = ws.wsv_and_contains(xsec_str)
                  ? ws.share(xsec_str)->share<xsec_t>()
                  : std::shared_ptr<xsec_t>{};
  auto predef = ws.wsv_and_contains(predef_str)
                    ? ws.share(predef_str)->share<predef_t>()
                    : std::shared_ptr<predef_t>{};

  spectral_radiance_operator =
      fwd::spectral_radiance(altitude_grid,
                             latitude,
                             longitude,
                             atmospheric_field,
                             surface_field,
                             std::move(lines),
                             std::move(cia),
                             std::move(xsec),
                             std::move(predef),
                             cia_extrapolation,
                             cia_robust);
}
