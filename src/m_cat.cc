#include <workspace.h>

void AbsorptionReadSpeciesSplitCatalogs(
    PredefinedModelData& absorption_predefined_model_data,
    ArrayOfXsecRecord& absorption_xsec_fit_data,
    ArrayOfCIARecord& absorption_cia_data,
    AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const String& basename) {
  absorption_bandsReadSpeciesSplitCatalog(
      absorption_bands, absorption_species, basename + "lines/", 1);

  absorption_cia_dataReadSpeciesSplitCatalog(
      absorption_cia_data, absorption_species, basename + "cia/", 0);

  absorption_xsec_fit_dataReadSpeciesSplitCatalog(
      absorption_xsec_fit_data, absorption_species, basename + "xsec/");

  absorption_predefined_model_dataReadSpeciesSplitCatalog(
      absorption_predefined_model_data,
      absorption_species,
      basename + "predef/",
      1);
}