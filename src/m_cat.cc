#include <debug.h>
#include <workspace.h>

void ReadCatalogData(PredefinedModelData& absorption_predefined_model_data,
                     ArrayOfXsecRecord& absorption_xsec_fit_data,
                     ArrayOfCIARecord& abs_cia_data,
                     AbsorptionBands& abs_bands,
                     const ArrayOfSpeciesTag& absorption_species,
                     const String& basename,
                     const Index& ignore_missing) try {
  ARTS_TIME_REPORT

  abs_bandsReadSpeciesSplitCatalog(abs_bands,
                                          absorption_species,
                                          basename + "lines/",
                                          ignore_missing);

  abs_cia_dataReadSpeciesSplitCatalog(abs_cia_data,
                                             absorption_species,
                                             basename + "cia/",
                                             ignore_missing);

  absorption_xsec_fit_dataReadSpeciesSplitCatalog(absorption_xsec_fit_data,
                                                  absorption_species,
                                                  basename + "xsec/",
                                                  ignore_missing);

  absorption_predefined_model_dataReadSpeciesSplitCatalog(
      absorption_predefined_model_data,
      absorption_species,
      basename + "predef/",
      1,
      ignore_missing);
}
ARTS_METHOD_ERROR_CATCH