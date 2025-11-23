#include <debug.h>
#include <workspace.h>

void ReadCatalogData(PredefinedModelData& abs_predef_data,
                     XsecRecords& abs_xfit_data,
                     CIARecords& abs_cia_data,
                     AbsorptionBands& abs_bands,
                     const ArrayOfSpeciesTag& abs_species,
                     const String& basename,
                     const Index& ignore_missing) try {
  ARTS_TIME_REPORT

  abs_bandsReadSpeciesSplitCatalog(
      abs_bands, abs_species, basename + "lines/", ignore_missing);

  abs_cia_dataReadSpeciesSplitCatalog(
      abs_cia_data, abs_species, basename + "cia/", ignore_missing);

  abs_xfit_dataReadSpeciesSplitCatalog(
      abs_xfit_data, abs_species, basename + "xsec/", ignore_missing);

  abs_predef_dataReadSpeciesSplitCatalog(
      abs_predef_data, abs_species, basename + "predef/", 1, ignore_missing);
}
ARTS_METHOD_ERROR_CATCH