#include <wigner_functions.h>
#include <workspace.h>

/* Workspace method: Doxygen documentation will be auto-generated */
void WignerInit(const Index& fast_wigner_stored_symbols,
                const Index& largest_wigner_symbol_parameter,
                const Index& symbol_type) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(symbol_type != 3 and symbol_type != 6,
                     "Invalid symbol type: {} (must be 3 or 6)",
                     symbol_type);

  WignerInformation(static_cast<int>(largest_wigner_symbol_parameter),
                    static_cast<int>(fast_wigner_stored_symbols),
                    symbol_type >= 3,
                    symbol_type >= 6);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void WignerUnload() {
  ARTS_TIME_REPORT
  WignerInformation{}.unload();
}
