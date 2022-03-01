#include <cstdlib>
#include <fstream>
#include <variant>

#include "global_data.h"

void define_wsv_group_names();

int main() {
  std::ofstream file("tokval.h");

  define_wsv_group_names();

  file << R"--(#ifndef auto_tokval_h
#define auto_tokval_h

#include "ppath_struct.h"
#include "absorptionlines.h"
#include "artstime.h"
#include "callback.h"
#include "cia.h"
#include "covariance_matrix.h"
#include "energylevelmap.h"
#include "gas_abs_lookup.h"
#include "hitran_xsec.h"
#include "linemixing.h"
#include "linemixing_hitran.h"
#include "matpackII.h"
#include "matpackVII.h"
#include "mc_antenna.h"
#include "optproperties.h"
#include "supergeneric.h"
#include "telsem.h"
#include "tessem.h"
#include "timer_struct.h"
#include "transmissionmatrix.h"

class TokVal {
  std::variant<
)--";

  bool first = true;
  for (auto& group : global_data::wsv_group_names) {
    if (group == "Agenda" or group == "ArrayOfAgenda") continue;
    if (not first) file << ",\n";
    first = false;
    file << "    std::unique_ptr<" << group << '>';
  }

  file << R"--(> value;

public:
)--";

  for (auto& group : global_data::wsv_group_names) {
    if (group == "Agenda" or group == "ArrayOfAgenda") continue;
    file << "  TokVal(" << group << " in) noexcept : value(std::make_unique<"
         << group << ">(std::move(in))) {}\n"
         << "  TokVal& operator=(" << group << " in) { value = std::make_unique<"
         << group << ">(std::move(in)); return *this; }\n"
         << "  operator " << group << R"--(() const { return *std::get<std::unique_ptr<)--"<< group << R"--(>>(value); }

)--";
  }

  file << R"--(  TokVal(const char * const c) : TokVal(String(c)) {}
    
  TokVal() noexcept : value(std::make_unique<Any>()) {} 
  TokVal(const TokVal& v) { std::visit([&](auto&& in) {*this = *in;}, v.value); }
  TokVal& operator=(const TokVal& v) { std::visit([&](auto&& in) {*this = *in;}, v.value); return *this; }
};

#endif
)--";

  return EXIT_SUCCESS;
}