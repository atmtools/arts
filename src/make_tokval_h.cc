#include <cstdlib>
#include <fstream>
#include <variant>

#include "global_data.h"

void define_wsv_groups();

int main() {
  std::ofstream file_h("tokval.h");
  std::ofstream file_cc("tokval.cc");

  define_wsv_groups();

  file_h << R"--(// auto-generated tokval interface

#ifndef auto_tokval_h
#define auto_tokval_h

#include "ppath_struct.h"
#include "absorptionlines.h"
#include "artstime.h"
#include "callback.h"
#include "cia.h"
#include "covariance_matrix.h"
#include "energylevelmap.h"
#include "gas_abs_lookup.h"
#include "xsec_fit.h"
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

struct TokVal {
  std::variant<
)--";

  bool first = true;
  for (auto& group : global_data::wsv_groups) {
    if (group == "Agenda" or group == "ArrayOfAgenda") continue;
    if (not first) file_h << ",\n";
    first = false;
    file_h << "    std::unique_ptr<" << group << '>';
  }

  file_h << R"--(> value;

)--";

  for (auto& group : global_data::wsv_groups) {
    if (group == "Agenda" or group == "ArrayOfAgenda") continue;
    file_h << "  TokVal(" << group << " in) noexcept;\n"
         << "  TokVal& operator=(" << group << " in);\n"
         << "  operator " << group << R"--(() const;

)--";
  }

  file_h << R"--(  TokVal(const char * const c);
    
  TokVal() noexcept;
  TokVal(const TokVal& v);
  TokVal& operator=(const TokVal& v);

  friend std::ostream& operator<<(std::ostream& os, const TokVal& t);
};

#endif
)--";



file_cc << "// auto-generated tokval implementation\n\n#include <tokval.h>\n\n";

  for (auto& group : global_data::wsv_groups) {
    if (group == "Agenda" or group == "ArrayOfAgenda") continue;
    file_cc << "TokVal::TokVal(" << group << " in) noexcept : value(std::make_unique<"
         << group << ">(std::move(in))) {}\n"
         << "TokVal& TokVal::operator=(" << group << " in) { value = std::make_unique<"
         << group << ">(std::move(in)); return *this; }\n"
         << "TokVal::operator " << group << R"--(() const { return *std::get<std::unique_ptr<)--"<< group << R"--(>>(value); }

)--";
  }

  file_cc << R"--(TokVal::TokVal(const char * const c) : TokVal(String(c)) {}
    
TokVal::TokVal() noexcept : value(std::make_unique<Any>()) {} 
TokVal::TokVal(const TokVal& v) { std::visit([&](auto&& in) {*this = *in;}, v.value); }
TokVal& TokVal::operator=(const TokVal& v) { std::visit([&](auto&& in) {*this = *in;}, v.value); return *this; }

std::ostream& operator<<(std::ostream& os, const TokVal& t) {
  return std::visit([&](auto&& val){return os << *val;}, t.value);
}
)--";

  return EXIT_SUCCESS;
}
