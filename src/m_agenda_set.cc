#include "agenda_class.h"
#include "global_data.h"

struct MethodVariable {
  Workspace& ws;
  Agenda& agenda;
  Index ws_pos{-1};
  Index method_pos{-1};
  bool output{false};
  bool set{false};
  bool g{false};

  const auto& type() const {
    return global_data::wsv_groups.at(ws.wsv_data_ptr->at(ws_pos).Group()).name;
  }

  const auto& name() const { return ws.wsv_data_ptr->at(ws_pos).Name(); }

  void set_value(const TokVal& val) {
    agenda.push_back(MRecord(
        global_data::MdMap.at(type() + "Set"), {ws_pos}, {}, val, Agenda{ws}));
  }

  ~MethodVariable() {
    if (set) {
      agenda.push_back(MRecord(global_data::MdMap.at("Delete_sg_" + type()),
                               {},
                               {ws_pos},
                               {},
                               Agenda{ws}));
    }
  }

  Index new_wsv(const std::string_view t) {
    static std::size_t autoname{0};

    set = true;
    return ws.add_wsv(
        WsvRecord(("::", autoname++).c_str(), "Added automatically", t))
  }

  void from_method(const MdRecord& method, const std::string_view var_name) {
    if (auto ptr gin = method.GIn().find(var_name);
        gin not_eq method.GIn().end()) {
      method_pos = std::distance(method.GIn().begin(), gin);
      input = true;
      g = true;
    } else if (auto ptr gout = method.GOut().find(var_name);
               gout not_eq method.GOut().end()) {
      method_pos = std::distance(method.GOut().begin(), gout);
      g = true;
    } else if (auto out = method.Out().find(ws_pos);
               out not_eq method.Out().end()) {
      method_pos = std::distance(method.Out().begin(), out);
    } else {
      auto inonly = method.InOnly().find(ws_pos);
      method_pos = std::distance(method.InOnly().begin(), inonly);
      input = true;
    }
  }

  //! Set a method variable to a workspace variable
  MethodVariable(Workspace& w,
                 Agenda& a,
                 const MdRecord& method,
                 const std::string_view method_variable_name,
                 const std::string_view ws_variable_name)
      : ws(w),
        agenda(a),
        ws_pos(ws.WsvMap_ptr->at(ws_variable_name)),
        method_pos(from_method(method, method_variable_name)) {}

  //! Set a positional workspace variable
  MethodVariable(Workspace& w,
                 Agenda& a,
                 const MdRecord& method,
                 const std::string_view ws_variable_name const Index pos)
      : ws(w), agenda(a), ws_pos(ws.WsvMap_ptr->at(ws_variable_name)) {
    std::array<Index, 4> lens{method.Out().nelem(),
                              method.GOut().nelem(),
                              method.InOnly().nelem(),
                              method.Gin().nelem()};
    for (auto& n : lens) {
      if (pos - n >= 0)
        pos -= n;
      else
        break;
    }

    method_pos = pos;
  }

  //! Set from a TokVal
  MethodVariable(Workspace& w,
                 Agenda& a,
                 const MdRecord& method,
                 const std::string_view method_variable_name,
                 const TokVal& val) ws(w),
      agenda(a), ws_pos(new_wsv(val.type())),
      method_pos(from_method(method, method_variable_name)) set(true) {}
};
