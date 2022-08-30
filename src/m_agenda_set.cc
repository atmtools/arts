#include <default_gins.h>

#include "agenda_class.h"
#include "arts_options.h"
#include "global_data.h"

//! A Complete Description of a Method Variable
struct AgendaMethodVariable {
  String name{};
  Index group{-1};
  Index ws_pos{-1};
  Index g_pos{-1};
  bool input{false};
  bool any{false};

  //! DEBUG
  friend std::ostream& operator<<(std::ostream& os,
                                  const AgendaMethodVariable& x) {
    return os << x.name << "@" << x.ws_pos;
  }
};

//! Return a full list of all the Method Variables with defaults on the workspace
Array<AgendaMethodVariable> sorted_mdrecord(Workspace& ws,
                                            const String& method_name) {
  static const Index any_pos = global_data::WsvGroupMap.at("Any");

  auto& method =
      global_data::md_data_raw.at(global_data::MdRawMap.at(method_name));

  Array<AgendaMethodVariable> var_order;

  for (auto& output : method.Out()) {
    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data_ptr->operator[](output).Name();
    var.group = ws.wsv_data_ptr->operator[](output).Group();
    var.ws_pos = output;
    var.input = false;
  }

  for (Index i = 0; i < method.GOut().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = method.GOut()[i];
    var.group = method.GOutType()[i];
    var.g_pos = i;
    var.input = false;
    var.any = var.group == any_pos;
  }

  for (auto& input : method.In()) {
    // Skip inout
    if (std::any_of(method.Out().begin(),
                    method.Out().end(),
                    [input](auto& output) { return input == output; }))
      continue;

    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data_ptr->operator[](input).Name();
    var.group = ws.wsv_data_ptr->operator[](input).Group();
    var.ws_pos = input;
    var.input = true;
  }

  for (Index i = 0; i < method.GIn().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = method.GIn()[i];
    var.group = method.GInType()[i];
    var.g_pos = i;
    var.input = true;
    var.any = var.group == any_pos;

    // Create defaults (on the workspace)
    if (method.GInDefault()[i] not_eq NODEF) {
      var.ws_pos =
          create_workspace_gin_default_internal(ws, method_name, var.name);
    }
  }

  return var_order;
}

//! Split input and output of method variables
std::pair<ArrayOfIndex, ArrayOfIndex> split_io(
    Array<AgendaMethodVariable>& var_order) {
  ArrayOfIndex in, out;
  for (auto& var : var_order) {
    if (var.input)
      in.push_back(var.ws_pos);
    else
      out.push_back(var.ws_pos);
  }
  return {in, out};
}

template <typename T>
concept ArtsType =
    not std::is_array_v<T> and
    not std::is_same_v<std::remove_cvref_t<std::remove_pointer_t<T>>, char> and
    requires(T a) {
  TokVal{a};
};

struct SetWsv {
  enum class opt : char { NameOnly, ValueOnly, NameAndValue };

  opt test;
  std::string str{""};
  TokVal val{};

  SetWsv(const std::string& x, const std::string& y)
      : test(opt::NameOnly), str(x + "=" + y) {}
  SetWsv(std::string_view x) : test(opt::NameOnly), str(x) {}
  SetWsv(const ArtsType auto& t) : test(opt::ValueOnly), val(t) {}
  SetWsv(std::string_view x, const ArtsType auto& t)
      : test(opt::NameAndValue), str(x), val(t) {}
};

struct MethodVariable {
  bool positional{false};

  Index method_pos{-1};
  Index ws_pos{-1};

  bool new_value{false};

  constexpr MethodVariable() {}

  //! Fully named parameters; fun(X) or fun(X = Y)
  MethodVariable(Workspace& ws,
                 const Array<AgendaMethodVariable>& list,
                 const SetWsv& wsv) {
    using enum SetWsv::opt;
    switch (wsv.test) {
      case NameOnly: {
        auto& expr = wsv.str;
        auto equal_sign = std::find(expr.begin(), expr.end(), '=');
        if (equal_sign == expr.end()) {
          positional = true;
          ws_pos = ws.WsvMap_ptr->at(expr);

        } else {
          positional = false;
          auto lhs = std::string_view(expr.begin(), equal_sign);
          auto rhs = std::string_view(equal_sign, expr.end());

          while (std::isspace(lhs.front()) or lhs.front() == '=')
            lhs.remove_prefix(1);
          while (std::isspace(lhs.back()) or lhs.back() == '=')
            lhs.remove_suffix(1);
          while (std::isspace(rhs.front()) or rhs.front() == '=')
            lhs.remove_prefix(1);
          while (std::isspace(rhs.back()) or rhs.back() == '=')
            lhs.remove_suffix(1);
          method_position(list, lhs);
          ws_pos = ws.WsvMap_ptr->at(rhs);
        }
      } break;
      case ValueOnly:
        positional = true;
        wsv_position(ws, wsv.val);
        break;
      case NameAndValue:
        positional = false;
        method_position(list, wsv.str);
        wsv_position(ws, wsv.val);
        break;
    }
  }

  void add_del(Workspace& ws, Agenda& a) const {
    if (new_value)
      a.push_back(MRecord(
          global_data::MdMap.at(var_string(
              "Delete_sg_",
              global_data::wsv_groups.at(ws.wsv_data_ptr->at(ws_pos).Group()))),
          {},
          {ws_pos},
          {},
          a));
  }

  void add_set(Workspace& ws, Agenda& a) const {
    if (new_value)
      a.push_back(MRecord(
          global_data::MdMap.at(var_string(
              global_data::wsv_groups.at(ws.wsv_data_ptr->at(ws_pos).Group())
                  .name,
              "Set")),
          {ws_pos},
          {},
          ws.wsv_data_ptr->at(ws_pos).default_value(),
          a));
  }

 private:
  void method_position(const Array<AgendaMethodVariable>& list,
                       const std::string_view rhs) {
    auto ptr = std::find_if(
        list.begin(), list.end(), [rhs](auto& x) { return x.name == rhs; });
    ARTS_ASSERT(
        ptr not_eq list.end(), "Wrongly named parameter selection: ", rhs);
    method_pos = std::distance(list.begin(), ptr);
  }

  void wsv_position(Workspace& ws, const TokVal& value) {
    new_value = true;
    static std::size_t n = 0;
    ws_pos = ws.add_wsv(WsvRecord(
        var_string("::wsv", n++).c_str(), "method value", value.type(), value));
  }
};

template <std::size_t N>
std::array<MethodVariable, N> input_data_array(
    Workspace& ws,
    const Array<AgendaMethodVariable>& list,
    const std::array<SetWsv, N>& vals) {
  std::array<MethodVariable, N> out{};
  for (std::size_t i = 0; i < N; i++)
    out[i] = MethodVariable(ws, list, vals[i]);
  return out;
}

template <std::size_t N>
const TokVal& to_tokval(Workspace& ws,
                        const std::array<MethodVariable, N>& input_data) {
  static TokVal any{};
  if (auto ptr = std::find_if(input_data.crbegin(),
                              input_data.crend(),
                              [](auto& mv) { return mv.new_value; });
      ptr not_eq input_data.crend())
    return ws.wsv_data_ptr->at(ptr->ws_pos).default_value();
  return any;
}

template <std::size_t N>
void add_method_and_setters(Workspace& ws,
                            Agenda& a,
                            const std::string_view method,
                            const std::array<SetWsv, N>& input) {
  auto list = sorted_mdrecord(ws, method);

  // List of all input data
  auto input_data = input_data_array(ws, list, input);

  // Adapt positional arguments
  for (Index i = 0; i < static_cast<Index>(input_data.size()); i++) {
    if (input_data[i].positional) input_data[i].method_pos = i;
  }

  for (auto& x : input_data) list[x.method_pos].ws_pos = x.ws_pos;

  for (auto& io : list)
    ARTS_USER_ERROR_IF(
        io.ws_pos < 0,
        "Not setting all input, this is a developer error.  Here is the data:\n",
        list)

  auto [in, out] = split_io(list);

  auto ptr = std::find_if(global_data::md_data.begin(),
                          global_data::md_data.end(),
                          [method](auto& x) { return method == x.Name(); });
  const auto m_id = std::distance(global_data::md_data.begin(), ptr);

  if (not ptr->SetMethod())
    for (auto& x : input_data) x.add_set(ws, a);
  a.push_back(MRecord(m_id, out, in, to_tokval(ws, input_data), a));
  if (not ptr->SetMethod())
    for (auto& x : input_data) x.add_del(ws, a);
}

//! Helper class to create an agenda
struct AgendaCreator {
  Workspace& ws;
  const Verbosity& v;
  Agenda agenda;

  AgendaCreator(Workspace& workspace,
                const char* name,
                const Verbosity& verbosity)
      : ws(workspace), v(verbosity), agenda(ws) {
    agenda.set_name(name);
  }

  //! Check the agenda and move it out of here
  Agenda finalize() {
    agenda.check(ws, v);
    return std::move(agenda);
  }

  //! Add a method with as many inputs as you want.  These inputs must be of type Wsv
  template <typename... Input>
  void add(const std::string_view method, Input... input) {
    if constexpr (sizeof...(Input) > 0)
      add_method_and_setters(
          ws, agenda, method, stdarrayify<SetWsv>(SetWsv(input)...));
    else
      add_method_and_setters(ws, agenda, method, std::array<SetWsv, 0>{});
  }
};

void iy_main_agendaSet(Workspace& ws,
                       Agenda& iy_main_agenda,
                       const String& option,
                       const Verbosity& verbosity) {
  AgendaCreator agenda(ws, "iy_main_agenda", verbosity);

  using enum Options::iy_main_agendaDefaultOptions;
  switch (Options::toiy_main_agendaDefaultOptionsOrThrow(option)) {
    case Emission:
      agenda.add("ppathCalc");
      agenda.add("iyEmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Clearsky:
      agenda.add("ppathCalc");
      agenda.add("iyClearsky");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Transmission:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case TransmissionUnitUnpolIntensity:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add(
          "MatrixUnitIntensity", "iy_transmitter", "stokes_dim", "f_grid");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case TransmissionUnitPolIntensity:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add("iy_transmitterSinglePol");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Freqloop:
      agenda.add("Ignore", "diy_dx");
      agenda.add("Ignore", "iy_id");
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "nlte_field");
      agenda.add("Ignore", "cloudbox_on");
      agenda.add("Ignore", "jacobian_do");
      agenda.add("iyLoopFrequencies");
      agenda.add("Touch", "ppath");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case ScattMC:
      agenda.add("Ignore", "rte_pos2");
      agenda.add("Ignore", "diy_dx");
      agenda.add("Ignore", "iy_id");
      agenda.add("Ignore", "nlte_field");
      agenda.add("iyMC");
      agenda.add("Touch", "ppath");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case FINAL:
      break;
  }

  iy_main_agenda = agenda.finalize();
}

void iy_loop_freqs_agendaSet(Workspace& ws,
                             Agenda& iy_loop_freqs_agenda,
                             const String& option,
                             const Verbosity& verbosity) {
  AgendaCreator agenda(ws, "iy_loop_freqs_agenda", verbosity);

  using enum Options::iy_loop_freqs_agendaDefaultOptions;
  switch (Options::toiy_loop_freqs_agendaDefaultOptionsOrThrow(option)) {
    case Emission:
      agenda.add("ppathCalc");
      agenda.add("iyEmissionStandard");
      break;
    case Transmission:
      agenda.add("Ignore", "iy_id");
      agenda.add("ppathCalc");
      agenda.add("iyTransmissionStandard");
      break;
    case FINAL:
      break;
  }

  iy_loop_freqs_agenda = agenda.finalize();
}

void iy_space_agendaSet(Workspace& ws,
                        Agenda& iy_space_agenda,
                        const String& option,
                        const Verbosity& verbosity) {
  AgendaCreator agenda(ws, "iy_space_agenda", verbosity);

  using enum Options::iy_space_agendaDefaultOptions;
  switch (Options::toiy_space_agendaDefaultOptionsOrThrow(option)) {
    case CosmicBackground:
      agenda.add("Ignore", "rtp_pos");
      agenda.add("Ignore", "rtp_los");
      agenda.add("MatrixCBR", "iy", "stokes_dim", "f_grid");
      break;
    case FINAL:
      break;
  }

  iy_space_agenda = agenda.finalize();
}

void iy_surface_agendaSet(Workspace& ws,
                          Agenda& iy_surface_agenda,
                          const String& option,
                          const Verbosity& verbosity) {
  AgendaCreator agenda(ws, "iy_surface_agenda", verbosity);

  using enum Options::iy_surface_agendaDefaultOptions;
  switch (Options::toiy_surface_agendaDefaultOptionsOrThrow(option)) {
    case UseSurfaceRtprop:
      agenda.add("SurfaceDummy");
      agenda.add("iySurfaceRtpropAgenda");
      break;
    case FINAL:
      break;
  }

  iy_surface_agenda = agenda.finalize();
}

void iy_cloudbox_agendaSet(Workspace& ws,
                           Agenda& iy_cloudbox_agenda,
                           const String& option,
                           const Verbosity& verbosity) {
  AgendaCreator agenda(ws, "iy_cloudbox_agenda", verbosity);

  using enum Options::iy_cloudbox_agendaDefaultOptions;
  switch (Options::toiy_cloudbox_agendaDefaultOptionsOrThrow(option)) {
    case LinInterpField:
      agenda.add("iyInterpCloudboxField");
      break;
    case QuarticInterpField:
      agenda.add("iyInterpCloudboxField", SetWsv{"za_interp_order", Index{4}});
      break;
    case FINAL:
      break;
  }

  iy_cloudbox_agenda = agenda.finalize();
}
