#pragma once

#include "agenda_class.h"
#include "array.h"
#include "global_data.h"
#include "matpack.h"
#include "mystring.h"
#include "workspace_ng.h"

namespace AgendaManip {
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
                                  const AgendaMethodVariable& x);
};

//! Return a full list of all the Method Variables with defaults on the workspace
Array<AgendaMethodVariable> sorted_mdrecord(Workspace& ws,
                                            const String& method_name);

//! Split input and output of method variables
std::pair<ArrayOfIndex, ArrayOfIndex> split_io(
    Array<AgendaMethodVariable>& var_order);

struct SetWsv {
  //! Checks for which type of information should be used later on
  enum class opt : char { NameOnly, ValueOnly, NameAndValue };

  opt test;
  std::string str{""};
  TokVal val{};

  //! For named wsv1=wsv2
  SetWsv(const std::string& x, const std::string& y);

  //! For either named "wsv1=wsv2" or positional wsv1
  SetWsv(std::string_view x);

  //! For position value
  SetWsv(ArtsType auto&& t) : test(opt::ValueOnly), val(std::forward<decltype(t)>(t)) {}

  //! For named value
  SetWsv(std::string_view x, ArtsType auto&& t) : test(opt::NameAndValue), str(x), val(std::forward<decltype(t)>(t)) {}
};

struct MethodVariable {
  bool positional{false};

  Index method_pos{-1};
  Index ws_pos{-1};

  bool new_value{false};

  constexpr MethodVariable() {}

  MethodVariable(Workspace& ws,
                 const Array<AgendaMethodVariable>& list,
                 const SetWsv& wsv);

  void add_del(Workspace& ws, Agenda& a) const;

  void add_set(Workspace& ws, Agenda& a) const;

 private:
  void method_position(const Array<AgendaMethodVariable>& list,
                       const std::string_view rhs);

  void wsv_position(Workspace& ws, const TokVal& value);
};

template <std::size_t N>
std::array<MethodVariable, N> input_data_array(
    Workspace& ws,
    const Array<AgendaMethodVariable>& list,
    const std::array<SetWsv, N>& vals) {
  std::array<MethodVariable, N> out{};
  if constexpr (N > 0)
    for (std::size_t i = 0; i < N; i++)
      out[i] = MethodVariable(ws, list, vals[i]);
  return out;
}

template <std::size_t N>
const TokVal& to_tokval(Workspace& ws,
                        const std::array<MethodVariable, N>& input_data) {
  const static TokVal any{};
  if constexpr (N > 0)
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
  if constexpr (N > 0) {
    for (Index i = 0; i < static_cast<Index>(N); i++) {
      if (input_data[i].positional) input_data[i].method_pos = i;
    }
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
  a.push_back(MRecord(m_id, out, in, to_tokval(ws, input_data), Agenda{ws}));
  if (not ptr->SetMethod())
    for (auto& x : input_data) x.add_del(ws, a);
}

//! Helper class to create an agenda
struct AgendaCreator {
  Workspace& ws;
  Agenda agenda;

  AgendaCreator(Workspace& workspace, const char* name);

  //! Check the agenda and return a copy of it (ignoring/touching all agenda input/output not dealt with)
  Agenda finalize();

  //! Add a method with as many inputs as you want.  These inputs must be of type Wsv
  template <typename... Input>
  void add(const std::string_view method, Input&&... input) {
    if constexpr (sizeof...(Input) > 0)
      add_method_and_setters(
          ws,
          agenda,
          method,
          stdarrayify<SetWsv>(SetWsv(std::forward<Input>(input))...));
    else
      add_method_and_setters(ws, agenda, method, std::array<SetWsv, 0>{});
  }

  //! Set a variable to a value
  void set(const std::string_view var, const TokVal& value);

  //! Ignores a variable (only call if the variable input is ignored despite being BOTH in- and output)
  void ignore(const std::string_view var);
};

Agenda get_iy_main_agenda(Workspace& ws, const String& option);
Agenda get_iy_loop_freqs_agenda(Workspace& ws, const String& option);
Agenda get_iy_space_agenda(Workspace& ws, const String& option);
Agenda get_iy_surface_agenda(Workspace& ws, const String& option);
Agenda get_iy_cloudbox_agenda(Workspace& ws, const String& option);
Agenda get_ppath_agenda(Workspace& ws, const String& option);
Agenda get_ppath_step_agenda(Workspace& ws, const String& option);
Agenda get_refr_index_air_agenda(Workspace& ws, const String& option);
Agenda get_water_p_eq_agenda(Workspace& ws, const String& option);
Agenda get_gas_scattering_agenda(Workspace& ws, const String& option);
Agenda get_surface_rtprop_agenda(Workspace& ws, const String& option);
Agenda get_g0_agenda(Workspace& ws, const String& option);
Agenda get_dobatch_calc_agenda(Workspace& ws, const String& option);
Agenda get_ybatch_calc_agenda(Workspace& ws, const String& option);
Agenda get_test_agenda(Workspace& ws, const String& option);
Agenda get_surface_rtprop_sub_agenda(Workspace& ws, const String& option);
Agenda get_spt_calc_agenda(Workspace& ws, const String& option);
Agenda get_sensor_response_agenda(Workspace& ws, const String& option);
Agenda get_propmat_clearsky_agenda(Workspace& ws, const String& option);
Agenda get_pha_mat_spt_agenda(Workspace& ws, const String& option);
Agenda get_met_profile_calc_agenda(Workspace& ws, const String& option);
Agenda get_main_agenda(Workspace& ws, const String& option);
Agenda get_jacobian_agenda(Workspace& ws, const String& option);
Agenda get_iy_radar_agenda(Workspace& ws, const String& option);
Agenda get_iy_independent_beam_approx_agenda(Workspace& ws, const String& option);
Agenda get_inversion_iterate_agenda(Workspace& ws, const String& option);
Agenda get_forloop_agenda(Workspace& ws, const String& option);
Agenda get_doit_scat_field_agenda(Workspace& ws, const String& option);
Agenda get_doit_rte_agenda(Workspace& ws, const String& option);
Agenda get_doit_mono_agenda(Workspace& ws, const String& option);
Agenda get_doit_conv_test_agenda(Workspace& ws, const String& option);
}  // namespace AgendaManip
