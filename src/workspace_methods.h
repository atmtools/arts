#pragma once

#include <wsv_value_wrapper.h>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceMethodInternalRecord {
  std::string                     desc;
  std::vector<std::string>        author;
  std::string                     return_type{"void"};
  std::string                     return_desc{};
  std::vector<std::string>        out{};
  std::vector<std::string>        gout{};
  std::vector<std::string>        gout_type{};
  std::vector<std::string>        gout_desc{};
  std::vector<std::string>        in{};
  std::vector<std::string>        gin{};
  std::vector<std::string>        gin_type{};
  ArrayOfArrayOfIndex             python_generic_sorting{};
  std::vector<std::optional<Wsv>> gin_value{};
  std::vector<std::string>        gin_desc{};
  bool                            pass_workspace{false};

  /*! Set to false to generate no size verification for this method.
   *
   * The generated checks read a variable's size as an invariant: what it is on
   * the way in is what its dimensions say.  A method that instead reads the
   * size as a message from its caller, e.g. *OEM*, where an empty
   * *model_state_vec* asks it to start from *model_state_vec_apriori*, has no
   * such invariant and cannot be verified this way.
   *
   * Turn it off for the whole method rather than for one variable, since a
   * method that negotiates sizes with its caller tends to do so for the group
   * of variables that travel together.
   */
  bool size_constraints{true};

  [[nodiscard]] static std::string generic_type(const std::string&, bool output = false);
  [[nodiscard]] std::string        header(const std::string& name) const;
  [[nodiscard]] std::string        call(const std::string& name) const;
};

const std::unordered_map<std::string, WorkspaceMethodInternalRecord>& internal_workspace_methods();

template <> struct std::formatter<WorkspaceMethodInternalRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const WorkspaceMethodInternalRecord& wsm, FmtContext& ctx) const {
    return tags.format(ctx,
                       "WorkspaceMethodInternalRecord{"sv,
                       "\n  .desc="sv,
                       wsm.desc,
                       "\n  .author="sv,
                       wsm.author,
                       "\n  .return_type="sv,
                       wsm.return_type,
                       "\n  .return_desc="sv,
                       wsm.return_desc,
                       "\n  .out="sv,
                       wsm.out,
                       "\n  .gout="sv,
                       wsm.gout,
                       "\n  .gout_type="sv,
                       wsm.gout_type,
                       "\n  .gout_desc="sv,
                       wsm.gout_desc,
                       "\n  .in="sv,
                       wsm.in,
                       "\n  .gin="sv,
                       wsm.gin,
                       "\n  .gin_type="sv,
                       wsm.gin_type,
                       "\n  .python_generic_sorting="sv,
                       wsm.python_generic_sorting,
                       "\n  .gin_value="sv,
                       wsm.gin_value,
                       "\n  .gin_desc="sv,
                       wsm.gin_desc,
                       "\n  .pass_workspace="sv,
                       wsm.pass_workspace ? "true"sv : "false"sv,
                       "\n  .size_constraints="sv,
                       wsm.size_constraints ? "true"sv : "false"sv,
                       "\n}"sv);
  }
};
