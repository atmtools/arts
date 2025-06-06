#pragma once

#include <string>
#include <vector>

struct EnumeratedOption {
  std::string name;
  std::string desc;
  std::vector<std::vector<std::string>> values_and_desc;
  std::size_t preferred_print{0};

  [[nodiscard]] std::string_view sz() const;
  [[nodiscard]] std::string head() const;
  [[nodiscard]] std::string tail() const;
  [[nodiscard]] std::string impl() const;
  [[nodiscard]] std::string docs() const;
};

const std::vector<EnumeratedOption>& internal_options();
