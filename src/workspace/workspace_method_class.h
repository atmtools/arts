#pragma once

#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "auto_wsg.h"

class Agenda;

class Method {
std::string name{};
std::vector<std::string> outargs{};
std::vector<std::string> inargs{};
std::optional<Wsv> setval{std::nullopt};
bool overwrite_setval{false};

public:
Method();
Method(const std::string& name, const std::vector<std::string>& args, const std::unordered_map<std::string, std::string>& kwargs);
Method(std::string name, const Wsv& wsv, bool=false);
Method(const std::string& name, const std::vector<std::string>& ins, const std::vector<std::string>& outs, const std::optional<Wsv>& wsv, bool overwrite);

[[nodiscard]] const std::string& get_name() const {return name;}
[[nodiscard]] const std::vector<std::string>& get_outs() const {return outargs;}
[[nodiscard]] const std::vector<std::string>& get_ins() const {return inargs;}
[[nodiscard]] const std::optional<Wsv>& get_setval() const {return setval;}
[[nodiscard]] bool overwrite() const {return overwrite_setval;}

void operator()(Workspace& ws) const;
void add_defaults_to_agenda(Agenda& agenda) const;

friend std::ostream& operator<<(std::ostream& os, const Method& m);
};
