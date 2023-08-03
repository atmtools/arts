#pragma once

#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "auto_wsg.h"

class Agenda;
class Workspace;

class Method {
std::string name;
std::vector<std::string> outargs{};
std::vector<std::string> inargs{};
std::function<void(Workspace&, const std::vector<std::string>&, const std::vector<std::string>&)> func{};
std::optional<Wsv> setval{};

public:
Method() = default;
Method(const std::string& name, const std::vector<std::string>& args={}, const std::unordered_map<std::string, std::string>& kwargs={});
Method(std::string name, const Wsv& wsv);

[[nodiscard]] const std::string& get_name() const {return name;}
[[nodiscard]] const std::vector<std::string>& get_outs() const {return outargs;}
[[nodiscard]] const std::vector<std::string>& get_ins() const {return inargs;}
[[nodiscard]] const std::optional<Wsv>& get_setval() const {return setval;}

void operator()(Workspace& ws) const;
void agenda_setvals(Agenda& agenda, bool del) const;

friend std::ostream& operator<<(std::ostream& os, const Method& m);
};
