#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

Agenda::Agenda(std::string n) : name(std::move(n)), methods{} {}
