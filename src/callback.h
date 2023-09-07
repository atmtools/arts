#ifndef callback_h
#define callback_h

class Workspace;

#include <functional>
#include <ostream>
#include <string>
#include <vector>

struct CallbackOperator {
  std::function<void(const std::shared_ptr<Workspace>&)> callback{};
  std::vector<std::string> inputs{};
  std::vector<std::string> outputs{};

  void operator()(Workspace& ws) const;

  friend std::ostream& operator<<(std::ostream& os, const CallbackOperator&);
};

#endif
