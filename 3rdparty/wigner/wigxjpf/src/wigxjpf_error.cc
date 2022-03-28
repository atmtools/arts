#include <stdexcept>

extern "C" {
void cpp_wigxjpf_throw(void) {
    throw std::runtime_error("\n"
	  "wigxjpf: Error detected! "
	  "** Library misuse?  See documentation. **\n"
	  "\n");
}
}