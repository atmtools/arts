#include "quantum.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char**argv) try {
    String s;
    for (Index i=1; i<argc; i++) {s += argv[i]; s += " ";}
    std::println("Read: {}", QuantumIdentifier(s));
    return EXIT_SUCCESS;
} catch(std::runtime_error&e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
}
