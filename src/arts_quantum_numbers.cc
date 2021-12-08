#include "quantum_numbers.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char**argv) try {
    //! Tests to be removed
    std::cout << Quantum::Number::Value(" J  1 3") << '\n';
    std::cout << "X    1   3123 1/3 1,,,,,2  " << '\n';
    std::cout << Quantum::Number::items<3>("X    1   3123 1/3 1,,,,,2  ", 0) << '\n';
    std::cout << Quantum::Number::items<3>("X    1   3123 1/3 1,,,,,2  ", 2) << '\n';
    std::cout << Quantum::Number::items<3>("X    1   3123 1/3 1,,,,,2  ", 4) << '\n';
    std::cout << Quantum::Number::items<3>("X    1   3123 1/3 1,,,,,2  ", 5) << '\n';

    String s;
    for (Index i=1; i<argc; i++) {s += argv[i]; s += " ";}
    std::cout << "Read: " << Quantum::Number::ValueList(s) << '\n';

auto g = QuantumIdentifier{1, Quantum::Number::ValueList(s)};
    std::cout << Quantum::Number::checkLocalGlobal(g, Quantum::Number::LocalState{Quantum::Number::ValueList{"Ka 1 3"}}, g, true) << '\n';
    return EXIT_SUCCESS;
} catch(std::runtime_error&e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
}
