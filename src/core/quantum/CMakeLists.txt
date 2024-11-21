add_library(quantum STATIC
            hitran_species.cc
            quantum_numbers.cc
            quantum_term_symbol.cc
            )
target_link_libraries(quantum PUBLIC species matpack binio)
target_include_directories(quantum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
