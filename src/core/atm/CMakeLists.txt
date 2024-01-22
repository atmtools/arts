add_library(atm STATIC atm.cpp)

target_link_libraries(atm PUBLIC matpack quantum species)
target_include_directories(atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)