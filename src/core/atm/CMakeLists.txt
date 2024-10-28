add_library(atm STATIC atm.cpp)

target_link_libraries(atm PUBLIC matpack operators quantum species physics)
target_include_directories(atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)