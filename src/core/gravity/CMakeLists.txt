add_library(gravity STATIC
  gravity.cpp
)

target_link_libraries(gravity PUBLIC matpack)
target_include_directories(atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})