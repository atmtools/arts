if (ENABLE_GUI)
  cmake_minimum_required(VERSION 3.10)
  project(imgui_interface)

  add_subdirectory(imgui)
  add_subdirectory(implot)

  add_executable(imgui_example imgui_example.cpp)
  target_include_directories(imgui_example PUBLIC ${PROJECT_SOURCE_DIR})
  target_link_libraries(imgui_example imgui)
endif()


if (ENABLE_GUI)
  message (STATUS "GUI enabled  (use -DENABLE_GUI=0 to disable)")
else()
  message (STATUS "GUI disabled (use -DENABLE_GUI=1 to enable)")
endif()
