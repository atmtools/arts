add_library(workspace_group STATIC
  workspace_groups.cpp
)

add_executable(make_auto_wsg make_auto_wsg.cpp)
add_custom_command(OUTPUT auto_wsg.h auto_wsg.cpp
                   COMMAND make_auto_wsg > auto_wsg.h 2> auto_wsg.cpp
                   DEPENDS make_auto_wsg)
target_link_libraries(make_auto_wsg PRIVATE workspace_group)

add_library(workspace_internals STATIC
  workspace_agendas.cpp
  workspace_methods.cpp
  workspace_variables.cpp
  auto_wsg.cpp
)
target_link_libraries(workspace_internals PUBLIC workspace_group matpack predef species fwd)
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(workspace_internals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(make_auto_wsv make_auto_wsv.cpp)
add_custom_command(OUTPUT auto_wsv.h auto_wsv.cpp
                   COMMAND make_auto_wsv > auto_wsv.h 2> auto_wsv.cpp
                   DEPENDS make_auto_wsv)
target_link_libraries(make_auto_wsv PRIVATE workspace_internals)

add_executable(make_auto_wsm make_auto_wsm.cpp auto_wsv.cpp)
add_custom_command(OUTPUT auto_wsm.h auto_wsm.cpp
                   COMMAND make_auto_wsm > auto_wsm.h 2> auto_wsm.cpp
                   DEPENDS make_auto_wsm)
target_link_libraries(make_auto_wsm PRIVATE workspace_internals)

add_executable(make_auto_wsa make_auto_wsa.cpp)
add_custom_command(OUTPUT auto_wsa.h auto_wsa.cpp
                   COMMAND make_auto_wsa > auto_wsa.h 2> auto_wsa.cpp
                   DEPENDS make_auto_wsa)
target_link_libraries(make_auto_wsa PRIVATE workspace_internals)

add_library(workspace STATIC
  workspace_agenda_class.cpp
  workspace_agenda_creator.cpp
  workspace_method_class.cpp
  workspace_class.cpp
  auto_wsv.cpp
  auto_wsm.cpp
  auto_wsa.cpp
)
target_link_libraries(workspace PUBLIC workspace_internals)
target_include_directories(workspace PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(workspace PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
