#ifndef gui_macros_h
#define gui_macros_h

#include <GL/glew.h>     // Initialize with glewInit()
#include <GLFW/glfw3.h>  // Include glfw3.h after our OpenGL definitions
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_stdlib.h>
#include <implot.h>
#include <implot_internal.h>
#include <stdio.h>

#include <iostream>

extern "C" {
/** A C error handling for GLFW */
inline static void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}
}

/** Begins an endless GUI loop (Opens a bracked)
 *
 * Expects:
 *  window: A working glfw window
 */
#define BeginWhileLoopGUI                  \
  while (!glfwWindowShouldClose(window)) { \
    glfwPollEvents();                      \
                                           \
    ImGui_ImplOpenGL3_NewFrame();          \
    ImGui_ImplGlfw_NewFrame();             \
    ImGui::NewFrame()

/** Ends an endless GUI loop (Closes a bracked)
 *
 * Expects:
 *  window: A working glfw window
 */
#define EndWhileLoopGUI                                   \
  ImGui::Render();                                        \
  int display_w, display_h;                               \
  glfwGetFramebufferSize(window, &display_w, &display_h); \
  glViewport(0, 0, display_w, display_h);                 \
  glClear(GL_COLOR_BUFFER_BIT);                           \
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); \
  glfwSwapBuffers(window);                                \
  }

/** Sets the glsl_version variable and hints for the windows
 *
 * System dependent...
 *
 * Sets:
 *  glsl_version: Version of GLSL
 */
#if __APPLE__
#define Selectglsl_versionGUI                                    \
  const char *glsl_version = "#version 150";                     \
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);                 \
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);                 \
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); \
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
#define Selectglsl_versionGUI                    \
  const char *glsl_version = "#version 130";     \
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); \
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif

/** Initialize the GUI
 *
 * Should set up all things needed to begin the main loop of drawing things
 *
 * Sets:
 *  window: A glfw window
 *  glsl_version: Version of GLSL
 */
#define InitializeGUI(NAME, WIDTH, HEIGHT)                                    \
  glfwSetErrorCallback(glfw_error_callback);                                  \
  if (!glfwInit()) throw std::runtime_error("Cannot initialize a window");    \
                                                                              \
  Selectglsl_versionGUI                                                       \
                                                                              \
      GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, NAME, NULL, NULL); \
  if (window == NULL) throw std::runtime_error("Cannot create a window");     \
  glfwMakeContextCurrent(window);                                             \
  glfwSwapInterval(1);                                                        \
                                                                              \
  if (glewInit() != GLEW_OK)                                                  \
    throw std::runtime_error("Cannot initialize OpenGL loader");              \
                                                                              \
  IMGUI_CHECKVERSION();                                                       \
  ImGui::CreateContext();                                                     \
  ImPlot::CreateContext();                                                    \
  ImGui::StyleColorsDark();                                                   \
                                                                              \
  ImGui_ImplGlfw_InitForOpenGL(window, true);                                 \
  ImGui_ImplOpenGL3_Init(glsl_version)

/** Cleanup for the GUI
 *
 * Should clean up all things InitializeGUI
 *
 * Destroys:
 *  window: The glfw window
 */
#define CleanupGUI              \
  ImGui_ImplOpenGL3_Shutdown(); \
  ImGui_ImplGlfw_Shutdown();    \
  ImPlot::DestroyContext();     \
  ImGui::DestroyContext();      \
                                \
  glfwDestroyWindow(window);    \
  glfwTerminate()

#endif  // gui_macros_h
