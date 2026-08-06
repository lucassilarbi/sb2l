/**
 * sb2l interactive editor: GLFW + OpenGL3 + Dear ImGui host.
 *
 * Run with no argument for the interactive editor. `--shot <file.ppm>` renders
 * the default scene and writes one PPM screenshot, then exits; `--shot3` does
 * the same on the 3D canvas (used to produce the documentation images).
 */
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "app.hpp"
#include "editor.hpp"

static void glfw_error(int e, const char* d) { std::fprintf(stderr, "GLFW %d: %s\n", e, d); }

// Dump the current framebuffer as a binary PPM (flipped: GL rows are bottom-up).
static void write_ppm(GLFWwindow* window, const std::string& path)
{
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    std::vector<unsigned char> px((size_t)w * h * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, px.data());
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", path.c_str()); return; }
    std::fprintf(f, "P6\n%d %d\n255\n", w, h);
    for (int y = h - 1; y >= 0; --y)
        std::fwrite(&px[(size_t)y * w * 3], 1, (size_t)w * 3, f);
    std::fclose(f);
    std::printf("wrote %s\n", path.c_str());
}

int main(int argc, char** argv)
{
    const char* shot_path = nullptr;
    bool shot_3d = false;
    if (argc >= 3 && std::strcmp(argv[1], "--shot") == 0) shot_path = argv[2];
    if (argc >= 3 && std::strcmp(argv[1], "--shot3") == 0) { shot_path = argv[2]; shot_3d = true; }

    glfwSetErrorCallback(glfw_error);
    if (!glfwInit()) return 1;

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(1280, 800, "sb2l B-spline editor", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    sb2gui::App app;
    if (shot_3d) app.set_dim(3);
    sb2gui::Editor editor(app);

    int frame = 0;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        editor.draw_controls_window();
        editor.draw_canvas_window();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.10f, 0.10f, 0.11f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);

        // A couple of frames let the layout settle: the canvas rect, hence the
        // view fit, is only known once the windows have been laid out.
        if (shot_path && ++frame >= 3) {
            write_ppm(window, shot_path);
            break;
        }
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
