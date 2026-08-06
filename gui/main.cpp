/**
 * sb2l interactive editor: GLFW + OpenGL3 + Dear ImGui host.
 *
 * Usage:
 *   sb2l_gui                          interactive editor
 *   sb2l_gui --shot <file.ppm> [opt]  render one screenshot of that configuration, then exit
 *   sb2l_gui --tour <dir> [opt]       render the scripted feature tour as dir/f%04d.ppm
 *
 * Options (all optional, applied before the first frame):
 *   --dim 2|3      --ps R|IR|Z    --cs R|IR|Z    --ct UR|UN|CR|CN
 *   --p <n>        --n <n>        --d <n>        --gen <n>
 *   --size <f>     scale every control set (half-widths and generators)
 *   --yaw <rad>    --pitch <rad>  --banner "<text>"
 *
 * Both capture modes exist so every image of the documentation is
 * reproducible from a single command (see docs/make_assets.sh).
 */
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
}

// ------------------------------------------------------------------ options ---

static sb2l::ParameterSet parse_set(const char* s)
{
    if (std::strcmp(s, "R") == 0) return sb2l::ParameterSet::R;
    if (std::strcmp(s, "Z") == 0) return sb2l::ParameterSet::Z;
    return sb2l::ParameterSet::IR;
}

static sb2l::CurveType parse_ct(const char* s)
{
    if (std::strcmp(s, "UR") == 0) return sb2l::CurveType::UNIFORM_RATIONAL;
    if (std::strcmp(s, "UN") == 0) return sb2l::CurveType::UNIFORM_NONRATIONAL;
    if (std::strcmp(s, "CR") == 0) return sb2l::CurveType::CLAMPED_RATIONAL;
    return sb2l::CurveType::CLAMPED_NONRATIONAL;
}

// ---------------------------------------------------------------- the tour ---

// One scripted stage: a caption, a state applied once on entry, and an
// optional per-frame animation (the local frame index is passed in).
struct Stage {
    int frames;
    const char* caption;
    void (*enter)(sb2gui::App&);
    void (*step)(sb2gui::App&, sb2gui::Editor&, int);
};

// The whole tour runs on the tightest evaluation the library offers: an
// affine parameter, so the many occurrences of u inside a basis function stay
// correlated instead of being enclosed independently, and the Taylor form,
// which expands the basis about the midpoint of each subdivision with a
// rigorous remainder. Measured on the scene below, that pair gives a tube
// about a third narrower than the same curve over an interval parameter.
static void st_boxes(sb2gui::App& a)
{
    a.p = 3; a.nCP = 6; a.d = 20;
    a.ct = sb2l::CurveType::CLAMPED_NONRATIONAL;
    a.f = sb2l::Form::TAYLOR;
    a.ps = sb2l::ParameterSet::Z;
    a.cs = sb2l::ParameterSet::IR;
    a.rebuild();
    a.want_fit = true;
}

// Walk a control set around a small circle: the tube follows, and only the
// p+1 segments it belongs to are re-evaluated.
static void an_drag(sb2gui::App& a, sb2gui::Editor&, int k)
{
    const int i = 2;
    if (i >= (int)a.cps.size()) return;
    const double th = 2.0 * 3.14159265358979 * k / 30.0;
    a.cps[i].cx += 0.055 * std::cos(th);
    a.cps[i].cy += 0.055 * std::sin(th);
    a.update_control_point(i);
}

static void st_zono(sb2gui::App& a)
{
    a.cs = sb2l::ParameterSet::Z;
    a.nGen = 2;
    a.rebuild();
}

// Turn every control zonotope: a generator is a free vector, so this
// reorients the sets instead of merely rescaling them.
static void an_turn(sb2gui::App& a, sb2gui::Editor&, int)
{
    const double c = std::cos(0.11), s = std::sin(0.11);
    for (int i = 0; i < (int)a.cps.size(); ++i) {
        for (sb2gui::Vec3& g : a.cps[i].gens) {
            const double x = g.x * c - g.y * s, y = g.x * s + g.y * c;
            g.x = x; g.y = y;
        }
        a.update_control_point(i);
    }
}

static void st_curve(sb2gui::App& a)
{
    a.cs = sb2l::ParameterSet::R;
    a.rebuild();
}

// Rational weights are shown on real *control points*: the rational basis is
// a quotient, so it loses the partition-of-unity cancellation that keeps the
// non-rational tube tight, and combining it with control boxes widens the
// enclosure by an order of magnitude whatever the subdivision. Collapsing the
// basis to its midpoint sidesteps that, and the affine parameter of the tour
// is kept -- it shifts the samples by half a subdivision and nothing more.
static void st_rational(sb2gui::App& a)
{
    a.cs = sb2l::ParameterSet::R;
    a.ps = sb2l::ParameterSet::R; // the quotient basis over any set overestimates sharply
    a.ct = sb2l::CurveType::CLAMPED_RATIONAL;
    for (size_t i = 0; i < a.wnum.size(); ++i) { a.wnum[i] = 1; a.wden[i] = 1; }
    if (a.wnum.size() > 2) a.wnum[2] = 6;
    a.rebuild();
}

static void st_3d(sb2gui::App& a)
{
    a.ct = sb2l::CurveType::CLAMPED_NONRATIONAL;
    for (size_t i = 0; i < a.wnum.size(); ++i) a.wnum[i] = 1;
    a.cs = sb2l::ParameterSet::IR;
    a.ps = sb2l::ParameterSet::Z; // back to the affine parameter after the rational stage
    a.rebuild();
    a.set_dim(3);
}

static void st_3dz(sb2gui::App& a)
{
    a.cs = sb2l::ParameterSet::Z;
    a.nGen = 4;   // three generators span a volume, four make its facets read
    a.d = 10;     // fewer, fatter elements along the curve
    a.rebuild();
    // Enlarge the control sets: at their editing size the result elements are
    // legible polyhedra rather than thin lenses.
    for (sb2gui::ControlPoint& c : a.cps)
        for (sb2gui::Vec3& g : c.gens) g = 1.4 * g;
    a.reeval();
    a.want_fit = true;
}

static void an_orbit(sb2gui::App&, sb2gui::Editor& e, int) { e.orbit_view(2.6f, 0.0f); }

static const Stage kTour[] = {
    {26, "Control boxes in, a guaranteed tube out", st_boxes, nullptr},
    {30, "Drag a set: only its p+1 segments are re-evaluated", nullptr, an_drag},
    {26, "Control points as zonotopes: oriented, sheared sets", st_zono, nullptr},
    {26, "Generators are free vectors: the sets turn with them", nullptr, an_turn},
    {20, "Real control points: thin zonotopes enclosing the curve", st_curve, nullptr},
    {22, "Rational weights pull the evaluated points onto a control point", st_rational, nullptr},
    {34, "3D: the same symbolic basis, on (x, y, z)", st_3d, an_orbit},
    {34, "3D zonotopes: exact silhouettes, true facet meshes", st_3dz, an_orbit},
};

int main(int argc, char** argv)
{
    const char* shot_path = nullptr;
    const char* tour_dir = nullptr;
    int opt_dim = 2, opt_p = -1, opt_n = -1, opt_d = -1, opt_gen = -1;
    const char* opt_ps = nullptr;
    const char* opt_cs = nullptr;
    const char* opt_ct = nullptr;
    const char* opt_banner = nullptr;
    double opt_yaw = 1e9, opt_pitch = 1e9, opt_size = 0.0;

    for (int i = 1; i < argc; ++i) {
        const char* k = argv[i];
        const bool has = (i + 1 < argc);
#define ARG(name, var) if (std::strcmp(k, name) == 0 && has) { var = argv[++i]; continue; }
#define ARGI(name, var) if (std::strcmp(k, name) == 0 && has) { var = std::atoi(argv[++i]); continue; }
#define ARGF(name, var) if (std::strcmp(k, name) == 0 && has) { var = std::atof(argv[++i]); continue; }
        ARG("--shot", shot_path)
        ARG("--tour", tour_dir)
        ARG("--ps", opt_ps)
        ARG("--cs", opt_cs)
        ARG("--ct", opt_ct)
        ARG("--banner", opt_banner)
        ARGI("--dim", opt_dim)
        ARGI("--p", opt_p)
        ARGI("--n", opt_n)
        ARGI("--d", opt_d)
        ARGI("--gen", opt_gen)
        ARGF("--yaw", opt_yaw)
        ARGF("--pitch", opt_pitch)
        ARGF("--size", opt_size)
#undef ARG
#undef ARGI
#undef ARGF
        std::fprintf(stderr, "unknown or incomplete option: %s\n", k);
        return 2;
    }

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
    sb2gui::Editor editor(app);

    if (shot_path) { // one-shot capture: apply the requested configuration
        if (opt_p > 0) app.p = opt_p;
        if (opt_n > 0) app.nCP = opt_n;
        if (opt_d > 0) app.d = opt_d;
        if (opt_gen > 0) app.nGen = opt_gen;
        if (opt_ct) app.ct = parse_ct(opt_ct);
        if (opt_ps) app.ps = parse_set(opt_ps);
        if (opt_cs) app.cs = parse_set(opt_cs);
        if (opt_ct && (app.ct == sb2l::CurveType::UNIFORM_RATIONAL ||
                       app.ct == sb2l::CurveType::CLAMPED_RATIONAL)) {
            // A visible weight, so the rational capture shows what it does.
            app.wnum.resize(app.nCP, 1);
            app.wden.resize(app.nCP, 1);
            if (app.nCP > 2) app.wnum[2] = 6;
        }
        app.rebuild();
        if (opt_dim == 3) app.set_dim(3);
        if (opt_size > 0.0) { // blow the control sets up, to show their shape
            for (sb2gui::ControlPoint& c : app.cps) {
                c.hx *= opt_size; c.hy *= opt_size; c.hz *= opt_size;
                for (sb2gui::Vec3& g : c.gens) g = opt_size * g;
            }
            app.reeval();
        }
        if (opt_yaw < 1e8 || opt_pitch < 1e8)
            editor.set_view(opt_yaw < 1e8 ? opt_yaw : 0.7, opt_pitch < 1e8 ? opt_pitch : 0.45);
        if (opt_banner) app.banner = opt_banner;
        app.want_fit = true;
    }

    const int n_stages = (int)(sizeof(kTour) / sizeof(kTour[0]));
    int frame = 0, stage = 0, stage_frame = 0, shot_no = 0;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (tour_dir) {
            if (stage >= n_stages) break;
            if (stage_frame == 0) {
                app.banner = kTour[stage].caption;
                if (kTour[stage].enter) kTour[stage].enter(app);
            }
            if (kTour[stage].step) kTour[stage].step(app, editor, stage_frame);
        }

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
        ++frame;
        if (shot_path && frame >= 3) {
            write_ppm(window, shot_path);
            std::printf("wrote %s\n", shot_path);
            break;
        }
        if (tour_dir && frame >= 3) {
            char path[512];
            std::snprintf(path, sizeof(path), "%s/f%04d.ppm", tour_dir, shot_no++);
            write_ppm(window, path);
            if (++stage_frame >= kTour[stage].frames) { stage_frame = 0; ++stage; }
        }
    }
    if (tour_dir) std::printf("wrote %d tour frames to %s\n", shot_no, tour_dir);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
