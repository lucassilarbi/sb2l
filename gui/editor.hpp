/**
 * ImGui editor: parameter panel + interactive canvas (draw + hit-testing/drag).
 *
 * Two canvases share the window: the 2D one (pan/zoom, world x-y) and the 3D
 * one (orbit camera), selected by App::dim. Both edit the same control points
 * through the same Selection kinds; the 3D drags happen in the view-parallel
 * plane through the grabbed feature (its view depth is frozen at click time).
 */
#ifndef SB2L_GUI_EDITOR_HPP_
#define SB2L_GUI_EDITOR_HPP_

#include "app.hpp"
#include "camera3d.hpp"
#include "canvas.hpp"

#include <vector>

namespace sb2gui {

class Editor {
public:
    explicit Editor(App& app) : app_(app) {}

    void draw_controls_window();
    void draw_canvas_window();

private:
    enum class Kind { None, Point, BoxBody, BoxCorner, ZonoBody, ZonoGen };
    struct Selection {
        Kind kind = Kind::None;
        int index = -1;
        int sx = 0, sy = 0, sz = 0; // grabbed-corner signs (+/-1); ZonoGen uses sx only
        int gen = -1;               // grabbed generator (ZonoGen)
        double depth = 0.0;         // view depth of the grab (3D drag plane)
    };

    App& app_;
    Canvas canvas_;   // 2D view
    Camera3D cam_;    // 3D view
    Selection sel_;
    bool canvas_ready_ = false;
    float panel_w_ = 0.0f; // controls-window width, sized to its labels each frame

    // Scratch buffer for the polygon being pushed to ImGui, reused across shapes
    // and across frames so that drawing N zonotopes is not N allocations a frame.
    mutable std::vector<ImVec2> pts_;
    // Scratch for projected zonotope generators (3D silhouette path).
    mutable std::vector<Vec2> pgens_;

    // ---- 2D canvas ----
    void fit_view();
    void handle_input();
    Selection hit_test(ImVec2 mouse) const;
    void drag(const Selection& s);
    void draw_scene(ImDrawList* dl) const;
    // World polygon -> screen polygon in pts_, dropping sub-pixel edges.
    void to_screen_poly(const std::vector<Vec2>& world) const;

    // ---- 3D canvas ----
    void fit_view3();
    void handle_input3();
    Selection hit_test3(ImVec2 mouse) const;
    void drag3(const Selection& s);
    void draw_scene3(ImDrawList* dl) const;
    // Exact screen silhouette of a zonotope (center + 3D generators) through
    // the orthographic camera: projected generators -> 2D boundary walk.
    // Fills pts_ (sub-pixel edges dropped). Returns false for a sub-polygon
    // result (point/segment on screen).
    bool silhouette3(const Vec3& c, const std::vector<Vec3>& gens) const;
};

} // namespace sb2gui

#endif // SB2L_GUI_EDITOR_HPP_
