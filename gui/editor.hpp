/**
 * ImGui editor: parameter panel + interactive canvas (draw + hit-testing/drag).
 */
#ifndef SB2L_GUI_EDITOR_HPP_
#define SB2L_GUI_EDITOR_HPP_

#include "app.hpp"
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
        int sx = 0, sy = 0; // grabbed-corner signs (+/-1); ZonoGen uses sx only
        int gen = -1;       // grabbed generator (ZonoGen)
    };

    App& app_;
    Canvas canvas_;
    Selection sel_;
    bool canvas_ready_ = false;

    // Scratch buffer for the polygon being pushed to ImGui, reused across shapes
    // and across frames so that drawing N zonotopes is not N allocations a frame.
    mutable std::vector<ImVec2> pts_;

    void fit_view();
    void handle_input();
    Selection hit_test(ImVec2 mouse) const;
    void drag(const Selection& s);
    void draw_scene(ImDrawList* dl) const;
    // World polygon -> screen polygon in pts_, dropping sub-pixel edges.
    void to_screen_poly(const std::vector<Vec2>& world) const;
};

} // namespace sb2gui

#endif // SB2L_GUI_EDITOR_HPP_
