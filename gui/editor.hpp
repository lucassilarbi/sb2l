/**
 * ImGui editor: parameter panel + interactive canvas (draw + hit-testing/drag).
 */
#ifndef SB2L_GUI_EDITOR_HPP_
#define SB2L_GUI_EDITOR_HPP_

#include "app.hpp"
#include "canvas.hpp"

namespace sb2gui {

class Editor {
public:
    explicit Editor(App& app) : app_(app) {}

    void draw_controls_window();
    void draw_canvas_window();

private:
    enum class Kind { None, Point, BoxBody, BoxCorner, ZonoBody, ZonoCorner };
    struct Selection {
        Kind kind = Kind::None;
        int index = -1;
        int sx = 0, sy = 0; // grabbed-corner signs (+/-1)
    };

    App& app_;
    Canvas canvas_;
    Selection sel_;
    bool canvas_ready_ = false;

    void fit_view();
    void handle_input();
    Selection hit_test(ImVec2 mouse) const;
    void drag(const Selection& s);
    void draw_scene(ImDrawList* dl) const;
};

} // namespace sb2gui

#endif // SB2L_GUI_EDITOR_HPP_
