/**
 * World <-> screen transform for the 2D canvas: pan, zoom, y-up.
 * Header-only. Depends only on ImGui's ImVec2.
 */
#ifndef SB2L_GUI_CANVAS_HPP_
#define SB2L_GUI_CANVAS_HPP_

#include "imgui.h"

namespace sb2gui {

struct Canvas {
    ImVec2 origin;       // screen pixel of world (0,0)
    float scale = 60.0f; // pixels per world unit
    ImVec2 tl, size;     // canvas rect on screen (top-left, size)

    // World y points up, screen y points down.
    ImVec2 to_screen(double wx, double wy) const
    {
        return ImVec2(origin.x + (float)wx * scale, origin.y - (float)wy * scale);
    }
    void to_world(ImVec2 s, double& wx, double& wy) const
    {
        wx = (s.x - origin.x) / scale;
        wy = (origin.y - s.y) / scale;
    }
    // Zoom about a fixed screen point (keeps that point under the cursor).
    void zoom_at(ImVec2 s, float factor)
    {
        double wx, wy;
        to_world(s, wx, wy);
        scale *= factor;
        if (scale < 1e-3f) scale = 1e-3f;
        origin.x = s.x - (float)wx * scale;
        origin.y = s.y + (float)wy * scale;
    }
    void pan(ImVec2 delta)
    {
        origin.x += delta.x;
        origin.y += delta.y;
    }
};

} // namespace sb2gui

#endif // SB2L_GUI_CANVAS_HPP_
