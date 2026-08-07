/**
 * Orbit camera and world <-> screen transform for the 3D canvas.
 * Header-only. Depends only on ImGui's ImVec2 (and zonotope3d.hpp for Vec3).
 *
 * The projection is orthographic, on purpose: the linear map "world -> screen
 * plane" then sends a zonotope to an exact 2D zonotope (project the center and
 * each generator, run the O(m log m) walk of zonotope2d.hpp), so the outline
 * drawn is the exact silhouette of the set -- no perspective would allow that.
 * What reaches the screen is that silhouette rounded to the pixel grid, the
 * vertices closer than half a pixel being dropped along the way.
 */
#ifndef SB2L_GUI_CAMERA3D_HPP_
#define SB2L_GUI_CAMERA3D_HPP_

#include "imgui.h"
#include "zonotope3d.hpp"

#include <cmath>

namespace sb2gui {

struct Camera3D {
    // Orbit state: the camera looks at `target` along -eye_dir(yaw, pitch),
    // world z up. scale is pixels per world unit (the orthographic "zoom").
    double yaw = 0.7;
    double pitch = 0.45;
    double scale = 60.0;
    Vec3 target;

    ImVec2 tl, size; // canvas rect on screen (top-left, size)

    // View basis, recomputed once per frame before any projection.
    Vec3 fwd, right, up;

    void recompute()
    {
        const double cp = std::cos(pitch), sp = std::sin(pitch);
        const double cy = std::cos(yaw), sy = std::sin(yaw);
        const Vec3 eye{cp * cy, cp * sy, sp}; // target -> camera direction
        fwd = -1.0 * eye;
        // pitch is clamped short of +/- pi/2, so fwd never aligns with world z.
        Vec3 r = cross(fwd, Vec3{0.0, 0.0, 1.0});
        right = (1.0 / norm(r)) * r;
        up = cross(right, fwd);
    }

    ImVec2 center() const { return ImVec2(tl.x + size.x * 0.5f, tl.y + size.y * 0.5f); }

    // Screen y points down, camera up points up.
    ImVec2 to_screen(const Vec3& w) const
    {
        const Vec3 v = w - target;
        const ImVec2 c = center();
        return ImVec2(c.x + (float)(dot(v, right) * scale), c.y - (float)(dot(v, up) * scale));
    }

    // View depth: grows away from the camera (painter's sort key).
    double depth(const Vec3& w) const { return dot(w - target, fwd); }

    // Linear part of the projection, in *world* units on the view plane
    // (right/up components). This is what maps a zonotope generator.
    void project_dir(const Vec3& v, double& px, double& py) const
    {
        px = dot(v, right);
        py = dot(v, up);
    }

    // Inverse of to_screen at a fixed view depth: the world point that
    // projects on s and sits `d` deep. Dragging a feature keeps its own
    // depth, so it moves in the view-parallel plane through itself.
    Vec3 unproject(ImVec2 s, double d) const
    {
        const ImVec2 c = center();
        const double px = (s.x - c.x) / scale;
        const double py = (c.y - s.y) / scale;
        return target + px * right + py * up + d * fwd;
    }

    void orbit(ImVec2 delta)
    {
        yaw -= delta.x * 0.010;
        pitch += delta.y * 0.010;
        const double lim = 1.55; // just short of the poles (keeps `right` defined)
        if (pitch > lim) pitch = lim;
        if (pitch < -lim) pitch = -lim;
    }

    void pan(ImVec2 delta)
    {
        target = target - (delta.x / scale) * right + (delta.y / scale) * up;
    }

    // Zoom about a fixed screen point (keeps that point under the cursor):
    // the world point at depth 0 under the cursor must reproject in place.
    void zoom_at(ImVec2 s, float factor)
    {
        const Vec3 w = unproject(s, 0.0);
        scale *= factor;
        if (scale < 1e-3) scale = 1e-3;
        const ImVec2 s2 = to_screen(w);
        target = target + ((s2.x - s.x) / scale) * right - ((s2.y - s.y) / scale) * up;
    }
};

} // namespace sb2gui

#endif // SB2L_GUI_CAMERA3D_HPP_
