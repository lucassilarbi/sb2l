#include "editor.hpp"

#include <algorithm>
#include <cmath>

namespace sb2gui {

// Every label of the controls panel. ImGui draws a label to the right of its
// widget and clips it at the window edge, so the panel is only ever as usable as
// the room left for the longest of them: both widths are derived from it below.
static const char* kLabels[] = {"degree p", "num control pts", "subdivisions d",
                                "taylor t (-1=auto)", "curve type", "form",
                                "curve parameter", "control points", "generators / point"};

// Width the labels need, and the panel width that leaves the widgets a usable
// share next to them (font-relative, so a larger font still fits).
static float label_width()
{
    float w = 0.0f;
    for (const char* l : kLabels) w = std::max(w, ImGui::CalcTextSize(l).x);
    return w + ImGui::GetStyle().ItemInnerSpacing.x;
}
static float panel_width()
{
    return label_width() + ImGui::GetFontSize() * 14.0f; // widget column: fits "CLAMPED_NONRATIONAL"
}

// AABB half-widths of a zonotope. For view fitting only: this loses the
// orientation, so it must not drive editing.
static void zono_extent(const ControlPoint& cp, double& ex, double& ey)
{
    ex = 0.0; ey = 0.0;
    for (const Vec2& g : cp.gens) { ex += std::fabs(g.x); ey += std::fabs(g.y); }
}

static float dist(ImVec2 a, ImVec2 b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// Inside test for a CCW convex polygon (the zonotope boundary is already CCW).
static bool inside_ccw(const std::vector<Vec2>& poly, double x, double y)
{
    const size_t n = poly.size();
    if (n < 3) return false;
    for (size_t i = 0; i < n; ++i) {
        const Vec2& a = poly[i];
        const Vec2& b = poly[(i + 1) % n];
        if ((b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x) < 0.0) return false;
    }
    return true;
}

// Project a world polygon to screen, skipping a vertex that lands within half a
// pixel of the one before it. A result zonotope carries one vertex pair per noise
// symbol -- hundreds of them once the degree and generator count go up -- and at
// screen scale nearly all of those edges are far shorter than a pixel: they cost
// ImGui vertices every frame without changing a single rendered pixel. Dropping
// them keeps the outline and stays convex, so the fill path still applies.
void Editor::to_screen_poly(const std::vector<Vec2>& world) const
{
    const float min_d2 = 0.25f; // (0.5 px)^2
    pts_.clear();
    pts_.reserve(world.size());
    for (const Vec2& q : world) {
        const ImVec2 p = canvas_.to_screen(q.x, q.y);
        if (!pts_.empty()) {
            const float dx = p.x - pts_.back().x, dy = p.y - pts_.back().y;
            if (dx * dx + dy * dy < min_d2) continue;
        }
        pts_.push_back(p);
    }
    // The polygon is closed, so a last vertex sitting on the first is redundant too.
    if (pts_.size() >= 2) {
        const float dx = pts_.back().x - pts_.front().x, dy = pts_.back().y - pts_.front().y;
        if (dx * dx + dy * dy < min_d2) pts_.pop_back();
    }
}

// ---------------------------------------------------------------- controls ---

void Editor::draw_controls_window()
{
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    panel_w_ = panel_width();
    ImGui::SetNextWindowSize(ImVec2(panel_w_, vp->WorkSize.y));
    ImGui::Begin("Controls", nullptr,
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
    App& a = app_;
    bool changed = false;    // structural: needs a new symbolic basis (rebuild)
    bool reeval = false;     // same basis, other control points (eval only)

    ImGui::PushItemWidth(-label_width()); // leave every label its full width

    changed |= ImGui::SliderInt("degree p", &a.p, 1, 8);
    changed |= ImGui::SliderInt("num control pts", &a.nCP, a.p + 1, 24);
    changed |= ImGui::SliderInt("subdivisions d", &a.d, 1, 60);
    changed |= ImGui::SliderInt("taylor t (-1=auto)", &a.t, -1, 11);

    const char* curve_types[] = {"UNIFORM_RATIONAL", "UNIFORM_NONRATIONAL",
                                 "CLAMPED_RATIONAL", "CLAMPED_NONRATIONAL"};
    const char* forms[] = {"NATURAL", "TAYLOR"};
    // The set the parameter u is taken in: it is what the basis is evaluated over.
    const char* psets[] = {"R (real u)", "IR (interval u)", "Z (affine u)"};
    // The set the control points are taken in, free of the one above.
    const char* csets[] = {"R (points)", "IR (boxes)", "Z (zonotopes)"};

    int ct = (int)a.ct, f = (int)a.f, ps = (int)a.ps, cs = (int)a.cs;
    if (ImGui::Combo("curve type", &ct, curve_types, 4)) { a.ct = (sb2l::CurveType)ct; changed = true; }
    if (ImGui::Combo("form", &f, forms, 2)) { a.f = (sb2l::Form)f; changed = true; }
    if (ImGui::Combo("curve parameter", &ps, psets, 3)) { a.ps = (sb2l::ParameterSet)ps; changed = true; }
    // Switching the control points reuses the basis as it is: no rebuild, just a
    // fresh eval_point / eval_box / eval_zonotope over the same B-spline.
    if (ImGui::Combo("control points", &cs, csets, 3)) { a.cs = (sb2l::ParameterSet)cs; reeval = true; }

    if (a.cs == sb2l::ParameterSet::Z)
        changed |= ImGui::SliderInt("generators / point", &a.nGen, 1, 6);

    bool rational = (a.ct == sb2l::CurveType::UNIFORM_RATIONAL ||
                     a.ct == sb2l::CurveType::CLAMPED_RATIONAL);
    if (rational) {
        ImGui::SeparatorText("rational weights (num / den)");
        for (int i = 0; i < a.nCP && i < (int)a.wnum.size(); ++i) {
            ImGui::PushID(i);
            ImGui::SetNextItemWidth(80);
            changed |= ImGui::InputInt("##num", &a.wnum[i]);
            if (a.wnum[i] < 1) a.wnum[i] = 1; // weights must stay positive
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80);
            changed |= ImGui::InputInt("##den", &a.wden[i]);
            if (a.wden[i] < 1) a.wden[i] = 1;
            ImGui::SameLine();
            ImGui::Text("w[%d]", i);
            ImGui::PopID();
        }
    }

    ImGui::PopItemWidth();

    if (changed) a.rebuild();
    else if (reeval) a.reeval();

    ImGui::Separator();
    if (ImGui::Button("Fit view")) fit_view();
    ImGui::TextWrapped("The curve parameter and the control points are chosen apart: any of the 9 "
                       "pairs is a B-spline (real points over an interval parameter give the tube "
                       "of the curve, and so on).");
    ImGui::TextWrapped("Drag control handles to edit. Box mode: drag inside to move, drag a corner "
                       "to resize. Zonotope mode: drag inside to move, drag a yellow generator tip "
                       "to set its direction and length (that is how you rotate/shear the shape). "
                       "Wheel zooms, drag empty space pans.");
    ImGui::Text("segments: %d   status: %s", a.nSegments(), a.status.c_str());
    ImGui::End();
}

// ------------------------------------------------------------------ canvas ---

void Editor::fit_view()
{
    App& a = app_;
    if (a.cps.empty()) return;
    double x0 = a.cps[0].cx, x1 = x0, y0 = a.cps[0].cy, y1 = y0;
    for (const ControlPoint& c : a.cps) {
        double ex = c.hx, ey = c.hy, zx, zy;
        zono_extent(c, zx, zy);
        if (zx > ex) ex = zx;
        if (zy > ey) ey = zy;
        x0 = std::min(x0, c.cx - ex); x1 = std::max(x1, c.cx + ex);
        y0 = std::min(y0, c.cy - ey); y1 = std::max(y1, c.cy + ey);
    }
    double w = std::max(1e-3, x1 - x0), h = std::max(1e-3, y1 - y0);
    float sx = canvas_.size.x * 0.85f / (float)w;
    float sy = canvas_.size.y * 0.85f / (float)h;
    canvas_.scale = std::min(sx, sy);
    double cx = 0.5 * (x0 + x1), cy = 0.5 * (y0 + y1);
    canvas_.origin.x = canvas_.tl.x + canvas_.size.x * 0.5f - (float)cx * canvas_.scale;
    canvas_.origin.y = canvas_.tl.y + canvas_.size.y * 0.5f + (float)cy * canvas_.scale;
}

Editor::Selection Editor::hit_test(ImVec2 m) const
{
    const App& a = app_;
    const float r = 7.0f;
    // Corner handles first (they sit on top and are small).
    if (a.cs == sb2l::ParameterSet::IR) {
        for (int i = 0; i < (int)a.cps.size(); ++i)
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2)
                    if (dist(m, canvas_.to_screen(a.cps[i].cx + sx * a.cps[i].hx,
                                                  a.cps[i].cy + sy * a.cps[i].hy)) < r)
                        return {Kind::BoxCorner, i, sx, sy};
    } else if (a.cs == sb2l::ParameterSet::Z) {
        // Generator tips: each one is free in direction and length, which is
        // what lets a zonotope be reoriented rather than merely rescaled.
        for (int i = 0; i < (int)a.cps.size(); ++i) {
            const ControlPoint& c = a.cps[i];
            for (int k = 0; k < (int)c.gens.size(); ++k)
                for (int s = -1; s <= 1; s += 2)
                    if (dist(m, canvas_.to_screen(c.cx + s * c.gens[k].x,
                                                  c.cy + s * c.gens[k].y)) < r)
                        return {Kind::ZonoGen, i, s, 0, k};
        }
    }
    // Bodies / centers.
    double wx, wy; canvas_.to_world(m, wx, wy);
    for (int i = 0; i < (int)a.cps.size(); ++i) {
        const ControlPoint& c = a.cps[i];
        if (a.cs == sb2l::ParameterSet::R) {
            if (dist(m, canvas_.to_screen(c.cx, c.cy)) < r + 2) return {Kind::Point, i, 0, 0};
        } else if (a.cs == sb2l::ParameterSet::IR) {
            if (wx >= c.cx - c.hx && wx <= c.cx + c.hx && wy >= c.cy - c.hy && wy <= c.cy + c.hy)
                return {Kind::BoxBody, i, 0, 0};
        } else {
            // Test the true polygon, not its bounding box: a rotated zonotope
            // is much smaller than its AABB and must not grab clicks outside it.
            if (inside_ccw(zonotope_boundary({c.cx, c.cy}, c.gens), wx, wy) ||
                dist(m, canvas_.to_screen(c.cx, c.cy)) < r + 2)
                return {Kind::ZonoBody, i, 0, 0};
        }
    }
    return {};
}

void Editor::drag(const Selection& s)
{
    App& a = app_;
    ImGuiIO& io = ImGui::GetIO();
    ControlPoint& c = a.cps[s.index];
    double mx, my; canvas_.to_world(io.MousePos, mx, my);
    double ddx = io.MouseDelta.x / canvas_.scale;
    double ddy = -io.MouseDelta.y / canvas_.scale;

    switch (s.kind) {
    case Kind::Point:
        c.cx = mx; c.cy = my;
        break;
    case Kind::BoxBody:
    case Kind::ZonoBody:
        c.cx += ddx; c.cy += ddy;
        break;
    case Kind::BoxCorner: {
        double opx = c.cx - s.sx * c.hx, opy = c.cy - s.sy * c.hy;
        c.hx = std::max(1e-4, std::fabs(mx - opx) * 0.5);
        c.hy = std::max(1e-4, std::fabs(my - opy) * 0.5);
        c.cx = 0.5 * (mx + opx); c.cy = 0.5 * (my + opy);
        break;
    }
    case Kind::ZonoGen: {
        // The grabbed tip is c + s.sx * g, so the mouse fixes the generator
        // outright: direction and length both follow the cursor.
        Vec2& g = c.gens[s.gen];
        double gx = s.sx * (mx - c.cx), gy = s.sx * (my - c.cy);
        double len = std::sqrt(gx * gx + gy * gy);
        if (len < 1e-4) { // never let a generator collapse to an unrecoverable zero
            double plen = std::sqrt(g.x * g.x + g.y * g.y);
            if (plen > 1e-12) { gx = 1e-4 * g.x / plen; gy = 1e-4 * g.y / plen; }
            else { gx = 1e-4; gy = 0.0; }
        }
        g.x = gx; g.y = gy;
        break;
    }
    default: return;
    }
    // Only the dragged control point moved, so only its p+1 segments need it.
    a.update_control_point(s.index);
}

void Editor::handle_input()
{
    if (!ImGui::IsItemHovered() && !ImGui::IsItemActive()) return;
    ImGuiIO& io = ImGui::GetIO();
    if (ImGui::IsItemHovered() && io.MouseWheel != 0.0f)
        canvas_.zoom_at(io.MousePos, std::pow(1.1f, io.MouseWheel));
    if (ImGui::IsItemActive()) {
        if (ImGui::IsMouseClicked(0)) sel_ = hit_test(io.MousePos);
        if (ImGui::IsMouseDown(0)) {
            if (sel_.kind == Kind::None) canvas_.pan(io.MouseDelta);
            else drag(sel_);
        }
    }
    if (ImGui::IsMouseReleased(0)) sel_ = Selection{};
}

void Editor::draw_scene(ImDrawList* dl) const
{
    const App& a = app_;
    const ImU32 col_axis = IM_COL32(90, 90, 90, 255);
    const ImU32 col_ctrl = IM_COL32(150, 150, 150, 180);
    const ImU32 col_curve = IM_COL32(80, 200, 120, 255);
    const ImU32 col_box = IM_COL32(230, 170, 60, 255);
    const ImU32 col_box_f = IM_COL32(230, 170, 60, 40);
    const ImU32 col_zono = IM_COL32(230, 90, 90, 255);
    const ImU32 col_zono_f = IM_COL32(230, 90, 90, 45);
    const ImU32 col_handle = IM_COL32(80, 150, 240, 255);
    const ImU32 col_gen = IM_COL32(240, 200, 90, 255);

    // Axes.
    ImVec2 o = canvas_.to_screen(0, 0);
    dl->AddLine(ImVec2(canvas_.tl.x, o.y), ImVec2(canvas_.tl.x + canvas_.size.x, o.y), col_axis);
    dl->AddLine(ImVec2(o.x, canvas_.tl.y), ImVec2(o.x, canvas_.tl.y + canvas_.size.y), col_axis);

    // Control polygon.
    if (a.cps.size() >= 2) {
        std::vector<ImVec2> poly;
        poly.reserve(a.cps.size());
        for (const ControlPoint& c : a.cps) poly.push_back(canvas_.to_screen(c.cx, c.cy));
        dl->AddPolyline(poly.data(), (int)poly.size(), col_ctrl, 0, 1.0f);
    }

    // Result geometry.
    for (const std::vector<Vec2>& line : a.polylines) {
        if (line.size() < 2) continue;
        std::vector<ImVec2> pts;
        pts.reserve(line.size());
        for (const Vec2& q : line) pts.push_back(canvas_.to_screen(q.x, q.y));
        dl->AddPolyline(pts.data(), (int)pts.size(), col_curve, 0, 2.0f);
    }
    for (const Box& b : a.boxes) {
        ImVec2 p0 = canvas_.to_screen(b.x0, b.y1); // top-left in screen
        ImVec2 p1 = canvas_.to_screen(b.x1, b.y0);
        dl->AddRectFilled(p0, p1, col_box_f);
        dl->AddRect(p0, p1, col_box);
    }
    for (const std::vector<Vec2>& z : a.zonos) {
        if (z.size() < 2) continue;
        to_screen_poly(z);
        if (pts_.size() >= 3) {
            dl->AddConvexPolyFilled(pts_.data(), (int)pts_.size(), col_zono_f);
            dl->AddPolyline(pts_.data(), (int)pts_.size(), col_zono, ImDrawFlags_Closed, 1.5f);
        } else if (!pts_.empty()) {
            // Whole zonotope is under a pixel wide: a dot is all it can render as.
            const ImVec2 p = pts_.front();
            dl->AddRectFilled(ImVec2(p.x, p.y), ImVec2(p.x + 1.0f, p.y + 1.0f), col_zono);
        }
    }

    // Control-point handles (on top).
    for (const ControlPoint& c : a.cps) {
        ImVec2 ctr = canvas_.to_screen(c.cx, c.cy);
        if (a.cs == sb2l::ParameterSet::IR) {
            ImVec2 tl = canvas_.to_screen(c.cx - c.hx, c.cy + c.hy);
            ImVec2 br = canvas_.to_screen(c.cx + c.hx, c.cy - c.hy);
            dl->AddRect(tl, br, col_handle);
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2) {
                    ImVec2 p = canvas_.to_screen(c.cx + sx * c.hx, c.cy + sy * c.hy);
                    dl->AddRectFilled(ImVec2(p.x - 3, p.y - 3), ImVec2(p.x + 3, p.y + 3), col_handle);
                }
        } else if (a.cs == sb2l::ParameterSet::Z) {
            const std::vector<Vec2> b = zonotope_boundary({c.cx, c.cy}, c.gens);
            if (b.size() >= 2) {
                to_screen_poly(b);
                if (pts_.size() >= 2)
                    dl->AddPolyline(pts_.data(), (int)pts_.size(), col_handle, ImDrawFlags_Closed, 1.0f);
            }
            // One spoke per generator, grabbable at either tip.
            for (const Vec2& g : c.gens) {
                ImVec2 tp = canvas_.to_screen(c.cx + g.x, c.cy + g.y);
                ImVec2 tm = canvas_.to_screen(c.cx - g.x, c.cy - g.y);
                dl->AddLine(tm, tp, col_gen, 1.0f);
                dl->AddCircleFilled(tp, 4.0f, col_gen);
                dl->AddCircleFilled(tm, 4.0f, col_gen);
            }
        }
        dl->AddCircleFilled(ctr, 4.0f, col_handle);
    }
}

void Editor::draw_canvas_window()
{
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    // The controls window is laid out first this frame, so its width is known.
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + panel_w_, vp->WorkPos.y));
    ImGui::SetNextWindowSize(ImVec2(vp->WorkSize.x - panel_w_, vp->WorkSize.y));
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Canvas", nullptr,
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoBringToFrontOnFocus);
    ImGui::PopStyleVar();
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (avail.x < 32) avail.x = 32;
    if (avail.y < 32) avail.y = 32;
    ImVec2 p0 = ImGui::GetCursorScreenPos();
    canvas_.tl = p0;
    canvas_.size = avail;
    if (!canvas_ready_) {
        canvas_.origin = ImVec2(p0.x + avail.x * 0.5f, p0.y + avail.y * 0.5f);
        fit_view();
        canvas_ready_ = true;
    }
    ImGui::InvisibleButton("canvas", avail,
                           ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(p0, ImVec2(p0.x + avail.x, p0.y + avail.y), true);
    dl->AddRectFilled(p0, ImVec2(p0.x + avail.x, p0.y + avail.y), IM_COL32(28, 28, 30, 255));
    handle_input();
    draw_scene(dl);
    dl->PopClipRect();
    ImGui::End();
}

} // namespace sb2gui
