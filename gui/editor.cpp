#include "editor.hpp"

#include <algorithm>
#include <cmath>

namespace sb2gui {

static const float kPanelWidth = 320.0f; // Controls window width

static void zono_extent(const ControlPoint& cp, double& ex, double& ey)
{
    ex = 0.0; ey = 0.0;
    for (const Vec2& g : cp.gens) { ex += std::fabs(g.x); ey += std::fabs(g.y); }
}

static float dist(ImVec2 a, ImVec2 b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// ---------------------------------------------------------------- controls ---

void Editor::draw_controls_window()
{
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(ImVec2(kPanelWidth, vp->WorkSize.y));
    ImGui::Begin("Controls", nullptr,
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
    App& a = app_;
    bool changed = false;

    changed |= ImGui::SliderInt("degree p", &a.p, 1, 8);
    changed |= ImGui::SliderInt("num control pts", &a.nCP, a.p + 1, 24);
    changed |= ImGui::SliderInt("subdivisions d", &a.d, 1, 60);
    changed |= ImGui::SliderInt("taylor t (-1=auto)", &a.t, -1, 11);

    const char* curve_types[] = {"UNIFORM_RATIONAL", "UNIFORM_NONRATIONAL",
                                 "CLAMPED_RATIONAL", "CLAMPED_NONRATIONAL"};
    const char* forms[] = {"NATURAL", "TAYLOR"};
    const char* psets[] = {"R (points)", "IR (boxes)", "Z (zonotopes)"};

    int ct = (int)a.ct, f = (int)a.f, ps = (int)a.ps;
    if (ImGui::Combo("curve type", &ct, curve_types, 4)) { a.ct = (sb2l::CurveType)ct; changed = true; }
    if (ImGui::Combo("form", &f, forms, 2)) { a.f = (sb2l::Form)f; changed = true; }
    if (ImGui::Combo("parameter set", &ps, psets, 3)) { a.ps = (sb2l::ParameterSet)ps; changed = true; }

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

    if (changed) a.rebuild();

    ImGui::Separator();
    if (ImGui::Button("Fit view")) fit_view();
    ImGui::TextWrapped("Drag control handles to edit. In box/zonotope modes, drag inside "
                       "to move and drag a corner to resize. Wheel zooms, drag empty space pans.");
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
    if (a.ps == sb2l::ParameterSet::IR) {
        for (int i = 0; i < (int)a.cps.size(); ++i)
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2)
                    if (dist(m, canvas_.to_screen(a.cps[i].cx + sx * a.cps[i].hx,
                                                  a.cps[i].cy + sy * a.cps[i].hy)) < r)
                        return {Kind::BoxCorner, i, sx, sy};
    } else if (a.ps == sb2l::ParameterSet::Z) {
        for (int i = 0; i < (int)a.cps.size(); ++i) {
            double ex, ey; zono_extent(a.cps[i], ex, ey);
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2)
                    if (dist(m, canvas_.to_screen(a.cps[i].cx + sx * ex,
                                                  a.cps[i].cy + sy * ey)) < r)
                        return {Kind::ZonoCorner, i, sx, sy};
        }
    }
    // Bodies / centers.
    double wx, wy; canvas_.to_world(m, wx, wy);
    for (int i = 0; i < (int)a.cps.size(); ++i) {
        const ControlPoint& c = a.cps[i];
        if (a.ps == sb2l::ParameterSet::R) {
            if (dist(m, canvas_.to_screen(c.cx, c.cy)) < r + 2) return {Kind::Point, i, 0, 0};
        } else if (a.ps == sb2l::ParameterSet::IR) {
            if (wx >= c.cx - c.hx && wx <= c.cx + c.hx && wy >= c.cy - c.hy && wy <= c.cy + c.hy)
                return {Kind::BoxBody, i, 0, 0};
        } else {
            double ex, ey; zono_extent(c, ex, ey);
            if (wx >= c.cx - ex && wx <= c.cx + ex && wy >= c.cy - ey && wy <= c.cy + ey)
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
    case Kind::ZonoCorner: {
        double ex, ey; zono_extent(c, ex, ey);
        double opx = c.cx - s.sx * ex, opy = c.cy - s.sy * ey;
        double nex = std::fabs(mx - opx) * 0.5, ney = std::fabs(my - opy) * 0.5;
        double fx = ex > 1e-9 ? nex / ex : 1.0;
        double fy = ey > 1e-9 ? ney / ey : 1.0;
        for (Vec2& g : c.gens) { g.x *= fx; g.y *= fy; }
        c.cx = 0.5 * (mx + opx); c.cy = 0.5 * (my + opy);
        break;
    }
    default: return;
    }
    a.reeval();
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
        std::vector<ImVec2> pts;
        pts.reserve(z.size());
        for (const Vec2& q : z) pts.push_back(canvas_.to_screen(q.x, q.y));
        if (pts.size() >= 3) dl->AddConvexPolyFilled(pts.data(), (int)pts.size(), col_zono_f);
        dl->AddPolyline(pts.data(), (int)pts.size(), col_zono, ImDrawFlags_Closed, 1.5f);
    }

    // Control-point handles (on top).
    for (const ControlPoint& c : a.cps) {
        ImVec2 ctr = canvas_.to_screen(c.cx, c.cy);
        if (a.ps == sb2l::ParameterSet::IR) {
            ImVec2 tl = canvas_.to_screen(c.cx - c.hx, c.cy + c.hy);
            ImVec2 br = canvas_.to_screen(c.cx + c.hx, c.cy - c.hy);
            dl->AddRect(tl, br, col_handle);
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2) {
                    ImVec2 p = canvas_.to_screen(c.cx + sx * c.hx, c.cy + sy * c.hy);
                    dl->AddRectFilled(ImVec2(p.x - 3, p.y - 3), ImVec2(p.x + 3, p.y + 3), col_handle);
                }
        } else if (a.ps == sb2l::ParameterSet::Z) {
            std::vector<Vec2> b = zonotope_boundary({c.cx, c.cy}, c.gens);
            if (b.size() >= 2) {
                std::vector<ImVec2> pts;
                for (const Vec2& q : b) pts.push_back(canvas_.to_screen(q.x, q.y));
                dl->AddPolyline(pts.data(), (int)pts.size(), col_handle, ImDrawFlags_Closed, 1.0f);
            }
            double ex, ey; zono_extent(c, ex, ey);
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2) {
                    ImVec2 p = canvas_.to_screen(c.cx + sx * ex, c.cy + sy * ey);
                    dl->AddRectFilled(ImVec2(p.x - 3, p.y - 3), ImVec2(p.x + 3, p.y + 3), col_handle);
                }
        }
        dl->AddCircleFilled(ctr, 4.0f, col_handle);
    }
}

void Editor::draw_canvas_window()
{
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + kPanelWidth, vp->WorkPos.y));
    ImGui::SetNextWindowSize(ImVec2(vp->WorkSize.x - kPanelWidth, vp->WorkSize.y));
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
