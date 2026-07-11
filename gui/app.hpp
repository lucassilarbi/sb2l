/**
 * Editor state: B-spline parameters, control points, and cached draw geometry.
 *
 * The SB2 object precomputes the (expensive) symbolic basis, so it is rebuilt
 * only when a structural parameter changes (rebuild()); moving a control point
 * only re-runs the cheap linear-combination eval (reeval()).
 */
#ifndef SB2L_GUI_APP_HPP_
#define SB2L_GUI_APP_HPP_

#include <sb2l.hpp>
#include <memory>
#include <string>
#include <vector>

#include "zonotope2d.hpp"

namespace sb2gui {

// One 2D control point. Interpretation depends on the parameter set:
//   R : center only.  IR : box center +/- (hx,hy).  Z : center + generators.
struct ControlPoint {
    double cx = 0.0, cy = 0.0;
    double hx = 0.1, hy = 0.1;         // half-widths (IR)
    std::vector<Vec2> gens;            // generators (Z)
};

struct Box {
    double x0, y0, x1, y1;
};

class App {
public:
    // Structural parameters (change => rebuild).
    int p = 3;
    int nCP = 6;
    int d = 20;
    int t = -1;
    sb2l::CurveType ct = sb2l::CurveType::CLAMPED_NONRATIONAL;
    sb2l::Form f = sb2l::Form::TAYLOR;
    sb2l::ParameterSet ps = sb2l::ParameterSet::IR;
    std::vector<int> wnum, wden;       // rational weights, per control point

    std::vector<ControlPoint> cps;

    // Cached draw geometry (world coordinates), refreshed by reeval().
    std::vector<std::vector<Vec2>> polylines; // R: one polyline per segment
    std::vector<Box> boxes;                    // IR: result boxes
    std::vector<std::vector<Vec2>> zonos;      // Z: result zonotope polygons

    std::string status;                        // last rebuild/eval message
    bool dirty_structural = true;

    App();

    // Rebuild the SB2 object (resizes/preserves control points and weights).
    void rebuild();
    // Re-run evaluation with the current control points (no rebuild).
    void reeval();

    int nSegments() const { return spline ? spline->get_nS() : 0; }

private:
    std::unique_ptr<sb2l::SB2> spline;

    void resize_control_points(); // keep cps/weights consistent with nCP/p
    void seed_default_scene();
};

} // namespace sb2gui

#endif // SB2L_GUI_APP_HPP_
