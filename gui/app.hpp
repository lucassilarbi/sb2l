/**
 * Editor state: B-spline parameters, control points, and cached draw geometry.
 *
 * The SB2 object precomputes the (expensive) symbolic basis, so it is rebuilt
 * only when a structural parameter changes (rebuild()); moving a control point
 * only re-runs the cheap linear-combination eval (reeval()).
 *
 * Dragging goes one step further: a control point only belongs to p+1 segments,
 * so update_control_point() re-evaluates and redraws those alone (see
 * SB2::impacted_segments) instead of the whole curve.
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
// A zonotope's shape lives entirely in its generators: each one is an arbitrary
// 2D vector, so both the extent *and* the orientation are edited by moving the
// generator tips (c +/- g_i), never through an axis-aligned box.
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
    int nGen = 2;                      // generators per control point (Z)
    sb2l::CurveType ct = sb2l::CurveType::CLAMPED_NONRATIONAL;
    sb2l::Form f = sb2l::Form::TAYLOR;
    sb2l::ParameterSet ps = sb2l::ParameterSet::IR;
    std::vector<int> wnum, wden;       // rational weights, per control point

    std::vector<ControlPoint> cps;

    // Cached draw geometry (world coordinates), refreshed by reeval().
    // Segment s owns a fixed slice of each container: the d elements at [s*d, s*d + d),
    // which is what lets a segment be redrawn in place (the R polyline holds one extra
    // sample at the very end, closing the curve).
    std::vector<std::vector<Vec2>> polylines; // R: the curve, as a single polyline
    std::vector<Box> boxes;                    // IR: result boxes
    std::vector<std::vector<Vec2>> zonos;      // Z: result zonotope polygons

    std::string status;                        // last rebuild/eval message
    bool dirty_structural = true;

    App();

    // Rebuild the SB2 object (resizes/preserves control points and weights).
    void rebuild();
    // Re-run evaluation with the current control points (no rebuild).
    void reeval();
    // Same, but only for the segments the control point i belongs to: this is the
    // drag path, and it costs p+1 segments instead of nCP - p.
    void update_control_point(int i);

    int nSegments() const { return spline ? spline->get_nS() : 0; }

private:
    std::unique_ptr<sb2l::SB2> spline;

    // Raw evaluation kept per segment, so that a moved control point can be pushed
    // through the impacted segments only. Exactly one of them is filled, per ps.
    std::vector<std::vector<std::vector<double>>> epoints;
    std::vector<std::vector<ibex::IntervalVector>> eboxes;
    std::vector<std::vector<ibex::Affine2Vector>> ezonos;

    void resize_control_points(); // keep cps/weights consistent with nCP/p
    void seed_default_scene();

    // Control points of the current scene, in the layout expected by SB2::eval_*.
    std::vector<std::vector<double>> control_points_R() const;
    std::vector<ibex::IntervalVector> control_points_IR() const;
    std::vector<ibex::Affine2Vector> control_points_Z() const;
    void alloc_geometry();                  // size the draw containers for the current spline
    void refresh_geometry(int s0, int s1);  // rebuild the geometry of segments [s0, s1] from the raw eval
};

} // namespace sb2gui

#endif // SB2L_GUI_APP_HPP_
