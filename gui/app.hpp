/**
 * Editor state: B-spline parameters, control points, and cached draw geometry.
 *
 * The SB2 object precomputes the (expensive) symbolic basis, so it is rebuilt
 * only when a structural parameter changes (rebuild()); moving a control point
 * only re-runs the cheap linear-combination eval (reeval()).
 *
 * SB2 is dimension-agnostic (it runs over P.size() coordinate rows), so the
 * same App drives the 2D and the 3D editor: `dim` only decides how many rows
 * the control points contribute and which canvas draws the result -- the
 * symbolic basis does not know the dimension, so switching costs a reeval(),
 * never a rebuild(). All cached geometry is stored in 3D (z = 0 in 2D mode);
 * only the zonotope elements keep two forms, because the two canvases consume
 * them differently (2D: the boundary polygon, cached; 3D: center + generators,
 * projected to an exact silhouette by the camera each frame).
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
#include "zonotope3d.hpp"

namespace sb2gui {

// One control point, stored in 3D (z components sit at 0 and are simply not
// sent to SB2 while dim == 2). Interpretation depends on the control set:
//   R : center only.  IR : box center +/- (hx,hy,hz).  Z : center + generators.
// A zonotope's shape lives entirely in its generators: each one is an arbitrary
// vector, so both the extent *and* the orientation are edited by moving the
// generator tips (c +/- g_i), never through an axis-aligned box.
struct ControlPoint {
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double hx = 0.1, hy = 0.1, hz = 0.1;   // half-widths (IR)
    std::vector<Vec3> gens;                // generators (Z)
};

struct Box {
    double x0, y0, z0, x1, y1, z1;
};

// A result zonotope in 3D form: center + generators, as extracted from the
// affine vector. The 3D canvas projects it through the (orthographic) camera
// and runs the 2D boundary walk on the projected generators -- the exact
// silhouette, recomputed per frame because it depends on the view.
struct Zono3 {
    Vec3 c;
    std::vector<Vec3> gens;
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
    // Set the B-spline parameter u is taken in: it decides how the basis itself is
    // evaluated (real samples, interval decomposition, or affine forms), so it is
    // baked into the SB2 object and changing it costs a rebuild().
    sb2l::ParameterSet ps = sb2l::ParameterSet::IR;
    std::vector<int> wnum, wden;       // rational weights, per control point

    // Set the control points are taken in: it selects eval_point / eval_box /
    // eval_zonotope on the *same* basis, so it is free of any rebuild (reeval()).
    // It is independent of ps -- every one of the 9 pairs is a valid B-spline, e.g.
    // real control points with an interval parameter give the tube of the curve.
    sb2l::ParameterSet cs = sb2l::ParameterSet::IR;

    // Spatial dimension of the curve: 2 or 3. SB2 takes it as the number of
    // coordinate rows, so switching costs a reeval(), never a rebuild().
    int dim = 2;

    std::vector<ControlPoint> cps;

    // Cached draw geometry (world coordinates), refreshed by reeval().
    // Segment s owns a fixed slice of each container: the d elements at [s*d, s*d + d),
    // which is what lets a segment be redrawn in place (an R *parameter* evaluates one
    // extra element at the very end, closing the curve, hence the trailing slot).
    std::vector<std::vector<Vec3>> polylines;  // cs = R: the curve, as a single polyline
    std::vector<Box> boxes;                    // cs = IR: result boxes
    std::vector<std::vector<Vec2>> zonos;      // cs = Z, dim 2: result zonotope polygons
    std::vector<Zono3> zonos3;                 // cs = Z, dim 3: center + generators

    std::string status;                        // last rebuild/eval message
    bool dirty_structural = true;
    bool want_fit = false;                     // consumed by the editor: refit view

    App();

    // Rebuild the SB2 object (resizes/preserves control points and weights).
    void rebuild();
    // Re-run evaluation with the current control points (no rebuild).
    void reeval();
    // Same, but only for the segments the control point i belongs to: this is the
    // drag path, and it costs p+1 segments instead of nCP - p.
    void update_control_point(int i);
    // Switch between the 2D and the 3D curve: the basis is untouched, so this
    // only re-runs the evaluation over a different number of coordinate rows.
    void set_dim(int nd);

    int nSegments() const { return spline ? spline->get_nS() : 0; }

private:
    std::unique_ptr<sb2l::SB2> spline;

    // Raw evaluation kept per segment, so that a moved control point can be pushed
    // through the impacted segments only. Exactly one of them is filled, per cs.
    std::vector<std::vector<std::vector<double>>> epoints;
    std::vector<std::vector<ibex::IntervalVector>> eboxes;
    std::vector<std::vector<ibex::Affine2Vector>> ezonos;

    void resize_control_points(); // keep cps/weights consistent with nCP/p
    void seed_default_scene();    // dim-aware: sine wave in 2D, helix in 3D

    // Control points of the current scene, in the layout expected by SB2::eval_*:
    // one row per coordinate, dim rows in all.
    std::vector<std::vector<double>> control_points_R() const;
    std::vector<ibex::IntervalVector> control_points_IR() const;
    std::vector<ibex::Affine2Vector> control_points_Z() const;
    void alloc_geometry();                  // size the draw containers for the current spline
    void refresh_geometry(int s0, int s1);  // rebuild the geometry of segments [s0, s1] from the raw eval
};

} // namespace sb2gui

#endif // SB2L_GUI_APP_HPP_
