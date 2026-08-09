#include "app.hpp"

#include <cmath>
#include <exception>
#include <utility>

namespace sb2gui {

static const double kPi = 3.14159265358979323846;

// Default generators for a fresh zonotope control point. 2D: n directions
// evenly spread over a half-turn (n = 2 gives the usual small square).
// 3D: n directions on a Fibonacci hemisphere (a zonotope is symmetric about
// its center, so half the sphere spans every orientation).
static std::vector<Vec3> default_gens(int n, int dim)
{
    std::vector<Vec3> g;
    if (n < 1) n = 1;
    if (dim == 3) {
        const double ga = kPi * (3.0 - std::sqrt(5.0)); // golden angle
        for (int k = 0; k < n; ++k) {
            const double zk = (k + 0.5) / n;
            const double r = std::sqrt(1.0 - zk * zk);
            g.push_back({0.12 * r * std::cos(ga * k), 0.12 * r * std::sin(ga * k), 0.12 * zk});
        }
    } else {
        for (int k = 0; k < n; ++k) {
            const double th = kPi * k / n;
            g.push_back({0.12 * std::cos(th), 0.12 * std::sin(th), 0.0});
        }
    }
    return g;
}

App::App()
{
    seed_default_scene();
    rebuild();
}

void App::seed_default_scene()
{
    scene_edited = false;
    cps.clear();
    for (int i = 0; i < nCP; ++i) {
        double u = (nCP > 1) ? (double)i / (nCP - 1) : 0.0;
        ControlPoint cp;
        if (dim == 3) {
            // A genuinely 3D scene: a helix, so the curve does not sit in any
            // coordinate plane and the orbit camera has something to show.
            const double th = 2.0 * kPi * u;
            cp.cx = 1.5 * std::cos(th);
            cp.cy = 1.5 * std::sin(th);
            cp.cz = 2.4 * u - 1.2;
        } else {
            cp.cx = 0.3 + 2.6 * u;
            cp.cy = 1.5 + 1.0 * std::sin(u * kPi * 1.5);
            cp.cz = 0.0;
        }
        cp.hx = cp.hy = cp.hz = 0.1;
        cp.gens = default_gens(nGen, dim);
        cps.push_back(cp);
    }
}

void App::resize_control_points()
{
    if ((int)cps.size() < nCP) {
        for (int i = (int)cps.size(); i < nCP; ++i) {
            ControlPoint cp;
            const ControlPoint& prev = cps.empty() ? ControlPoint{} : cps.back();
            cp.cx = prev.cx + 0.3;
            cp.cy = prev.cy;
            cp.cz = prev.cz;
            cp.hx = cp.hy = cp.hz = 0.1;
            cp.gens = default_gens(nGen, dim);
            cps.push_back(cp);
        }
    } else if ((int)cps.size() > nCP) {
        cps.resize(nCP);
    }
    // Keep every point at nGen generators, preserving the ones already edited.
    for (ControlPoint& cp : cps) {
        if ((int)cp.gens.size() == nGen) continue;
        std::vector<Vec3> fresh = default_gens(nGen, dim);
        for (int k = 0; k < nGen && k < (int)cp.gens.size(); ++k) fresh[k] = cp.gens[k];
        cp.gens = std::move(fresh);
    }
    wnum.resize(nCP, 1);
    wden.resize(nCP, 1);
}

void App::rebuild()
{
    if (p < 1) p = 1;
    if (nCP < p + 1) nCP = p + 1;
    if (d < 1) d = 1;
    if (nGen < 1) nGen = 1;
    if (dim != 3) dim = 2;
    resize_control_points();

    try {
        std::vector<SymEngine::Expression> W;
        bool rational = (ct == sb2l::CurveType::UNIFORM_RATIONAL ||
                         ct == sb2l::CurveType::CLAMPED_RATIONAL);
        if (rational) {
            for (int i = 0; i < nCP; ++i) {
                int num = wnum[i], den = wden[i] == 0 ? 1 : wden[i];
                W.push_back(SymEngine::Expression(num) / SymEngine::Expression(den));
            }
        }
        // One independent noise symbol per generator of every control point,
        // plus room for the ones the evaluation introduces itself.
        sb2l::SB2::set_affine_noise_number(nCP * nGen + 16);

        spline.reset(new sb2l::SB2(p, nCP, ct, f, ps, d, t, W));
        status = "OK";
        reeval();
    } catch (const std::exception& e) {
        spline.reset();
        epoints.clear();
        eboxes.clear();
        ezonos.clear();
        points.clear();
        boxes.clear();
        zonos.clear();
        zonos3.clear();
        status = std::string("rebuild failed: ") + e.what();
    }
}

void App::set_dim(int nd)
{
    if (nd != 2 && nd != 3) return;
    if (nd == dim) return;
    dim = nd;
    // A zonotope with m generators spans at most m dimensions, so the usual 2D
    // default of two generators would make every control set a flat
    // parallelogram in space. Three is the least that gives a solid.
    if (dim == 3 && nGen < 3) nGen = 3;
    // An edited scene is kept: the control points already hold their three
    // coordinates, and the symbolic basis does not know the dimension, so
    // there is nothing to throw away. A scene still untouched is the default
    // one, and the default is written for its dimension: the plane curve is
    // replaced by the helix, which is the one that leaves every coordinate
    // plane and gives the orbit camera something to show.
    if (!scene_edited)
        seed_default_scene();
    resize_control_points();
    reeval();
    want_fit = true;
}

void App::reset_scene()
{
    seed_default_scene();
    resize_control_points();
    reeval();
    want_fit = true;
}

void App::set_generator_count(int n)
{
    if (n < 1) n = 1;
    if (n == nGen) return;
    nGen = n;
    // The number of generators changes the control sets, never the basis, so
    // the symbolic rebuild is not needed: only the noise budget and one fresh
    // evaluation.
    resize_control_points();
    sb2l::SB2::set_affine_noise_number(nCP * nGen + 16);
    reeval();
}

// Boundary polygon of a 2D affine vector (center + generators + error box).
// Takes the form already compacted: it is read by reference, so that drawing a
// zonotope does not deep-copy its noise-term lists every frame.
// The element is read through sb2l::zonotope_of, the one place which knows how
// a set is taken out of an affine vector: reading it by hand is where the
// noise terms get lost (see the warning on that function). Only the layout
// changes here, into what the two canvases draw.
static std::vector<Vec2> zono_from_affine(const ibex::Affine2Vector& v)
{
    const sb2l::Zonotope z = sb2l::zonotope_of(v);
    const Vec2 c{z.center[0], z.center[1]};
    std::vector<Vec2> gens;
    gens.reserve(z.m);
    for (int k = 0; k < z.m; ++k) gens.push_back({z.gen(k, 0), z.gen(k, 1)});
    return zonotope_boundary(c, gens);
}

// 3D center + generators. The boundary itself is not built here -- the camera
// projects these generators and the 2D walk runs on the projection, so the
// exact silhouette costs O(m log m) per frame and the O(m^2) facet mesh is
// never needed for result elements.
static Zono3 zono3_from_affine(const ibex::Affine2Vector& v)
{
    const sb2l::Zonotope zo = sb2l::zonotope_of(v);
    Zono3 z;
    z.c = {zo.center[0], zo.center[1], zo.center[2]};
    z.gens.reserve(zo.m);
    for (int k = 0; k < zo.m; ++k) z.gens.push_back({zo.gen(k, 0), zo.gen(k, 1), zo.gen(k, 2)});
    return z;
}

std::vector<std::vector<double>> App::control_points_R() const
{
    std::vector<std::vector<double>> P(dim, std::vector<double>(nCP));
    for (int i = 0; i < nCP; ++i) {
        P[0][i] = cps[i].cx;
        P[1][i] = cps[i].cy;
        if (dim == 3) P[2][i] = cps[i].cz;
    }
    return P;
}

std::vector<ibex::IntervalVector> App::control_points_IR() const
{
    // Real control points are degenerate intervals: zero half-widths.
    const bool degen = (cs == sb2l::ParameterSet::R);
    std::vector<ibex::IntervalVector> P(dim, ibex::IntervalVector(nCP));
    for (int i = 0; i < nCP; ++i) {
        const double hx = degen ? 0.0 : cps[i].hx;
        const double hy = degen ? 0.0 : cps[i].hy;
        const double hz = degen ? 0.0 : cps[i].hz;
        P[0][i] = ibex::Interval(cps[i].cx - hx, cps[i].cx + hx);
        P[1][i] = ibex::Interval(cps[i].cy - hy, cps[i].cy + hy);
        if (dim == 3) P[2][i] = ibex::Interval(cps[i].cz - hz, cps[i].cz + hz);
    }
    return P;
}

std::vector<ibex::Affine2Vector> App::control_points_Z() const
{
    // Real control points are constant affine forms: no generator at all.
    const bool degen = (cs == sb2l::ParameterSet::R);
    std::vector<ibex::Affine2Vector> P(dim, ibex::Affine2Vector(nCP));
    for (int i = 0; i < nCP; ++i) {
        ibex::Affine2 ax(cps[i].cx), ay(cps[i].cy), az(cps[i].cz);
        if (!degen)
            for (const Vec3& g : cps[i].gens) {
                ibex::Affine2 e(ibex::Interval(-1, 1)); // fresh independent noise
                ax += g.x * e;
                ay += g.y * e;
                if (dim == 3) az += g.z * e;
            }
        P[0][i] = ax;
        P[1][i] = ay;
        if (dim == 3) P[2][i] = az;
    }
    return P;
}

void App::alloc_geometry()
{
    const int nS = spline->get_nS(), dd = spline->get_d();
    // Only a real parameter evaluates the closing element of the last segment
    // (SB2::du_count), whatever set the control points live in.
    const int n = nS * dd + (ps == sb2l::ParameterSet::R ? 1 : 0);
    const sb2l::ParameterSet ec = effective_cs();
    points.clear();
    boxes.clear();
    zonos.clear();
    zonos3.clear();
    if (ec == sb2l::ParameterSet::R) {
        // The evaluated points of the real-real pair, one slot per sample.
        points.assign(n, Vec3{});
    } else if (ec == sb2l::ParameterSet::IR) {
        boxes.assign(n, Box{});
    } else if (dim == 3) {
        zonos3.assign(n, Zono3{});
    } else {
        zonos.assign(n, std::vector<Vec2>());
    }
}

void App::refresh_geometry(int s0, int s1)
{
    const int nS = spline->get_nS(), dd = spline->get_d();
    const sb2l::ParameterSet ec = effective_cs();
    if (s0 < 0) s0 = 0;
    if (s1 > nS - 1) s1 = nS - 1;
    for (int s = s0; s <= s1; ++s) {
        if (ec == sb2l::ParameterSet::R) {
            const std::vector<std::vector<double>>& seg = epoints[s];
            for (size_t du = 0; du < seg.size(); ++du)
                points[s * dd + du] = {seg[du][0], seg[du][1], dim == 3 ? seg[du][2] : 0.0};
        } else if (ec == sb2l::ParameterSet::IR) {
            const std::vector<ibex::IntervalVector>& seg = eboxes[s];
            for (size_t du = 0; du < seg.size(); ++du) {
                const ibex::IntervalVector& b = seg[du];
                boxes[s * dd + du] = {b[0].lb(), b[1].lb(), dim == 3 ? b[2].lb() : 0.0,
                                      b[0].ub(), b[1].ub(), dim == 3 ? b[2].ub() : 0.0};
            }
        } else {
            // Compact in place: eval/update_zonotope always rewrites these forms
            // wholesale before they are read again, so dropping their negligible
            // noise terms here costs nothing and keeps the polygons small.
            std::vector<ibex::Affine2Vector>& seg = ezonos[s];
            for (size_t du = 0; du < seg.size(); ++du) {
                seg[du].compact();
                if (dim == 3)
                    zonos3[s * dd + du] = zono3_from_affine(seg[du]);
                else
                    zonos[s * dd + du] = zono_from_affine(seg[du]);
            }
        }
    }
}

void App::reeval()
{
    if (!spline) return;
    epoints.clear();
    eboxes.clear();
    ezonos.clear();
    const sb2l::ParameterSet ec = effective_cs();
    try {
        if (ec == sb2l::ParameterSet::R)
            epoints = spline->eval_point(control_points_R());
        else if (ec == sb2l::ParameterSet::IR)
            eboxes = spline->eval_box(control_points_IR());
        else
            ezonos = spline->eval_zonotope(control_points_Z());
        alloc_geometry();
        refresh_geometry(0, spline->get_nS() - 1);
        status = "OK";
    } catch (const std::exception& e) {
        epoints.clear();
        eboxes.clear();
        ezonos.clear();
        points.clear();
        boxes.clear();
        zonos.clear();
        zonos3.clear();
        status = std::string("eval failed: ") + e.what();
    }
}

void App::update_control_point(int i)
{
    if (!spline) return;
    scene_edited = true; // the scene stops being the default one it was seeded with
    const int nS = spline->get_nS();
    const sb2l::ParameterSet ec = effective_cs();
    const size_t cached = (ec == sb2l::ParameterSet::R)  ? epoints.size()
                        : (ec == sb2l::ParameterSet::IR) ? eboxes.size()
                                                         : ezonos.size();
    if ((int)cached != nS) { // nothing to patch: mode just changed, or the last eval failed
        reeval();
        return;
    }
    // The raw eval must carry the current dimension too: the containers keep
    // their per-segment layout across a dim switch, but not their row count.
    const int cached_dim =
        (ec == sb2l::ParameterSet::R)  ? (epoints[0].empty() ? 0 : (int)epoints[0][0].size())
      : (ec == sb2l::ParameterSet::IR) ? (eboxes[0].empty() ? 0 : eboxes[0][0].size())
                                       : (ezonos[0].empty() ? 0 : ezonos[0][0].size());
    if (cached_dim != dim) {
        reeval();
        return;
    }
    try {
        const std::pair<int, int> seg = spline->impacted_segments(i);
        if (ec == sb2l::ParameterSet::R)
            spline->update_point(control_points_R(), i, epoints);
        else if (ec == sb2l::ParameterSet::IR)
            spline->update_box(control_points_IR(), i, eboxes);
        else
            spline->update_zonotope(control_points_Z(), i, ezonos);
        refresh_geometry(seg.first, seg.second);
        status = "OK";
    } catch (const std::exception& e) {
        status = std::string("eval failed: ") + e.what();
    }
}

} // namespace sb2gui
