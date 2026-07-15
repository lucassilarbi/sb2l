#include "app.hpp"

#include <cmath>
#include <list>
#include <stdexcept>
#include <utility>

namespace sb2gui {

// Default generators for a fresh zonotope control point: n directions evenly
// spread over a half-turn (n = 2 gives the usual small square).
static std::vector<Vec2> default_gens(int n)
{
    std::vector<Vec2> g;
    for (int k = 0; k < n; ++k) {
        double th = 3.14159265358979 * k / (n > 0 ? n : 1);
        g.push_back({0.12 * std::cos(th), 0.12 * std::sin(th)});
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
    cps.clear();
    for (int i = 0; i < nCP; ++i) {
        double u = (nCP > 1) ? (double)i / (nCP - 1) : 0.0;
        ControlPoint cp;
        cp.cx = 0.3 + 2.6 * u;
        cp.cy = 1.5 + 1.0 * std::sin(u * 3.14159265 * 1.5);
        cp.hx = cp.hy = 0.1;
        cp.gens = default_gens(nGen);
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
            cp.hx = cp.hy = 0.1;
            cp.gens = default_gens(nGen);
            cps.push_back(cp);
        }
    } else if ((int)cps.size() > nCP) {
        cps.resize(nCP);
    }
    // Keep every point at nGen generators, preserving the ones already edited.
    for (ControlPoint& cp : cps) {
        if ((int)cp.gens.size() == nGen) continue;
        std::vector<Vec2> fresh = default_gens(nGen);
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
        ibex::AF_fAFFullI::setAffineNoiseNumber(nCP * nGen + 16);

        spline.reset(new sb2l::SB2(p, nCP, ct, f, ps, d, t, W));
        dirty_structural = false;
        status = "OK";
        reeval();
    } catch (const std::exception& e) {
        spline.reset();
        epoints.clear();
        eboxes.clear();
        ezonos.clear();
        polylines.clear();
        boxes.clear();
        zonos.clear();
        status = std::string("rebuild failed: ") + e.what();
    }
}

// Boundary polygon of a 2D affine vector (center + generators + error box).
// Takes the form already compacted: it is read by reference, so that drawing a
// zonotope does not deep-copy its noise-term lists every frame.
static std::vector<Vec2> zono_from_affine(const ibex::Affine2Vector& v)
{
    Vec2 c{v[0].val(0), v[1].val(0)};

    // Both coordinates draw on the same noise symbols but each stores only its
    // own nonzero terms, so walk the two index-sorted lists together and emit
    // one generator per symbol either of them uses.
    //
    // Do not iterate k = 1 .. size() calling val(k) instead: noise indices come
    // from a global counter that every affine multiplication bumps and nothing
    // ever resets, so size() is a symbol index in the millions after a few
    // seconds of dragging, not a term count, and that loop grows without bound.
    const std::list<std::pair<int, double> >& rx = v[0].rays();
    const std::list<std::pair<int, double> >& ry = v[1].rays();
    std::vector<Vec2> gens;
    gens.reserve(rx.size() + ry.size() + 2);
    std::list<std::pair<int, double> >::const_iterator ix = rx.begin(), iy = ry.begin();
    while (ix != rx.end() || iy != ry.end()) {
        if (iy == ry.end() || (ix != rx.end() && ix->first < iy->first))
            gens.push_back({(ix++)->second, 0.0});
        else if (ix == rx.end() || iy->first < ix->first)
            gens.push_back({0.0, (iy++)->second});
        else
            gens.push_back({(ix++)->second, (iy++)->second});
    }
    // Accumulated error terms as axis-aligned generators (keeps over-approx).
    if (v[0].err() > 0.0) gens.push_back({v[0].err(), 0.0});
    if (v[1].err() > 0.0) gens.push_back({0.0, v[1].err()});
    return zonotope_boundary(c, gens);
}

std::vector<std::vector<double>> App::control_points_R() const
{
    std::vector<std::vector<double>> P(2, std::vector<double>(nCP));
    for (int i = 0; i < nCP; ++i) { P[0][i] = cps[i].cx; P[1][i] = cps[i].cy; }
    return P;
}

std::vector<ibex::IntervalVector> App::control_points_IR() const
{
    std::vector<ibex::IntervalVector> P(2, ibex::IntervalVector(nCP));
    for (int i = 0; i < nCP; ++i) {
        P[0][i] = ibex::Interval(cps[i].cx - cps[i].hx, cps[i].cx + cps[i].hx);
        P[1][i] = ibex::Interval(cps[i].cy - cps[i].hy, cps[i].cy + cps[i].hy);
    }
    return P;
}

std::vector<ibex::Affine2Vector> App::control_points_Z() const
{
    std::vector<ibex::Affine2Vector> P(2, ibex::Affine2Vector(nCP));
    for (int i = 0; i < nCP; ++i) {
        ibex::Affine2 ax(cps[i].cx), ay(cps[i].cy);
        for (const Vec2& g : cps[i].gens) {
            ibex::Affine2 e(ibex::Interval(-1, 1)); // fresh independent noise
            ax += g.x * e;
            ay += g.y * e;
        }
        P[0][i] = ax;
        P[1][i] = ay;
    }
    return P;
}

void App::alloc_geometry()
{
    const int nS = spline->get_nS(), d = spline->get_d();
    // Only a real parameter evaluates the closing element of the last segment
    // (SB2::du_count), whatever set the control points live in.
    const int n = nS * d + (ps == sb2l::ParameterSet::R ? 1 : 0);
    polylines.clear();
    boxes.clear();
    zonos.clear();
    if (cs == sb2l::ParameterSet::R) {
        // One continuous polyline: joining the segments end-to-end removes the
        // gaps between per-segment sample runs.
        polylines.assign(1, std::vector<Vec2>(n));
    } else if (cs == sb2l::ParameterSet::IR) {
        boxes.assign(n, Box{});
    } else {
        zonos.assign(n, std::vector<Vec2>());
    }
}

void App::refresh_geometry(int s0, int s1)
{
    const int nS = spline->get_nS(), d = spline->get_d();
    if (s0 < 0) s0 = 0;
    if (s1 > nS - 1) s1 = nS - 1;
    for (int s = s0; s <= s1; ++s) {
        if (cs == sb2l::ParameterSet::R) {
            const std::vector<std::vector<double>>& seg = epoints[s];
            for (size_t du = 0; du < seg.size(); ++du)
                polylines[0][s * d + du] = {seg[du][0], seg[du][1]};
        } else if (cs == sb2l::ParameterSet::IR) {
            const std::vector<ibex::IntervalVector>& seg = eboxes[s];
            for (size_t du = 0; du < seg.size(); ++du) {
                const ibex::IntervalVector& b = seg[du];
                boxes[s * d + du] = {b[0].lb(), b[1].lb(), b[0].ub(), b[1].ub()};
            }
        } else {
            // Compact in place: eval/update_zonotope always rewrites these forms
            // wholesale before they are read again, so dropping their negligible
            // noise terms here costs nothing and keeps the polygons small.
            std::vector<ibex::Affine2Vector>& seg = ezonos[s];
            for (size_t du = 0; du < seg.size(); ++du) {
                seg[du].compact();
                zonos[s * d + du] = zono_from_affine(seg[du]);
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
    try {
        if (cs == sb2l::ParameterSet::R)
            epoints = spline->eval_point(control_points_R());
        else if (cs == sb2l::ParameterSet::IR)
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
        polylines.clear();
        boxes.clear();
        zonos.clear();
        status = std::string("eval failed: ") + e.what();
    }
}

void App::update_control_point(int i)
{
    if (!spline) return;
    const int nS = spline->get_nS();
    const size_t cached = (cs == sb2l::ParameterSet::R)  ? epoints.size()
                        : (cs == sb2l::ParameterSet::IR) ? eboxes.size()
                                                         : ezonos.size();
    if ((int)cached != nS) { // nothing to patch: mode just changed, or the last eval failed
        reeval();
        return;
    }
    try {
        const std::pair<int, int> seg = spline->impacted_segments(i);
        if (cs == sb2l::ParameterSet::R)
            spline->update_point(control_points_R(), i, epoints);
        else if (cs == sb2l::ParameterSet::IR)
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
