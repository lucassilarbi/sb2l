#include "app.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace sb2gui {

// Default generators for a fresh zonotope control point (a small parallelogram).
static std::vector<Vec2> default_gens()
{
    return {{0.12, 0.04}, {-0.04, 0.12}};
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
        cp.gens = default_gens();
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
            cp.gens = default_gens();
            cps.push_back(cp);
        }
    } else if ((int)cps.size() > nCP) {
        cps.resize(nCP);
    }
    wnum.resize(nCP, 1);
    wden.resize(nCP, 1);
}

void App::rebuild()
{
    if (p < 1) p = 1;
    if (nCP < p + 1) nCP = p + 1;
    if (d < 1) d = 1;
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
        // Track generators precisely in affine mode (fast boundary removes the
        // old 2^m concern, so we can afford a generous noise budget).
        ibex::AF_fAFFullI::setAffineNoiseNumber(2 * nCP + 16);

        spline.reset(new sb2l::SB2(p, nCP, ct, f, ps, d, t, W));
        dirty_structural = false;
        status = "OK";
        reeval();
    } catch (const std::exception& e) {
        spline.reset();
        polylines.clear();
        boxes.clear();
        zonos.clear();
        status = std::string("rebuild failed: ") + e.what();
    }
}

// Boundary polygon of a 2D affine vector (center + generators + error box).
static std::vector<Vec2> zono_from_affine(ibex::Affine2Vector v)
{
    v.compact();
    Vec2 c{v[0].val(0), v[1].val(0)};
    int n = v[0].size();
    if (v[1].size() > n) n = v[1].size();
    std::vector<Vec2> gens;
    gens.reserve(n + 2);
    for (int k = 1; k <= n; ++k) {
        double gx = (k <= v[0].size()) ? v[0].val(k) : 0.0;
        double gy = (k <= v[1].size()) ? v[1].val(k) : 0.0;
        gens.push_back({gx, gy});
    }
    // Accumulated error terms as axis-aligned generators (keeps over-approx).
    if (v[0].err() > 0.0) gens.push_back({v[0].err(), 0.0});
    if (v[1].err() > 0.0) gens.push_back({0.0, v[1].err()});
    return zonotope_boundary(c, gens);
}

void App::reeval()
{
    if (!spline) return;
    polylines.clear();
    boxes.clear();
    zonos.clear();
    try {
        if (ps == sb2l::ParameterSet::R) {
            std::vector<std::vector<double>> P(2, std::vector<double>(nCP));
            for (int i = 0; i < nCP; ++i) { P[0][i] = cps[i].cx; P[1][i] = cps[i].cy; }
            auto pts = spline->eval_point(P);
            // One continuous polyline: joining segments end-to-end removes the
            // gaps between per-segment sample runs.
            std::vector<Vec2> line;
            for (auto& seg : pts)
                for (auto& q : seg) line.push_back({q[0], q[1]});
            if (line.size() >= 2) polylines.push_back(std::move(line));
        } else if (ps == sb2l::ParameterSet::IR) {
            std::vector<ibex::IntervalVector> P(2, ibex::IntervalVector(nCP));
            for (int i = 0; i < nCP; ++i) {
                P[0][i] = ibex::Interval(cps[i].cx - cps[i].hx, cps[i].cx + cps[i].hx);
                P[1][i] = ibex::Interval(cps[i].cy - cps[i].hy, cps[i].cy + cps[i].hy);
            }
            auto bx = spline->eval_box(P);
            for (auto& seg : bx)
                for (auto& b : seg)
                    boxes.push_back({b[0].lb(), b[1].lb(), b[0].ub(), b[1].ub()});
        } else { // Z
            std::vector<ibex::Affine2Vector> aP(2, ibex::Affine2Vector(nCP));
            for (int i = 0; i < nCP; ++i) {
                ibex::Affine2 ax(cps[i].cx), ay(cps[i].cy);
                for (const Vec2& g : cps[i].gens) {
                    ibex::Affine2 e(ibex::Interval(-1, 1)); // fresh independent noise
                    ax += g.x * e;
                    ay += g.y * e;
                }
                aP[0][i] = ax;
                aP[1][i] = ay;
            }
            auto zs = spline->eval_zonotope(aP);
            for (auto& seg : zs)
                for (auto& z : seg)
                    zonos.push_back(zono_from_affine(z));
        }
        status = "OK";
    } catch (const std::exception& e) {
        status = std::string("eval failed: ") + e.what();
    }
}

} // namespace sb2gui
