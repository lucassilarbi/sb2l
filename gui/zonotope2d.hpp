/**
 * Fast 2D zonotope boundary (center + generators -> convex polygon).
 *
 * Header-only, pure double: no ibex / OpenGL / sb2l dependency, so it can be
 * unit-tested and reused freely. Replaces the exponential 2^m enumeration +
 * QHull used in the examples by the standard O(m log m) Minkowski-sum walk: a
 * 2D zonotope with m nonzero generators has exactly 2m boundary vertices,
 * produced here already in CCW order.
 */
#ifndef SB2L_GUI_ZONOTOPE2D_HPP_
#define SB2L_GUI_ZONOTOPE2D_HPP_

#include <algorithm>
#include <cmath>
#include <vector>

namespace sb2gui {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

/**
 * Boundary vertices of the zonotope { c + sum_i s_i * g_i : s_i in [-1,1] }.
 * Returns the vertices in counter-clockwise order (2m for m nonzero
 * generators; a single point for m == 0; a segment's 2 endpoints for m == 1).
 */
inline std::vector<Vec2> zonotope_boundary(const Vec2& c, const std::vector<Vec2>& gens)
{
    // Fold every generator into the upper half-plane [0, pi); the zonotope is
    // symmetric about c so orientation of a generator does not change the set.
    std::vector<Vec2> g;
    g.reserve(gens.size());
    for (const Vec2& e : gens) {
        if (e.x == 0.0 && e.y == 0.0) continue;
        if (e.y < 0.0 || (e.y == 0.0 && e.x < 0.0))
            g.push_back({-e.x, -e.y});
        else
            g.push_back({e.x, e.y});
    }
    if (g.empty()) return {c};

    std::sort(g.begin(), g.end(), [](const Vec2& a, const Vec2& b) {
        return std::atan2(a.y, a.x) < std::atan2(b.y, b.x);
    });

    // Start at the bottom-most vertex: c minus the sum of all generators.
    Vec2 p = c;
    for (const Vec2& e : g) { p.x -= e.x; p.y -= e.y; }

    std::vector<Vec2> pts;
    pts.reserve(2 * g.size());
    for (const Vec2& e : g) {                 // upper chain
        pts.push_back(p);
        p.x += 2.0 * e.x;
        p.y += 2.0 * e.y;
    }
    for (const Vec2& e : g) {                 // lower chain (point-symmetric)
        pts.push_back(p);
        p.x -= 2.0 * e.x;
        p.y -= 2.0 * e.y;
    }
    return pts;
}

} // namespace sb2gui

#endif // SB2L_GUI_ZONOTOPE2D_HPP_
