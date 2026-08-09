/**
 * Fast 3D zonotope boundary (center + generators -> facet mesh).
 *
 * Header-only, pure double: no ibex / OpenGL / sb2l dependency, so it can be
 * unit-tested and reused freely -- the 3D sibling of zonotope2d.hpp.
 *
 * No general convex-hull run is ever needed: every facet of a 3D zonotope is
 * a parallelogram spanned by a *pair* of (non-parallel) generators, translated
 * to the boundary by putting every other generator on its extreme side. So
 * the m(m-1) facets are enumerated directly, in O(m^2) facets and O(m^3)
 * arithmetic -- exact, and far below a QHull call on the 2^m vertex cloud.
 *
 * That cost is kept interactive by reduce_generators(): parallel generators
 * are merged (exact -- a zonotope is a Minkowski sum, so collinear segments
 * add lengthwise) and the smallest ones are outer-approximated by one
 * axis-aligned box (3 generators), the standard conservative order reduction.
 * The many *result* zonotopes of an evaluation never even need the mesh: an
 * orthographic projection of a zonotope is again a zonotope, so their screen
 * silhouette is the O(m log m) 2D walk of zonotope2d.hpp on the projected
 * generators.
 */
#ifndef SB2L_GUI_ZONOTOPE3D_HPP_
#define SB2L_GUI_ZONOTOPE3D_HPP_

#include <algorithm>
#include <cmath>
#include <vector>

namespace sb2gui {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, const Vec3& a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

/** One boundary facet: a parallelogram, vertices in CCW order seen from
 *  outside, with its outward unit normal (used for culling and shading). */
struct ZFace {
    Vec3 n;
    Vec3 v[4];
};

/**
 * Conservative order reduction of a generator list.
 *
 * 1. Zero generators are dropped and (near-)parallel ones merged: g_i and g_j
 *    with |g_i x g_j| <= eps |g_i||g_j| collapse to g_i + s g_j (s the sign of
 *    their dot product), which is *exact* for the zonotope they span.
 * 2. If more than max_g remain, the smallest are replaced by the axis-aligned
 *    box of their sum of absolute components (up to 3 generators): the result
 *    encloses the original zonotope, never shrinks it.
 */
inline std::vector<Vec3> reduce_generators(const std::vector<Vec3>& gens, int max_g)
{
    std::vector<Vec3> g;
    g.reserve(gens.size());
    for (const Vec3& e : gens)
        if (e.x != 0.0 || e.y != 0.0 || e.z != 0.0) g.push_back(e);

    // Largest first, so the merge keeps the dominant direction as pivot.
    const auto longest_first = [](const Vec3& a, const Vec3& b) { return dot(a, a) > dot(b, b); };
    std::sort(g.begin(), g.end(), longest_first);

    const double eps = 1e-12;
    std::vector<Vec3> kept;
    kept.reserve(g.size());
    for (const Vec3& e : g) {
        bool merged = false;
        for (Vec3& k : kept) {
            const Vec3 c = cross(k, e);
            if (dot(c, c) <= eps * dot(k, k) * dot(e, e)) {
                k = k + (dot(k, e) >= 0.0 ? 1.0 : -1.0) * e;
                merged = true;
                break;
            }
        }
        if (!merged) kept.push_back(e);
    }

    // Merging lengthens the generator it merges into, so the order above no
    // longer holds. Restore it, otherwise the box below swallows generators
    // which are not the smallest ones and the enclosure is looser than needed.
    std::sort(kept.begin(), kept.end(), longest_first);

    if (max_g >= 3 && (int)kept.size() > max_g) {
        Vec3 box{0.0, 0.0, 0.0};
        for (size_t i = max_g - 3; i < kept.size(); ++i) {
            box.x += std::fabs(kept[i].x);
            box.y += std::fabs(kept[i].y);
            box.z += std::fabs(kept[i].z);
        }
        kept.resize(max_g - 3);
        if (box.x > 0.0) kept.push_back({box.x, 0.0, 0.0});
        if (box.y > 0.0) kept.push_back({0.0, box.y, 0.0});
        if (box.z > 0.0) kept.push_back({0.0, 0.0, box.z});
    }
    return kept;
}

/**
 * Boundary facets of the zonotope { c + sum_i s_i * g_i : s_i in [-1,1] }.
 *
 * The generators are assumed pairwise non-parallel (run reduce_generators
 * first); parallel leftovers are skipped, never mis-drawn. Degenerate inputs
 * fall out naturally: 0 or 1 generator (or all parallel) give no facet -- the
 * zonotope is a point or a segment, which the caller draws directly; exactly
 * 2 give the two sides of one flat parallelogram.
 *
 * For every unordered pair (i, j), the facet normal is n = g_i x g_j and the
 * facet center is c + sum_{k != i,j} sign(n . g_k) g_k, the choice that puts
 * all remaining generators on the far side -- with its antipodal twin at the
 * mirrored center. When some other g_k is coplanar with the pair (n.g_k = 0)
 * the true facet is a hexagon or larger that the coplanar pairs tile between
 * them; the tie is broken deterministically (+) so the tiles stay consistent.
 */
inline std::vector<ZFace> zonotope_mesh(const Vec3& c, const std::vector<Vec3>& gens)
{
    const int m = (int)gens.size();
    std::vector<ZFace> faces;
    if (m < 2) return faces;
    faces.reserve((size_t)m * (m - 1));

    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            const Vec3& gi = gens[i];
            const Vec3& gj = gens[j];
            Vec3 n = cross(gi, gj);
            const double ln = norm(n);
            if (ln <= 1e-12 * norm(gi) * norm(gj)) continue; // parallel pair
            n = (1.0 / ln) * n;

            Vec3 off{0.0, 0.0, 0.0};
            for (int k = 0; k < m; ++k) {
                if (k == i || k == j) continue;
                off = off + (dot(n, gens[k]) >= 0.0 ? 1.0 : -1.0) * gens[k];
            }

            // CCW seen from the +n side, for a right-handed (gi, gj, n).
            const Vec3 fc = c + off;
            ZFace fp;
            fp.n = n;
            fp.v[0] = fc - gi - gj;
            fp.v[1] = fc + gi - gj;
            fp.v[2] = fc + gi + gj;
            fp.v[3] = fc - gi + gj;
            faces.push_back(fp);

            // The antipodal facet: mirrored center, opposite normal, and the
            // vertex order reversed so it stays CCW from *its* outside.
            const Vec3 fm = c - off;
            ZFace fn;
            fn.n = -1.0 * n;
            fn.v[0] = fm - gi - gj;
            fn.v[1] = fm - gi + gj;
            fn.v[2] = fm + gi + gj;
            fn.v[3] = fm + gi - gj;
            faces.push_back(fn);
        }
    }
    return faces;
}

} // namespace sb2gui

#endif // SB2L_GUI_ZONOTOPE3D_HPP_
