/**
 * Checks of the fast 3D zonotope boundary (gui/zonotope3d.hpp) and of its
 * consistency with the 2D walk (gui/zonotope2d.hpp) used for silhouettes.
 *
 * Pure double, no ibex/sb2l dependency: the mesh is a geometric object and is
 * checked as one.
 *
 *   - the facet mesh is *exact*: along every facet normal, the plane of the
 *     facet supports the zonotope (offset = sum of |n . g_i|);
 *   - every vertex of the zonotope (all sign choices) is enclosed by the mesh;
 *   - face count: m(m-1) facets for m pairwise non-parallel generators, 6 for
 *     the box;
 *   - reduce_generators merges parallels exactly and never shrinks the set
 *     (support function only ever grows);
 *   - the orthographic-projection silhouette (2D walk on projected
 *     generators) encloses the projection of every mesh vertex -- the exact
 *     property the 3D canvas rendering rests on.
 */
#include "../gui/zonotope2d.hpp"
#include "../gui/zonotope3d.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using sb2gui::Vec2;
using sb2gui::Vec3;
using sb2gui::ZFace;

static int failures = 0;

static void check(bool ok, const char* what)
{
    std::printf("%s: %s\n", ok ? "ok" : "FAIL", what);
    if (!ok) ++failures;
}

// Deterministic pseudo-random doubles in [-1, 1] (no <random>: reproducible).
static double rnd()
{
    static unsigned long long s = 88172645463325252ull;
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return 2.0 * ((double)(s % 1000003ull) / 1000003.0) - 1.0;
}

// Support function of the zonotope along n: sum_i |n . g_i|.
static double support(const std::vector<Vec3>& g, const Vec3& n)
{
    double s = 0.0;
    for (const Vec3& e : g) s += std::fabs(dot(n, e));
    return s;
}

// Signed distance of p to the outside of the mesh (max over facet planes).
static double outside_by(const std::vector<ZFace>& mesh, const Vec3& c, const Vec3& p)
{
    double worst = -1e30;
    for (const ZFace& f : mesh) {
        const double d = dot(f.n, p - f.v[0]);
        if (d > worst) worst = d;
    }
    (void)c;
    return worst;
}

int main()
{
    const double eps = 1e-9;

    // ---- the box: 3 axis generators ----
    {
        const Vec3 c{1.0, -2.0, 0.5};
        const std::vector<Vec3> g = {{0.5, 0, 0}, {0, 0.7, 0}, {0, 0, 0.9}};
        const std::vector<ZFace> mesh = sb2gui::zonotope_mesh(c, g);
        check(mesh.size() == 6, "box mesh has 6 facets");
        bool planes = true;
        for (const ZFace& f : mesh) {
            // The facet plane must support the zonotope exactly.
            const double off = dot(f.n, f.v[0] - c);
            if (std::fabs(off - support(g, f.n)) > eps) planes = false;
            // All four vertices in the plane.
            for (int k = 1; k < 4; ++k)
                if (std::fabs(dot(f.n, f.v[k] - f.v[0])) > eps) planes = false;
        }
        check(planes, "box facet planes are exact support planes");
    }

    // ---- generic generators: face count and exact support planes ----
    {
        const Vec3 c{0.0, 0.0, 0.0};
        std::vector<Vec3> g;
        for (int i = 0; i < 6; ++i) g.push_back({rnd(), rnd(), rnd()});
        const std::vector<ZFace> mesh = sb2gui::zonotope_mesh(c, g);
        check((int)mesh.size() == 6 * 5, "generic m=6 mesh has m(m-1) facets");

        bool planes = true;
        for (const ZFace& f : mesh)
            if (std::fabs(dot(f.n, f.v[0] - c) - support(g, f.n)) > eps) planes = false;
        check(planes, "generic facet planes are exact support planes");

        // Every vertex of the zonotope (all sign words) is on or inside.
        bool inside = true;
        for (int w = 0; w < (1 << 6); ++w) {
            Vec3 p = c;
            for (int i = 0; i < 6; ++i) p = p + ((w >> i) & 1 ? 1.0 : -1.0) * g[i];
            if (outside_by(mesh, c, p) > eps) inside = false;
        }
        check(inside, "all 2^m sign vertices are enclosed by the mesh");
    }

    // ---- parallel merge is exact ----
    {
        const Vec3 a{0.3, 0.4, -0.2};
        const std::vector<Vec3> g = {a, -2.0 * a, {0.1, -0.5, 0.7}};
        const std::vector<Vec3> red = sb2gui::reduce_generators(g, 99);
        check(red.size() == 2, "parallel generators merge to one");
        bool same = true;
        for (int t = 0; t < 200; ++t) {
            const Vec3 n{rnd(), rnd(), rnd()};
            if (std::fabs(support(g, n) - support(red, n)) > eps) same = false;
        }
        check(same, "parallel merge keeps the exact support function");
    }

    // ---- order reduction is conservative ----
    {
        std::vector<Vec3> g;
        for (int i = 0; i < 24; ++i) g.push_back({rnd(), rnd(), rnd()});
        const std::vector<Vec3> red = sb2gui::reduce_generators(g, 8);
        check((int)red.size() <= 8, "reduction respects the generator budget");
        bool conservative = true;
        for (int t = 0; t < 500; ++t) {
            const Vec3 n{rnd(), rnd(), rnd()};
            if (support(red, n) < support(g, n) - eps) conservative = false;
        }
        check(conservative, "reduction never shrinks the zonotope (support only grows)");
    }

    // ---- projected silhouette encloses the projected mesh ----
    {
        const Vec3 c{0.0, 0.0, 0.0};
        std::vector<Vec3> g;
        for (int i = 0; i < 5; ++i) g.push_back({rnd(), rnd(), rnd()});
        const std::vector<ZFace> mesh = sb2gui::zonotope_mesh(c, g);

        // An arbitrary orthographic view: right/up from a tilted basis.
        const Vec3 r0{0.8, 0.6, 0.0};
        const Vec3 u0{-0.36, 0.48, 0.8};
        std::vector<Vec2> pg;
        for (const Vec3& e : g) pg.push_back({dot(e, r0), dot(e, u0)});
        const std::vector<Vec2> sil = sb2gui::zonotope_boundary({0.0, 0.0}, pg);
        check(sil.size() == 2 * pg.size(), "silhouette has 2m vertices");

        // Inside test against the CCW silhouette, with a small slack.
        bool enclosed = true;
        for (const ZFace& f : mesh)
            for (int k = 0; k < 4; ++k) {
                const double px = dot(f.v[k] - c, r0), py = dot(f.v[k] - c, u0);
                for (size_t i = 0; i < sil.size(); ++i) {
                    const Vec2& s0 = sil[i];
                    const Vec2& s1 = sil[(i + 1) % sil.size()];
                    if ((s1.x - s0.x) * (py - s0.y) - (s1.y - s0.y) * (px - s0.x) < -1e-9)
                        enclosed = false;
                }
            }
        check(enclosed, "projected mesh vertices lie inside the projected silhouette");
    }

    std::printf(failures == 0 ? "\nall zonotope3d checks passed\n"
                              : "\n%d zonotope3d check(s) FAILED\n", failures);
    return failures == 0 ? 0 : 1;
}
