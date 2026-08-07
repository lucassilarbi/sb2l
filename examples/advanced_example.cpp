/**
 * \date 03/2026
 * \author CEA/DRT/LIST/DIASI/SRI/LCSR
 * \author L.Si Larbi
 *
 * @par Licence
 * Copyright © 2026 CEA
 */

#include <sb2l.hpp>
#include "vibes.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <stdexcept>
#include <utility>
#include <vector>

std::pair<std::vector<double>, std::vector<double>> compute_zonotope2D(ibex::Affine2Vector &v) // compute the zonotope resulting from two affine forms
{
    v.compact(); // shortens the list of noise terms, nothing else
    if (v.size() != 2)
        throw std::runtime_error("The zonotope is drawn in the plane: exactly two affine forms are expected");
    // The center and the generators of the element, read through the library:
    // sb2l::zonotope_of is the one function which knows how a set comes out of
    // an affine vector, and its documentation says why reading it by hand,
    // through size() and val(i), loses the generators without a word.
    const sb2l::Zonotope z = sb2l::zonotope_of(v);
    const std::vector<double> c = z.center;
    std::vector<std::vector<double>> G;
    for (int k = 0; k < z.m; k++)
        G.push_back({z.gen(k, 0), z.gen(k, 1)});
    // Boundary of { c + sum s_i*G_i : s_i in [-1,1] } by the O(m log m)
    // Minkowski walk: exactly 2*m vertices, already counter-clockwise.
    // Replaces the 2^m vertex enumeration + convex hull.

    // Fold every generator into the upper half-plane [0,pi); the zonotope is
    // symmetric about c, so flipping a generator leaves the set unchanged.
    std::vector<std::vector<double>> g;
    for (const auto &e : G)
    {
        if (e[0] == 0 && e[1] == 0)
            continue;
        if (e[1] < 0 || (e[1] == 0 && e[0] < 0))
            g.push_back({-e[0], -e[1]});
        else
            g.push_back({e[0], e[1]});
    }
    // Sort by angle. Each generator gets a number computed once and the sort
    // compares numbers: comparing two generators by the sign of their cross
    // product can be intransitive in rounding on nearly parallel ones, which
    // std::sort does not tolerate. Over the half-plane above, x / (|x| + y)
    // decreases from 1 to -1 as the angle goes from 0 to pi.
    std::vector<std::pair<double, std::vector<double>>> keyed;
    for (const auto &e : g)
        keyed.push_back(std::make_pair(e[0] / (std::fabs(e[0]) + e[1]), e));
    std::sort(keyed.begin(), keyed.end(),
              [](const std::pair<double, std::vector<double>> &a, const std::pair<double, std::vector<double>> &b)
              { return a.first > b.first; });
    for (size_t i = 0; i < g.size(); i++)
        g[i] = keyed[i].second;

    // Start at the bottom vertex c - sum(g), then trace both point-symmetric halves.
    std::vector<double> X({}), Y({});
    double px(c[0]), py(c[1]);
    for (const auto &e : g)
    {
        px -= e[0];
        py -= e[1];
    }
    for (const auto &e : g) // upper chain
    {
        X.push_back(px);
        Y.push_back(py);
        px += 2 * e[0];
        py += 2 * e[1];
    }
    for (const auto &e : g) // lower chain (point-symmetric)
    {
        X.push_back(px);
        Y.push_back(py);
        px -= 2 * e[0];
        py -= 2 * e[1];
    }
    if (g.empty()) // degenerate: single point
    {
        X.push_back(c[0]);
        Y.push_back(c[1]);
    }
    return std::pair<std::vector<double>, std::vector<double>>(X, Y);
}

int main()
{
    // Parameters
    int p(3);
    int nCP(12);
    sb2l::CurveType ct = sb2l::CurveType::CLAMPED_RATIONAL;
    sb2l::Form f = sb2l::Form::TAYLOR;
    sb2l::ParameterSet ps = sb2l::ParameterSet::Z;
    int d(10);
    int t(-1); // -1: Automatic
    std::vector<SymEngine::Expression> rw({SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(2),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(3) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(4) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1),
                                           SymEngine::Expression(1) / SymEngine::Expression(1)});

    // An affine form keeps this many noise symbols and merges the smallest
    // ones into its error term beyond that, so this is what decides how much
    // of the shape of a control zonotope survives the evaluation. One symbol
    // per control point and per coordinate, plus room for the ones the
    // evaluation introduces itself. The boundary walk below is linearithmic in
    // the number of symbols, so raising it costs little.
    sb2l::SB2::set_affine_noise_number(2 * nCP + 16);
    ibex::AF_fAFFullI::setAffineTolerance(2.5e-15);

    vibes::newFigure("BSPLINE");

    // Example of control points
    //================================================== CONTROL POINTS EXAMPLE 1 ==========================================================//
    // double r(0.01); // radius
    // ibex::IntervalVector Px = ibex::IntervalVector(12);
    // Px[0] = ibex::Interval(1.0 - r, 1.0 + r);
    // Px[1] = ibex::Interval(0.5 - r, 0.5 + r);
    // Px[2] = ibex::Interval(0.1 - r * 10, 0.1 + r * 10);
    // Px[3] = ibex::Interval(0.4 - r, 0.4 + r);
    // Px[4] = ibex::Interval(1.9 - r * 20, 1.9 + r * 20);
    // Px[5] = ibex::Interval(0.8 - r, 0.8 + r);
    // Px[6] = ibex::Interval(1.9 - r, 1.9 + r);
    // Px[7] = ibex::Interval(2.8 - r * 15, 2.8 + r * 15);
    // Px[8] = ibex::Interval(1.4 - r, 1.4 + r);
    // Px[9] = ibex::Interval(2.2 - r, 2.2 + r);
    // Px[10] = ibex::Interval(2.5 - r, 2.5 + r);
    // Px[11] = ibex::Interval(2.9 - r * 7, 2.9 + r * 7);
    // ibex::IntervalVector Py = ibex::IntervalVector(12);
    // Py[0] = ibex::Interval(1.0 - r, 1.0 + r);
    // Py[1] = ibex::Interval(1.5 - r, 1.5 + r);
    // Py[2] = ibex::Interval(0.9 - r * 13, 0.9 + r * 13);
    // Py[3] = ibex::Interval(0.1 - r, 0.1 + r);
    // Py[4] = ibex::Interval(0.5 - r * 20, 0.5 + r * 20);
    // Py[5] = ibex::Interval(1.9 - r, 1.9 + r);
    // Py[6] = ibex::Interval(2.2 - r, 2.2 + r);
    // Py[7] = ibex::Interval(1.5 - r * 5, 1.5 + r * 5);
    // Py[8] = ibex::Interval(1.8 - r, 1.8 + r);
    // Py[9] = ibex::Interval(0.5 - r, 0.5 + r);
    // Py[10] = ibex::Interval(1.0 - r, 1.0 + r);
    // Py[11] = ibex::Interval(0.7 - r * 7, 0.7 + r * 7);
    // std::vector<ibex::IntervalVector> P{Px, Py};
    // std::vector<std::vector<double>> rP({{}, {}}); // real control points
    // for (int i = 0; i < P[0].size(); i++)
    // {
    //     rP[0].push_back(Px[i].mid());
    //     rP[1].push_back(Py[i].mid());
    // }
    // std::vector<ibex::Affine2Vector> aP(2, ibex::Affine2Vector(Px.size())); // affine control points
    // for (int i = 0; i < P[0].size(); i++)
    // {
    //     aP[0][i] = Px[i];
    //     aP[1][i] = Py[i];
    //     ibex::Affine2Vector av(2);
    //     av[0] = aP[0][i];
    //     av[1] = aP[1][i];
    //     std::pair<std::vector<double>, std::vector<double>> zonotope(compute_zonotope2D(av));
    //     vibes::drawPolygon(zonotope.first, zonotope.second, "grey[grey]");
    // }
    //================================================== CONTROL POINTS EXAMPLE 2 ==========================================================//
    int reserved(2); // reserved allow us to reserve some epsilon for control points. for example, if reserved=2, each control points is suppose to be defind with at most 2 affine form. 
                     // moreover, it is really important to correcly setup epsilon numbers (the first element of std::pair<int,double>(1, 0.5)) to be consistent with reserved.
    ibex::Affine2Vector aPx = ibex::Affine2Vector(12, ibex::Affine2(0.0));
    for(int i=0;i<reserved;i++) aPx[0]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[0].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(1, 0.01), std::pair<int,double>(2, -0.01), std::pair<int,double>(25, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[1]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[1].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(3, 0.01), std::pair<int,double>(4, -0.01), std::pair<int,double>(26, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[2]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[2].initialize(0.1, std::list<std::pair<int,double>>({std::pair<int,double>(5, 0.01), std::pair<int,double>(6, -0.01), std::pair<int,double>(27, 0.0), std::pair<int,double>(39, 0.1)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[3]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[3].initialize(0.4, std::list<std::pair<int,double>>({std::pair<int,double>(7, 0.01), std::pair<int,double>(8, -0.01), std::pair<int,double>(28, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[4]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[4].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(9, 0.1), std::pair<int,double>(10, -0.05), std::pair<int,double>(29, 0.05)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[5]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[5].initialize(0.8, std::list<std::pair<int,double>>({std::pair<int,double>(11, 0.01), std::pair<int,double>(12, -0.01), std::pair<int,double>(30, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[6]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[6].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(13, 0.01), std::pair<int,double>(14, -0.01), std::pair<int,double>(31, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[7]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[7].initialize(2.8, std::list<std::pair<int,double>>({std::pair<int,double>(15, -0.05), std::pair<int,double>(16, -0.1), std::pair<int,double>(32, 0.15), std::pair<int,double>(44, 0.1)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[8]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[8].initialize(1.4, std::list<std::pair<int,double>>({std::pair<int,double>(17, 0.01), std::pair<int,double>(18, -0.01), std::pair<int,double>(33, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[9]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[9].initialize(2.2, std::list<std::pair<int,double>>({std::pair<int,double>(19, 0.01), std::pair<int,double>(20, -0.01), std::pair<int,double>(34, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[10]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[10].initialize(2.5, std::list<std::pair<int,double>>({std::pair<int,double>(21, 0.01), std::pair<int,double>(22, -0.01), std::pair<int,double>(35, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPx[11]+=ibex::Affine2(ibex::Interval(-1,1));
    aPx[11].initialize(2.9, std::list<std::pair<int,double>>({std::pair<int,double>(23, 0.07), std::pair<int,double>(24, -0.01), std::pair<int,double>(36, -0.1)}), ibex::Interval(0.0));
    ibex::Affine2Vector aPy = ibex::Affine2Vector(12, ibex::Affine2(0.0));
    for(int i=0;i<reserved;i++) aPy[0]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[0].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(1, 0.01), std::pair<int,double>(2, 0.01), std::pair<int,double>(37, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[1]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[1].initialize(1.5, std::list<std::pair<int,double>>({std::pair<int,double>(3, 0.01), std::pair<int,double>(4, 0.01), std::pair<int,double>(38, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[2]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[2].initialize(0.9, std::list<std::pair<int,double>>({std::pair<int,double>(5, 0.01), std::pair<int,double>(6, 0.01), std::pair<int,double>(39, 0.3)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[3]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[3].initialize(0.1, std::list<std::pair<int,double>>({std::pair<int,double>(7, 0.01), std::pair<int,double>(8, 0.01), std::pair<int,double>(40, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[4]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[4].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(9, 0.1), std::pair<int,double>(10, 0.1), std::pair<int,double>(41, 0.05)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[5]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[5].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(11, 0.01), std::pair<int,double>(12, 0.01), std::pair<int,double>(42, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[6]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[6].initialize(2.2, std::list<std::pair<int,double>>({std::pair<int,double>(13, 0.01), std::pair<int,double>(14, 0.01), std::pair<int,double>(43, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[7]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[7].initialize(1.5, std::list<std::pair<int,double>>({std::pair<int,double>(15, 0.05), std::pair<int,double>(16, 0.15), std::pair<int,double>(44, 0.1), std::pair<int,double>(32, 0.03)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[8]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[8].initialize(1.8, std::list<std::pair<int,double>>({std::pair<int,double>(17, 0.01), std::pair<int,double>(18, 0.01), std::pair<int,double>(45, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[9]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[9].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(19, 0.01), std::pair<int,double>(20, 0.01), std::pair<int,double>(46, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[10]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[10].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(21, 0.01), std::pair<int,double>(22, 0.01), std::pair<int,double>(47, 0.0)}), ibex::Interval(0.0));
    for(int i=0;i<reserved;i++) aPy[11]+=ibex::Affine2(ibex::Interval(-1,1));
    aPy[11].initialize(0.7, std::list<std::pair<int,double>>({std::pair<int,double>(23, 0.01), std::pair<int,double>(24, 0.03), std::pair<int,double>(48, 0.0)}), ibex::Interval(0.0));
    std::vector<ibex::Affine2Vector> aP(2, ibex::Affine2Vector(aPx.size())); // affine control points
    for (int i=0; i<aPx.size(); i++)
    {
        aP[0][i] = aPx[i];
        aP[1][i] = aPy[i];
        ibex::Affine2Vector av(2);
        av[0]=aP[0][i]; av[1]=aP[1][i];
        std::pair<std::vector<double>,std::vector<double>> zonotope(compute_zonotope2D(av));
        vibes::drawPolygon(zonotope.first, zonotope.second, "grey[grey]");
    }
    std::vector<ibex::IntervalVector> P(2, ibex::IntervalVector(aPx.size()));
    for (int i=0; i<P[0].size(); i++)
    {
        P[0][i] = aPx[i].itv();
        P[1][i] = aPy[i].itv();
    }
    std::vector<std::vector<double>> rP({{}, {}}); // real control points
    for (int i = 0; i < P[0].size(); i++)
    {
        rP[0].push_back(aPx[i].mid());
        rP[1].push_back(aPy[i].mid());
    }
    //============================================================================================================//

    std::cout << "Interval B-spline parameters:" << std::endl;
    std::cout << "  - Degree: " << p << std::endl;
    std::cout << "  - Number of control points: " << nCP << std::endl;
    std::cout << "  - Curve type: " << ct << std::endl;
    std::cout << "  - Equation form: " << f << std::endl;
    std::cout << "  - Parameter set: " << ps << std::endl;
    std::cout << "  - Number of evaluation per segment: " << d << std::endl;
    if (f == sb2l::Form::TAYLOR)
    {
        if (t == -1)
        {
            std::cout << "  - Taylor order: automatic" << std::endl;
        }
        else
        {
            std::cout << "  - Taylor order: " << t << std::endl;
        }
    }
    std::cout << "  - Rational weight vector: [";
    for (size_t i = 0; i < rw.size(); i++)
    {
        if (i != rw.size() - 1)
            std::cout << rw[i] << ", ";
        else
            std::cout << rw[i] << "]";
    }
    std::cout << "" << "" << std::endl;

    // B-spline generation
    sb2l::SB2 My_bspline(p, nCP, ct, f, ps, d, t, rw);

    // Eval boxes
    std::vector<std::vector<ibex::IntervalVector>> boxes = My_bspline.eval_box(P);
    // Plot the interval B-spline using Vibes
    vibes::selectFigure("BSPLINE");
    for (unsigned int s = 0; s < boxes.size(); s++)
    {
        for (unsigned int du = 0; du < boxes[s].size(); du++)
        {
            vibes::drawBox(boxes[s][du][0].lb(), boxes[s][du][0].ub(), boxes[s][du][1].lb(), boxes[s][du][1].ub(), "black[]");
        }
    }
    // Plot the basis using Vibes
    std::vector<std::vector<ibex::Interval>> idu = My_bspline.get_idu();
    std::vector<std::vector<std::vector<ibex::Interval>>> ieBf = My_bspline.get_ieBf();
    vibes::newFigure("BASIS"); // created before it is drawn into
    for (unsigned int s = 0; s < ieBf.size(); s++)
    {
        for (unsigned int i = 0; i < ieBf[s].size(); i++)
        {
            for (unsigned int du = 0; du < ieBf[s][i].size(); du++)
            {
                vibes::drawBox(idu[s][du].lb(), idu[s][du].ub(), ieBf[s][i][du].lb(), ieBf[s][i][du].ub(), "black[]");
            }
        }
    }

    // Affine B-spline evaluation
    std::vector<std::vector<ibex::Affine2Vector>> zonotopes = My_bspline.eval_zonotope(aP);
    // Plot the affine B-spline using Vibes
    vibes::selectFigure("BSPLINE");
    for (unsigned int s = 0; s < zonotopes.size(); s++)
    {
        for (unsigned int du = 0; du < zonotopes[s].size(); du++)
        {
            std::pair<std::vector<double>, std::vector<double>> zonotope(compute_zonotope2D(zonotopes[s][du]));
            vibes::drawPolygon(zonotope.first, zonotope.second, "red[]");
        }
    }
    // Plot the basis using Vibes
    vibes::selectFigure("BASIS");
    std::vector<std::vector<ibex::Affine2>> adu = My_bspline.get_adu();
    std::vector<std::vector<std::vector<ibex::Affine2>>> aeBf = My_bspline.get_aeBf();
    for (unsigned int s = 0; s < aeBf.size(); s++)
    {
        for (unsigned int i = 0; i < aeBf[s].size(); i++)
        {
            for (unsigned int du = 0; du < aeBf[s][i].size(); du++)
            {
                ibex::Affine2Vector av(2);
                av[0] = adu[s][du];
                av[1] = aeBf[s][i][du];
                std::pair<std::vector<double>, std::vector<double>> zonotope(compute_zonotope2D(av));
                vibes::drawPolygon(zonotope.first, zonotope.second, "red[]");
            }
        }
    }

    // Eval points
    std::vector<std::vector<std::vector<double>>> points = My_bspline.eval_point(rP);
    // Plot the interval B-spline using Vibes
    vibes::selectFigure("BSPLINE");
    for (unsigned int s = 0; s < points.size(); s++)
    {
        for (unsigned int du = 0; du < points[s].size(); du++)
        {
            vibes::drawCircle(points[s][du][0], points[s][du][1], 0.01, "blue[]");
        }
    }
    // Plot the basis using Vibes
    std::vector<std::vector<ibex::Interval>> rdu = My_bspline.get_rdu();
    std::vector<std::vector<std::vector<double>>> reBf = My_bspline.get_reBf();
    vibes::selectFigure("BASIS");
    // over x
    for (unsigned int s = 0; s < reBf.size(); s++)
    {
        for (unsigned int i = 0; i < reBf[s].size(); i++)
        {
            for (unsigned int du = 0; du < reBf[s][i].size(); du++)
            {
                vibes::drawCircle(rdu[s][du].mid(), reBf[s][i][du], 0.01, "blue[]");
            }
        }
    }
}