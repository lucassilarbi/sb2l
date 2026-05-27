/**
 * \date 03/2026
 * \author CEA/DRT/LIST/DIASI/SRI/LCSR
 * \author L.Si Larbi
 *
 * @par Licence
 * Copyright © 2026 CEA
 */

#include <sb2l.hpp>
#include <libqhullcpp/Qhull.h>
#include "vibes.h"

std::pair<std::vector<double>, std::vector<double>> compute_zonotope2D(ibex::Affine2Vector &v) // compute the zonotope resulting from two affine forms
{
    v.compact();
    int m(v.size()); // number of affine forms
    if (m > 2)
    {
        throw std::runtime_error("The number of affine form can not exceed 2");
        return std::pair<std::vector<double>, std::vector<double>>(std::vector<double>({}), std::vector<double>({}));
    }
    // compute the center points (also works for m>2)
    std::vector<double> c(m, 0); // center points
    for (int i = 0; i < m; i++)
    {
        c[i] = v[i].val(0);
    }
    // compute genrator vectors (also works for m>2)
    int G_size(v[0].size());
    for (int i = 1; i < m; i++)
    {
        if (v[i].size() > G_size)
            G_size = v[i].size();
    }
    // std::vector<std::vector<double>> G(G_size, std::vector<double>(m, 0));
    std::vector<std::vector<double>> uG(G_size, std::vector<double>(m, 0)); // uncompact generator vector
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < G_size; j++)
        {
            uG[j][i] = v[i].val(j + 1);
        }
    }
    // compact the generator vector (from here on, the code only works in 2D)
    std::vector<std::vector<double>> G({}); // compact generator vector
    for (int j = 0; j < G_size; j++)
    {
        if (uG[j][0] != 0 || uG[j][1] != 0)
        {
            G.push_back(std::vector<double>({uG[j][0], uG[j][1]}));
        }
    }
    G_size = G.size();
    // compute zontope points
    std::vector<std::vector<double>> vertices; // vertices of the zonotopes
    for (int b = 0; b < (1 << G_size); b++)    // 1 << G_size = 2**(G_size)
    {
        vertices.push_back(std::vector<double>(m, 0));
        for (int i = 0; i < m; i++)
        {
            vertices[vertices.size() - 1][i] = v[i].val(0) + [&]()
            {
                double r(0.0);
                for (int j = 0; j < G_size; j++)
                {
                    r += ((b >> j) & 1 ? 1 : -1) * G[j][i]; // b>>j return the jth bit of b. &1 only save this bit. b>>j)&1 with b from 0 to 2**(G_size) represent the binary system evolution with each bit of b identified and usable.
                }
                return r;
            }();
        }
    }
    // report
    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "maximum affine noises per affine form: " << ibex::AF_fAFFullI::getAffineNoiseNumber() << std::endl;
    // std::cout << "affine tolerance of the compaction: " << ibex::AF_fAFFullI::getAffineTolerance() << std::endl;
    // std::cout << "dim: " << m << std::endl;
    // std::cout << "number of generator vector: " << G_size << std::endl;
    // std::cout << "number of verticies: " << vertices.size() << std::endl;
    // flatten points (from here on, the code only works in 2D)
    std::vector<double> pts;
    for (const auto &p : vertices)
    {
        pts.push_back(p[0]);
        pts.push_back(p[1]);
    }
    // convex hull
    orgQhull::Qhull qh;
    std::vector<std::pair<double, double>> zonotope({});
    if (vertices.size() > 2)
    {
        qh.runQhull("", 2, vertices.size(), pts.data(), "Qt");
        // extract hull points
        for (const orgQhull::QhullVertex &v : qh.vertexList())
        {
            const double *p = v.point().coordinates();
            zonotope.emplace_back(p[0], p[1]);
        }
        // sort around center c
        std::sort(zonotope.begin(), zonotope.end(), [&](const auto &p1, const auto &p2)
                  { return std::atan2(p1.second - c[1], p1.first - c[0]) < std::atan2(p2.second - c[1], p2.first - c[0]); }); // return true if p1 is before p2 around the center (before in a counter clockwise sens)
    }
    else
    {
        for (const auto &p : vertices)
        {
            zonotope.emplace_back(p[0], p[1]);
        }
    }
    // to pair
    std::vector<double> X({}), Y({});
    for (const auto &p : zonotope)
    {
        X.push_back(p.first);
        Y.push_back(p.second);
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

    ibex::AF_fAFFullI::setAffineNoiseNumber(10); // need to be lower than 15 because in 2D, the worth case involving two affine forms with 15 independant affine noises will generate 30 generator vectors and then 2^30 vertices -> close to the signed int (required by qHull) numerical limit
    ibex::AF_fAFFullI::setAffineTolerance(2.5e-15);

    vibes::newFigure("BSPLINE");

    // Example of control points
    //================================================== CONTROL POINTS EXAMPLE 1 ==========================================================//
    double r(0.01); // radius
    ibex::IntervalVector Px = ibex::IntervalVector(12);
    Px[0] = ibex::Interval(1.0 - r, 1.0 + r);
    Px[1] = ibex::Interval(0.5 - r, 0.5 + r);
    Px[2] = ibex::Interval(0.1 - r * 10, 0.1 + r * 10);
    Px[3] = ibex::Interval(0.4 - r, 0.4 + r);
    Px[4] = ibex::Interval(1.9 - r * 20, 1.9 + r * 20);
    Px[5] = ibex::Interval(0.8 - r, 0.8 + r);
    Px[6] = ibex::Interval(1.9 - r, 1.9 + r);
    Px[7] = ibex::Interval(2.8 - r * 15, 2.8 + r * 15);
    Px[8] = ibex::Interval(1.4 - r, 1.4 + r);
    Px[9] = ibex::Interval(2.2 - r, 2.2 + r);
    Px[10] = ibex::Interval(2.5 - r, 2.5 + r);
    Px[11] = ibex::Interval(2.9 - r * 7, 2.9 + r * 7);
    ibex::IntervalVector Py = ibex::IntervalVector(12);
    Py[0] = ibex::Interval(1.0 - r, 1.0 + r);
    Py[1] = ibex::Interval(1.5 - r, 1.5 + r);
    Py[2] = ibex::Interval(0.9 - r * 13, 0.9 + r * 13);
    Py[3] = ibex::Interval(0.1 - r, 0.1 + r);
    Py[4] = ibex::Interval(0.5 - r * 20, 0.5 + r * 20);
    Py[5] = ibex::Interval(1.9 - r, 1.9 + r);
    Py[6] = ibex::Interval(2.2 - r, 2.2 + r);
    Py[7] = ibex::Interval(1.5 - r * 5, 1.5 + r * 5);
    Py[8] = ibex::Interval(1.8 - r, 1.8 + r);
    Py[9] = ibex::Interval(0.5 - r, 0.5 + r);
    Py[10] = ibex::Interval(1.0 - r, 1.0 + r);
    Py[11] = ibex::Interval(0.7 - r * 7, 0.7 + r * 7);
    std::vector<ibex::IntervalVector> P{Px, Py};
    std::vector<std::vector<double>> rP({{}, {}}); // real control points
    for (int i = 0; i < P[0].size(); i++)
    {
        rP[0].push_back(Px[i].mid());
        rP[1].push_back(Py[i].mid());
    }
    std::vector<ibex::Affine2Vector> aP(2, ibex::Affine2Vector(Px.size())); // affine control points
    for (int i = 0; i < P[0].size(); i++)
    {
        aP[0][i] = Px[i];
        aP[1][i] = Py[i];
        ibex::Affine2Vector av(2);
        av[0] = aP[0][i];
        av[1] = aP[1][i];
        std::pair<std::vector<double>, std::vector<double>> zonotope(compute_zonotope2D(av));
        vibes::drawPolygon(zonotope.first, zonotope.second, "grey[grey]");
    }
    //================================================== CONTROL POINTS EXAMPLE 2 ==========================================================//
    // int reserved(2); // reserved allow us to reserve some epsilon for control points. for example, if reserved=2, each control points is suppose to be defind with at most 2 affine form. 
    //                  // moreover, it is really important to correcly setup epsilon numbers (the first element of std::pair<int,double>(1, 0.5)) to be consistent with reserved.
    // ibex::Affine2Vector aPx = ibex::Affine2Vector(12, ibex::Affine2(0.0));
    // for(int i=0;i<reserved;i++) aPx[0]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[0].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(1, 0.01), std::pair<int,double>(2, -0.01), std::pair<int,double>(25, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[1]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[1].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(3, 0.01), std::pair<int,double>(4, -0.01), std::pair<int,double>(26, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[2]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[2].initialize(0.1, std::list<std::pair<int,double>>({std::pair<int,double>(5, 0.01), std::pair<int,double>(6, -0.01), std::pair<int,double>(27, 0.0), std::pair<int,double>(39, 0.1)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[3]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[3].initialize(0.4, std::list<std::pair<int,double>>({std::pair<int,double>(7, 0.01), std::pair<int,double>(8, -0.01), std::pair<int,double>(28, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[4]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[4].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(9, 0.1), std::pair<int,double>(10, -0.05), std::pair<int,double>(29, 0.05)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[5]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[5].initialize(0.8, std::list<std::pair<int,double>>({std::pair<int,double>(11, 0.01), std::pair<int,double>(12, -0.01), std::pair<int,double>(30, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[6]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[6].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(13, 0.01), std::pair<int,double>(14, -0.01), std::pair<int,double>(31, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[7]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[7].initialize(2.8, std::list<std::pair<int,double>>({std::pair<int,double>(15, -0.05), std::pair<int,double>(16, -0.1), std::pair<int,double>(32, 0.15), std::pair<int,double>(44, 0.1)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[8]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[8].initialize(1.4, std::list<std::pair<int,double>>({std::pair<int,double>(17, 0.01), std::pair<int,double>(18, -0.01), std::pair<int,double>(33, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[9]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[9].initialize(2.2, std::list<std::pair<int,double>>({std::pair<int,double>(19, 0.01), std::pair<int,double>(20, -0.01), std::pair<int,double>(34, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[10]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[10].initialize(2.5, std::list<std::pair<int,double>>({std::pair<int,double>(21, 0.01), std::pair<int,double>(22, -0.01), std::pair<int,double>(35, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPx[11]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPx[11].initialize(2.9, std::list<std::pair<int,double>>({std::pair<int,double>(23, 0.07), std::pair<int,double>(24, -0.01), std::pair<int,double>(36, -0.1)}), ibex::Interval(0.0));
    // ibex::Affine2Vector aPy = ibex::Affine2Vector(12, ibex::Affine2(0.0));
    // for(int i=0;i<reserved;i++) aPy[0]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[0].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(1, 0.01), std::pair<int,double>(2, 0.01), std::pair<int,double>(37, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[1]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[1].initialize(1.5, std::list<std::pair<int,double>>({std::pair<int,double>(3, 0.01), std::pair<int,double>(4, 0.01), std::pair<int,double>(38, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[2]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[2].initialize(0.9, std::list<std::pair<int,double>>({std::pair<int,double>(5, 0.01), std::pair<int,double>(6, 0.01), std::pair<int,double>(39, 0.3)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[3]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[3].initialize(0.1, std::list<std::pair<int,double>>({std::pair<int,double>(7, 0.01), std::pair<int,double>(8, 0.01), std::pair<int,double>(40, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[4]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[4].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(9, 0.1), std::pair<int,double>(10, 0.1), std::pair<int,double>(41, 0.05)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[5]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[5].initialize(1.9, std::list<std::pair<int,double>>({std::pair<int,double>(11, 0.01), std::pair<int,double>(12, 0.01), std::pair<int,double>(42, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[6]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[6].initialize(2.2, std::list<std::pair<int,double>>({std::pair<int,double>(13, 0.01), std::pair<int,double>(14, 0.01), std::pair<int,double>(43, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[7]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[7].initialize(1.5, std::list<std::pair<int,double>>({std::pair<int,double>(15, 0.05), std::pair<int,double>(16, 0.15), std::pair<int,double>(44, 0.1), std::pair<int,double>(32, 0.03)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[8]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[8].initialize(1.8, std::list<std::pair<int,double>>({std::pair<int,double>(17, 0.01), std::pair<int,double>(18, 0.01), std::pair<int,double>(45, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[9]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[9].initialize(0.5, std::list<std::pair<int,double>>({std::pair<int,double>(19, 0.01), std::pair<int,double>(20, 0.01), std::pair<int,double>(46, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[10]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[10].initialize(1.0, std::list<std::pair<int,double>>({std::pair<int,double>(21, 0.01), std::pair<int,double>(22, 0.01), std::pair<int,double>(47, 0.0)}), ibex::Interval(0.0));
    // for(int i=0;i<reserved;i++) aPy[11]+=ibex::Affine2(ibex::Interval(-1,1));
    // aPy[11].initialize(0.7, std::list<std::pair<int,double>>({std::pair<int,double>(23, 0.01), std::pair<int,double>(24, 0.03), std::pair<int,double>(48, 0.0)}), ibex::Interval(0.0));
    // std::vector<ibex::Affine2Vector> aP(2, ibex::Affine2Vector(aPx.size())); // affine control points
    // for (int i=0; i<aPx.size(); i++)
    // {
    //     aP[0][i] = aPx[i];
    //     aP[1][i] = aPy[i];
    //     ibex::Affine2Vector av(2);
    //     av[0]=aP[0][i]; av[1]=aP[1][i];
    //     std::pair<std::vector<double>,std::vector<double>> zonotope(compute_zonotope2D(av));
    //     vibes::drawPolygon(zonotope.first, zonotope.second, "grey[grey]");
    // }
    // std::vector<ibex::IntervalVector> P(2, ibex::IntervalVector(aPx.size()));
    // for (int i=0; i<P[0].size(); i++)
    // {
    //     P[0][i] = aPx[i].itv();
    //     P[1][i] = aPy[i].itv();
    // }
    // std::vector<std::vector<double>> rP({{}, {}}); // real control points
    // for (int i = 0; i < P[0].size(); i++)
    // {
    //     rP[0].push_back(aPx[i].mid());
    //     rP[1].push_back(aPy[i].mid());
    // }
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
    vibes::selectFigure("BASIS");
    std::vector<std::vector<ibex::Interval>> idu = My_bspline.get_idu();
    std::vector<std::vector<std::vector<ibex::Interval>>> ieBf = My_bspline.get_ieBf();
    vibes::newFigure("BASIS");
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