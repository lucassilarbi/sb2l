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

#include <iostream>
#include <vector>

int main()
{
    // Parameters. The curve is non-rational, so it takes no weight: the
    // rational basis, and the weights that go with it, are in the advanced
    // example.
    int p(3);
    int nCP(5);
    sb2l::CurveType ct = sb2l::CurveType::UNIFORM_NONRATIONAL;
    sb2l::Form f = sb2l::Form::TAYLOR;
    sb2l::ParameterSet ps = sb2l::ParameterSet::IR;
    int d(100);
    int t(1); // -1: Automatic

    vibes::newFigure("BSPLINE");

    // Example of control boxes
    ibex::IntervalVector Px = ibex::IntervalVector(5);
    Px[0] = ibex::Interval(0, 1);
    Px[1] = ibex::Interval(50, 60);
    Px[2] = ibex::Interval(100, 101);
    Px[3] = ibex::Interval(130, 131);
    Px[4] = ibex::Interval(180, 181);
    ibex::IntervalVector Py = ibex::IntervalVector(5);
    Py[0] = ibex::Interval(60, 61);
    Py[1] = ibex::Interval(100, 120);
    Py[2] = ibex::Interval(60, 61);
    Py[3] = ibex::Interval(22, 23);
    Py[4] = ibex::Interval(30, 31);
    std::vector<ibex::IntervalVector> P{Px, Py};

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

    // B-spline generation
    sb2l::SB2 My_bspline(p, nCP, ct, f, ps, d, t);

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
}
