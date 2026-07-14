/**
 * \date 03/2026
 * \author CEA/DRT/LIST/DIASI/SRI/LCSR
 * \author L.Si Larbi
 *
 * @par Licence
 * Copyright © 2026 CEA
 */

#include <sb2l.hpp>

#include <algorithm>

#define LWE true // limit wrapping effect when evaluating the B-spline using Z for the parameter set and IR for control points but increase the computation time.

namespace sb2l {

SB2::SB2(const int p,
           const int nCP,
           const CurveType ct,
           const Form f,
           const ParameterSet ps,
           const int d,
           const int t,
           const std::vector<SymEngine::Expression> W)
    : p_(p),
      n_(nCP - 1),
      ct_(ct),
      f_(f),
      ps_(ps),
      d_(d),
      t_(t),
      W_(W),
      k_(nCP + p),
      nS_(nCP - p)
{
    if (p_ >= n_ + 1)
        throw std::runtime_error("The degree of the curve p must be lower than the number of control points n+1");
    if (f_ == Form::TAYLOR)
    {
        if (t_ == -1) // by default, the taylor order is chosen to minimize the wraping effect
        {
            if (p_ <= 12)
                t_ = p_ - 1; // numerical limit for the factorial
            else
                t_ = 11;
        }
        else if (t_ > 11)
            t_ = 11;
    }
    else
        t_ = -1;

    if (ct == CurveType::UNIFORM_RATIONAL || ct == CurveType::CLAMPED_RATIONAL)
    {
        if (W_.size() == 0) // the user does not provide rational weight. The default ones are therfore setup to 1
        {
            for (int i = 0; i < n_ + 1; i++)
            {
                W_.push_back(SymEngine::Expression(1));
            }
        }
        else if (W_.size() != n_ + 1)
            throw std::runtime_error("The chosen rational weight vector is wrong");
    }
    compute_U_();
    compute_N_();
    compute_iN_();
    switch (ps)
    {
    case ParameterSet::R:
        compute_rdu_();
        compute_reBf();
        break;
    case ParameterSet::IR:
        compute_idu_();
        compute_ieBf();
        break;
    case ParameterSet::Z:
        compute_adu_();
        compute_aeBf();
        break;
    default:
        throw std::runtime_error("Unknown ParameterSet");
    }
}
int SB2::get_p() { return p_; }
int SB2::get_n() { return n_; }
int SB2::get_nS() { return nS_; }
int SB2::get_d() { return d_; }
std::vector<std::vector<std::vector<double>>> SB2::get_reBf() 
{
    switch (ps_)
    {
    case ParameterSet::IR:
        reBf_ = std::vector<std::vector<std::vector<double>>>(nS_, std::vector<std::vector<double>>(n_ + 1, std::vector<double>(d_, 0.0)));
        for (int s = 0; s < nS_; s++)
        {
            for (int i = s; i < s + p_ + 1; i++)
            {
                for (int du = 0; du < d_; du++)
                {
                    reBf_[s][i][du] = ieBf_[s][i][du].mid();
                }
            }
        }
        break;
    case ParameterSet::Z:
        reBf_ = std::vector<std::vector<std::vector<double>>>(nS_, std::vector<std::vector<double>>(n_ + 1, std::vector<double>(d_, 0.0)));
        for (int s = 0; s < nS_; s++)
        {
            for (int i = s; i < s + p_ + 1; i++)
            {
                for (int du = 0; du < d_; du++)
                {
                    reBf_[s][i][du] = aeBf_[s][i][du].mid();
                }
            }
        }
        break;
    }
    return reBf_; 
}
std::vector<std::vector<std::vector<ibex::Interval>>> SB2::get_ieBf()
{
    switch (ps_)
    {
    case ParameterSet::R:
        throw std::runtime_error("Can not build an interval basis from a real one");
        break;
    case ParameterSet::Z:
        ieBf_ = std::vector<std::vector<std::vector<ibex::Interval>>>(nS_, std::vector<std::vector<ibex::Interval>>(n_ + 1, std::vector<ibex::Interval>(d_, ibex::Interval(0.0))));
        for (int s = 0; s < nS_; s++)
        {
            for (int i = s; i < s + p_ + 1; i++)
            {
                for (int du = 0; du < d_; du++)
                {
                    ieBf_[s][i][du] = aeBf_[s][i][du].itv();
                }
            }
        }
        break;
    }
    return ieBf_; 
}
std::vector<std::vector<std::vector<ibex::Affine2>>> SB2::get_aeBf()
{
    switch (ps_)
    {
    case ParameterSet::R:
        throw std::runtime_error("Can not build an affine basis from a real one");
        break;
    case ParameterSet::IR:
        throw std::runtime_error("Can not build an affine basis from an interval one");
        break;
    }
    return aeBf_;
}
std::vector<std::vector<std::vector<double>>> SB2::get_reBf_diff(const int order)
{
    auto to_ibex = [](std::string expr) { for (size_t pos = 0; (pos = expr.find("**", pos)) != std::string::npos; expr.replace(pos, 2, "^"), ++pos); return expr; };
    std::vector<std::vector<std::vector<double>>> reBf_diff(nS_, std::vector<std::vector<double>>(n_ + 1, std::vector<double>(d_, 0.0)));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            SymEngine::Expression expr = N_[0][s][i];
            for (int t = 0; t < order; t++) expr = expr.diff(u_);
            ibex::Function iN_diff("u", to_ibex(expr.get_basic()->__str__()).c_str());
            for (int du = 0; du < d_; du++) // evaluated on a single parameter value: neither the Taylor form nor its remainder is needed
            {
                const double u = (ps_ == ParameterSet::Z) ? adu_[s][du].mid() : ((ps_ == ParameterSet::IR) ? idu_[s][du].mid() : rdu_[s][du].mid());
                reBf_diff[s][i][du] = iN_diff.eval(ibex::IntervalVector(1, u)).mid();
            }
        }
    }
    return reBf_diff;
}
std::vector<std::vector<std::vector<ibex::Interval>>> SB2::get_ieBf_diff(const int order)
{
    if (ps_ == ParameterSet::R)
        throw std::runtime_error("Can not build an interval basis from a real one");
    auto to_ibex = [](std::string expr) { for (size_t pos = 0; (pos = expr.find("**", pos)) != std::string::npos; expr.replace(pos, 2, "^"), ++pos); return expr; };
    const int nT = (f_ == Form::TAYLOR) ? t_ + 1 : 0; // last derivative level used by the Taylor remainder
    std::vector<std::vector<std::vector<ibex::Interval>>> ieBf_diff(nS_, std::vector<std::vector<ibex::Interval>>(n_ + 1, std::vector<ibex::Interval>(d_, ibex::Interval(0.0))));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            std::vector<std::shared_ptr<ibex::Function>> iN_diff(nT + 1);
            SymEngine::Expression expr = N_[0][s][i];
            for (int t = 0; t < order; t++) expr = expr.diff(u_);
            for (int t = 0; t <= nT; t++)
            {
                iN_diff[t] = std::make_shared<ibex::Function>("u", to_ibex(expr.get_basic()->__str__()).c_str());
                expr = expr.diff(u_);
            }
            for (int du = 0; du < d_; du++)
            {
                if (f_ == Form::TAYLOR)
                {
                    if (ps_ == ParameterSet::Z)
                    {
                        ibex::Affine2 buffer(0);
                        for (int t = 0; t <= t_; t++)
                            buffer += iN_diff[t]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du].mid())) / factorial(t) * pow(adu_[s][du] - adu_[s][du].mid(), t);
                        buffer += iN_diff[t_ + 1]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du])) / factorial(t_ + 1) * pow(adu_[s][du] - adu_[s][du].mid(), t_ + 1);
                        ieBf_diff[s][i][du] = buffer.itv();
                    }
                    else
                    {
                        ibex::Interval buffer(0, 0);
                        for (int t = 0; t <= t_; t++)
                            buffer += iN_diff[t]->eval(ibex::IntervalVector(1, idu_[s][du].mid())) / factorial(t) * pow(idu_[s][du] - idu_[s][du].mid(), t);
                        buffer += iN_diff[t_ + 1]->eval(ibex::IntervalVector(1, idu_[s][du])) / factorial(t_ + 1) * pow(idu_[s][du] - idu_[s][du].mid(), t_ + 1);
                        ieBf_diff[s][i][du] = buffer;
                    }
                }
                else
                {
                    ieBf_diff[s][i][du] = (ps_ == ParameterSet::Z) ? iN_diff[0]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du])).itv()
                                                                    : iN_diff[0]->eval(ibex::IntervalVector(1, idu_[s][du]));
                }
            }
        }
    }
    return ieBf_diff;
}
std::vector<std::vector<std::vector<ibex::Affine2>>> SB2::get_aeBf_diff(const int order)
{
    if (ps_ != ParameterSet::Z)
        throw std::runtime_error("Can not build an affine basis derivative outside of the Z parameter set");
    auto to_ibex = [](std::string expr) { for (size_t pos = 0; (pos = expr.find("**", pos)) != std::string::npos; expr.replace(pos, 2, "^"), ++pos); return expr; };
    const int nT = (f_ == Form::TAYLOR) ? t_ + 1 : 0; // last derivative level used by the Taylor remainder
    std::vector<std::vector<std::vector<ibex::Affine2>>> aeBf_diff(nS_, std::vector<std::vector<ibex::Affine2>>(n_ + 1, std::vector<ibex::Affine2>(d_, ibex::Affine2(0.0))));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            std::vector<std::shared_ptr<ibex::Function>> iN_diff(nT + 1);
            SymEngine::Expression expr = N_[0][s][i];
            for (int t = 0; t < order; t++) expr = expr.diff(u_);
            for (int t = 0; t <= nT; t++)
            {
                iN_diff[t] = std::make_shared<ibex::Function>("u", to_ibex(expr.get_basic()->__str__()).c_str());
                expr = expr.diff(u_);
            }
            for (int du = 0; du < d_; du++)
            {
                if (f_ == Form::TAYLOR)
                {
                    ibex::Affine2 buffer(0);
                    for (int t = 0; t <= t_; t++)
                        buffer += iN_diff[t]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du].mid())) / factorial(t) * pow(adu_[s][du] - adu_[s][du].mid(), t);
                    buffer += iN_diff[t_ + 1]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du])) / factorial(t_ + 1) * pow(adu_[s][du] - adu_[s][du].mid(), t_ + 1);
                    aeBf_diff[s][i][du] = buffer;
                }
                else
                {
                    aeBf_diff[s][i][du] = iN_diff[0]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du]));
                }
            }
        }
    }
    return aeBf_diff;
}
std::vector<std::vector<ibex::Interval>> SB2::get_rdu()
{
    switch (ps_)
    {
    case ParameterSet::IR:
        rdu_ = std::vector<std::vector<ibex::Interval>>(nS_);
        for (int s = 0; s < nS_; s++)
            for (int du = 0; du < d_; du++)
                rdu_[s].push_back(ibex::Interval(idu_[s][du].mid()));
        break;
    case ParameterSet::Z:
        rdu_ = std::vector<std::vector<ibex::Interval>>(nS_);
        for (int s = 0; s < nS_; s++)
            for (int du = 0; du < d_; du++)
                rdu_[s].push_back(ibex::Interval(adu_[s][du].mid()));
        break;
    }
    return rdu_;
}
std::vector<std::vector<ibex::Interval>> SB2::get_idu()
{
    switch (ps_)
    {
    case ParameterSet::R:
        throw std::runtime_error("Can not build an interval decomposition from a real one");
        break;
    case ParameterSet::Z:
        idu_ = std::vector<std::vector<ibex::Interval>>(nS_);
        for (int s = 0; s < nS_; s++)
            for (int du = 0; du < d_; du++)
                idu_[s].push_back(adu_[s][du].itv());
        break;
    }
    return idu_;
}
std::vector<std::vector<ibex::Affine2>> SB2::get_adu()
{
    switch (ps_)
    {
    case ParameterSet::R:
        throw std::runtime_error("Can not build an affine decomposition from a real one");
        break;
    case ParameterSet::IR:
        throw std::runtime_error("Can not build an affine decomposition from an interval one");
        break;
    }
    return adu_;
}
void SB2::compute_U_()
{
    if (ct_ == CurveType::UNIFORM_RATIONAL || ct_ == CurveType::UNIFORM_NONRATIONAL)
    {
        U_ = {};
        for (int i = 0; i < k_ + 1; i++)
        {
            U_.push_back(i);
        }
    }
    else if (ct_ == CurveType::CLAMPED_RATIONAL || ct_ == CurveType::CLAMPED_NONRATIONAL)
    {
        U_ = {};
        for (int i = 0; i < k_ + 1; i++)
        {
            if (i < p_ + 1)
            {
                U_.push_back(p_);
            }
            else if (i >= k_ - p_)
            {
                U_.push_back(nS_ + p_);
            }
            else
            {
                U_.push_back(i - p_ + p_);
            }
        }
    }
    else
        throw std::runtime_error("The chosen curve type is not allowed");
}
void SB2::compute_N_()
{
    N_ = std::vector<std::vector<std::vector<SymEngine::Expression>>>(t_ + 2, std::vector<std::vector<SymEngine::Expression>>(nS_, std::vector<SymEngine::Expression>(n_ + 2, 0)));
    std::vector<SymEngine::Expression> buffer; // Save the previous recurrence step
    for (int s = 0; s < nS_; s++)
    {
        N_[0][s][s + p_] = 1;
        for (int d = 1; d < p_ + 1; d++)
        {
            buffer = N_[0][s];
            for (int i = s + p_ - d; i < s + p_ + 1; i++)
            {
                if ((U_[i + d] - U_[i]) == 0 && (U_[i + d + 1] - U_[i + 1]) == 0)
                {
                    N_[0][s][i] = 0;
                }
                else if ((U_[i + d] - U_[i]) == 0)
                {
                    N_[0][s][i] = SymEngine::expand(((U_[i + d + 1] - u_) / (U_[i + d + 1] - U_[i + 1])) * buffer[i + 1]);
                }
                else if ((U_[i + d + 1] - U_[i + 1]) == 0)
                {
                    N_[0][s][i] = SymEngine::expand(((u_ - U_[i]) / (U_[i + d] - U_[i])) * buffer[i]);
                }
                else
                {
                    N_[0][s][i] = SymEngine::expand(((u_ - U_[i]) / (U_[i + d] - U_[i])) * buffer[i] + ((U_[i + d + 1] - u_) / (U_[i + d + 1] - U_[i + 1])) * buffer[i + 1]);
                }
                if (d == p_)
                {
                    compute_horner(N_[0][s][i]);
                }
            }
        }
    }
    if (ct_ == CurveType::UNIFORM_RATIONAL || ct_ == CurveType::CLAMPED_RATIONAL)
    {
        buffer = std::vector<SymEngine::Expression>(nS_, 0);
        for (int s = 0; s < nS_; s++)
        {
            for (int i = s; i < s + p_ + 1; i++)
            {
                buffer[s] += expand(N_[0][s][i] * W_[i]);
                compute_horner(buffer[s]);
            }
            for (int i = s; i < s + p_ + 1; i++)
            {
                N_[0][s][i] = N_[0][s][i] * W_[i] / buffer[s];
            }
        }
    }
    if (f_ == Form::TAYLOR)
    {
        for (int t = 1; t <= t_ + 1; t++)
        {
            for (int s = 0; s < nS_; s++)
            {
                for (int i = s; i < s + p_ + 1; i++)
                {
                    N_[t][s][i] = (N_[t - 1][s][i].diff(u_));
                    if (ct_ == CurveType::UNIFORM_NONRATIONAL || ct_ == CurveType::CLAMPED_NONRATIONAL)
                        compute_horner(N_[t][s][i]);
                }
            }
        }
    }
}
void SB2::compute_horner(SymEngine::Expression &expr)
{
    std::vector<SymEngine::Expression> Coefficients;
    Coefficients.push_back(expr.subs({{u_, SymEngine::Expression(0)}}));
    expr = expr.diff(u_);
    for (int i = 1; i < p_ + 1; i++)
    {
        Coefficients.push_back(expr.subs({{u_, SymEngine::Expression(0)}}) / SymEngine::factorial(i));
        expr = expr.diff(u_);
    }
    expr = (Coefficients[Coefficients.size() - 1]);
    for (int i = 2; i < p_ + 2; i++)
    {
        expr = Coefficients[Coefficients.size() - i] + u_ * expr;
    }
}
void SB2::compute_iN_()
{
    ibex::Variable u("u");
    iN_ = {{}};
    for (int s = 0; s < nS_; s++)
    {
        iN_[0].push_back({});
        for (int i = 0; i < n_ + 1; i++)
        {
            iN_[0][s].push_back(std::make_shared<ibex::Function>(
                "u", ([](std::string expr)
                      { for(size_t pos=0; (pos=expr.find("**",pos))!=std::string::npos; expr.replace(pos,2,"^"), ++pos); return expr; }(N_[0][s][i].get_basic()->__str__()))
                         .c_str()));
        }
    }
    if (f_ == Form::TAYLOR)
    {
        for (int t = 1; t <= t_ + 1; t++)
        {
            iN_.push_back({});
            for (int s = 0; s < nS_; s++)
            {
                iN_[t].push_back({});
                for (int i = 0; i < n_ + 1; i++)
                {
                    iN_[t][s].push_back(std::make_shared<ibex::Function>(
                        "u", ([](std::string expr)
                              { for(size_t pos=0; (pos=expr.find("**",pos))!=std::string::npos; expr.replace(pos,2,"^"), ++pos); return expr; }(N_[t][s][i].get_basic()->__str__()))
                                 .c_str()));
                }
            }
        }
    }
}
void SB2::compute_rdu_()
{
    rdu_ = {};
    for (int s = 0; s < nS_; s++)
    {
        rdu_.push_back({});
        for (int du = 0; du < d_; du++)
        {
            rdu_[s].push_back(ibex::Interval(std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du)), -INFINITY),
                                             std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du)), INFINITY)));
        }
    }
    rdu_[nS_-1].push_back(ibex::Interval(std::nextafter(SymEngine::eval_double(U_[n_] + SymEngine::Expression(1)), -INFINITY),
                                        std::nextafter(SymEngine::eval_double(U_[n_] + SymEngine::Expression(1)), INFINITY))); // only for real-based B-spline: a last point must be evaluated
}
void SB2::compute_idu_()
{
    idu_ = {};
    for (int s = 0; s < nS_; s++)
    {
        idu_.push_back({});
        for (int du = 0; du < d_; du++)
        {
            idu_[s].push_back(ibex::Interval(std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du)), -INFINITY),
                                             std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du + 1)), INFINITY)));
        }
    }
}
void SB2::compute_adu_()
{
    adu_ = {};
    for (int s = 0; s < nS_; s++)
    {
        adu_.push_back({});
        for (int du = 0; du < d_; du++)
        {
            adu_[s].push_back(ibex::Affine2(ibex::Interval(std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du)), -INFINITY),
                                                           std::nextafter(SymEngine::eval_double(U_[p_ + s] + SymEngine::Expression(1) / SymEngine::Expression(d_) * SymEngine::Expression(du + 1)), INFINITY))));
        }
    }
}
int SB2::factorial(int n)
{
    if (n > 12)
    {
        throw std::runtime_error("factorial(" + std::to_string(n) + ") is numerically too big");
    }
    int r = 1;
    for (int i = 2; i <= n; i++)
    {
        r *= i;
    }
    return r;
}
void SB2::compute_reBf()
{
    reBf_ = std::vector<std::vector<std::vector<double>>>(nS_, std::vector<std::vector<double>>(n_ + 1, std::vector<double>(d_, 0.0)));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            for (int du = 0; du < d_; du++)
            {
                if (f_ == Form::TAYLOR)
                {
                    ibex::Interval buffer(0.0);
                    for (int t = 0; t <= t_; t++)
                    {
                        buffer += (iN_[t][s][i]->eval(ibex::IntervalVector(1, rdu_[s][du].mid())) / factorial(t) * pow(rdu_[s][du] - rdu_[s][du].mid(), t)).mid();
                    }
                    buffer += (iN_[t_ + 1][s][i]->eval(ibex::IntervalVector(1, rdu_[s][du])) / factorial(t_ + 1) * pow(rdu_[s][du] - rdu_[s][du].mid(), t_ + 1)).mid();
                    reBf_[s][i][du] = buffer.mid();
                }
                else
                {
                    reBf_[s][i][du] = (iN_[0][s][i]->eval(ibex::IntervalVector(1, rdu_[s][du]))).mid();
                }
            }
        }
    }
    if (f_ == Form::TAYLOR) // only for real-based B-spline: a last point must be evaluated
    {
        for (int i = nS_-1; i<n_+1; i++)
        {
            ibex::Interval buffer(0.0);
            for (int t = 0; t <= t_; t++)
            {
                buffer += (iN_[t][nS_-1][i]->eval(ibex::IntervalVector(1, rdu_[nS_-1][d_].mid())) / factorial(t) * pow(rdu_[nS_-1][d_] - rdu_[nS_-1][d_].mid(), t)).mid();
            }
            buffer += (iN_[t_ + 1][nS_-1][i]->eval(ibex::IntervalVector(1, rdu_[nS_-1][d_])) / factorial(t_ + 1) * pow(rdu_[nS_-1][d_] - rdu_[nS_-1][d_].mid(), t_ + 1)).mid();
            reBf_[nS_-1][i].push_back(buffer.mid());
        }
    }
    else
    {
        for (int i = nS_-1; i<n_+1; i++)
        {
            reBf_[nS_-1][i].push_back((iN_[0][nS_-1][i]->eval(ibex::IntervalVector(1, rdu_[nS_-1][d_]))).mid());
        }
    }
}
void SB2::compute_ieBf()
{
    ieBf_ = std::vector<std::vector<std::vector<ibex::Interval>>>(nS_, std::vector<std::vector<ibex::Interval>>(n_ + 1, std::vector<ibex::Interval>(d_, ibex::Interval(0.0))));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            for (int du = 0; du < d_; du++)
            {
                if (f_ == Form::TAYLOR)
                {
                    ibex::Interval buffer(0, 0);
                    for (int t = 0; t <= t_; t++)
                    {
                        buffer += iN_[t][s][i]->eval(ibex::IntervalVector(1, idu_[s][du].mid())) / factorial(t) * pow(idu_[s][du] - idu_[s][du].mid(), t);
                    }
                    buffer += iN_[t_ + 1][s][i]->eval(ibex::IntervalVector(1, idu_[s][du])) / factorial(t_ + 1) * pow(idu_[s][du] - idu_[s][du].mid(), t_ + 1);
                    ieBf_[s][i][du] = buffer & ibex::Interval(0, 1);
                }
                else
                {
                    ieBf_[s][i][du] = iN_[0][s][i]->eval(ibex::IntervalVector(1, idu_[s][du])) & ibex::Interval(0, 1);
                }
            }
        }
    }
}
void SB2::compute_aeBf()
{
    aeBf_ = std::vector<std::vector<std::vector<ibex::Affine2>>>(nS_, std::vector<std::vector<ibex::Affine2>>(n_ + 1, std::vector<ibex::Affine2>(d_, ibex::Affine2(0.0))));
    for (int s = 0; s < nS_; s++)
    {
        for (int i = s; i < s + p_ + 1; i++)
        {
            for (int du = 0; du < d_; du++)
            {
                if (f_ == Form::TAYLOR)
                {
                    ibex::Affine2 buffer(0);
                    for (int t = 0; t <= t_; t++)
                    {
                        buffer += iN_[t][s][i]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du].mid())) / factorial(t) * pow(adu_[s][du] - adu_[s][du].mid(), t);
                    }
                    buffer += iN_[t_ + 1][s][i]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du])) / factorial(t_ + 1) * pow(adu_[s][du] - adu_[s][du].mid(), t_ + 1);
                    aeBf_[s][i][du] = buffer;
                }
                else
                {
                    aeBf_[s][i][du] = iN_[0][s][i]->eval_affine2(ibex::Affine2Vector(1, adu_[s][du]));
                }
            }
        }
    }
}
namespace
{
template <typename T> // std::vector<double>, ibex::IntervalVector or ibex::Affine2Vector
void check_control_points(const std::vector<T> &P, const int n)
{
    for (size_t dim = 0; dim < P.size(); dim++)
    {
        if ((int)P[dim].size() < n + 1)
            throw std::runtime_error("A control point dimension has fewer than n+1 control points");
    }
}
}
int SB2::du_count(const int s) const
{
    if (ps_ == ParameterSet::R && s == nS_ - 1) // only for real-based B-spline: a last point must be evaluated
        return d_ + 1;
    return d_;
}
std::pair<int, int> SB2::impacted_segments(const int i) const
{
    if (i < 0 || i > n_)
        return std::pair<int, int>(0, -1); // empty range: i is not a control point index
    return std::pair<int, int>(std::max(0, i - p_), std::min(nS_ - 1, i));
}
void SB2::eval_point_segment(const std::vector<std::vector<double>> &P, const int s, std::vector<std::vector<double>> &points)
{
    const int nb = du_count(s);
    points.assign(nb, std::vector<double>(P.size(), 0.0));
    for (size_t dim = 0; dim < P.size(); dim++)
    {
        for (int du = 0; du < nb; du++)
        {
            double buffer(0.0);
            switch (ps_)
            {
            case ParameterSet::R:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * reBf_[s][i][du];
                }
                break;
            case ParameterSet::IR:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * ieBf_[s][i][du].mid();
                }
                break;
            case ParameterSet::Z:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * aeBf_[s][i][du].mid();
                }
                break;
            default:
                throw std::runtime_error("Unknown ParameterSet");
            }
            points[du][dim] = buffer;
        }
    }
}
void SB2::eval_box_segment(const std::vector<ibex::IntervalVector> &P, const int s, std::vector<ibex::IntervalVector> &boxes)
{
    const int nb = du_count(s);
    boxes.assign(nb, ibex::IntervalVector(P.size()));
    for (size_t dim = 0; dim < P.size(); dim++)
    {
        for (int du = 0; du < nb; du++)
        {
            switch (ps_)
            {
            case ParameterSet::R:
            {
                ibex::Interval buffer(0.0);
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * reBf_[s][i][du];
                }
                boxes[du][dim] = buffer;
                break;
            }
            case ParameterSet::IR:
            {
                ibex::Interval buffer(0.0);
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * ieBf_[s][i][du];
                }
                boxes[du][dim] = buffer;
                break;
            }
            case ParameterSet::Z:
            {
                if (LWE) // here the basis and the final B-spline are computed using affine arithmetic and then the hull is extracted (reduce the wrapping effect)
                {
                    ibex::Affine2 buffer(0.0);
                    for (int i = s; i < s + p_ + 1; i++)
                    {
                        buffer += P[dim][i] * aeBf_[s][i][du];
                    }
                    boxes[du][dim] = buffer.itv(); // hull extraction
                }
                else // here the basis is computed using affine arithmetic and then the hull is extracted to compute the final B-spline
                {
                    ibex::Interval buffer(0.0);
                    for (int i = s; i < s + p_ + 1; i++)
                    {
                        buffer += P[dim][i] * aeBf_[s][i][du].itv(); // hull extraction
                    }
                    boxes[du][dim] = buffer;
                }
                break;
            }
            default:
                throw std::runtime_error("Unknown ParameterSet");
            }
        }
    }
}
void SB2::eval_zonotope_segment(const std::vector<ibex::Affine2Vector> &P, const int s, std::vector<ibex::Affine2Vector> &zonotopes)
{
    const int nb = du_count(s);
    zonotopes.assign(nb, ibex::Affine2Vector(P.size()));
    for (size_t dim = 0; dim < P.size(); dim++)
    {
        for (int du = 0; du < nb; du++)
        {
            ibex::Affine2 buffer(0.0);
            switch (ps_)
            {
            case ParameterSet::R:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * reBf_[s][i][du];
                }
                break;
            case ParameterSet::IR:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * ieBf_[s][i][du];
                }
                break;
            case ParameterSet::Z:
                for (int i = s; i < s + p_ + 1; i++)
                {
                    buffer += P[dim][i] * aeBf_[s][i][du];
                }
                break;
            default:
                throw std::runtime_error("Unknown ParameterSet");
            }
            zonotopes[du][dim] = buffer;
        }
    }
}
std::vector<std::vector<std::vector<double>>> SB2::eval_point(const std::vector<std::vector<double>> &P)
{
    check_control_points(P, n_);
    std::vector<std::vector<std::vector<double>>> points(nS_);
    for (int s = 0; s < nS_; s++)
    {
        eval_point_segment(P, s, points[s]);
    }
    return points;
}
std::vector<std::vector<ibex::IntervalVector>> SB2::eval_box(const std::vector<ibex::IntervalVector> &P)
{
    check_control_points(P, n_);
    std::vector<std::vector<ibex::IntervalVector>> boxes(nS_);
    for (int s = 0; s < nS_; s++)
    {
        eval_box_segment(P, s, boxes[s]);
    }
    return boxes;
}
std::vector<std::vector<ibex::Affine2Vector>> SB2::eval_zonotope(const std::vector<ibex::Affine2Vector> &P)
{
    check_control_points(P, n_);
    std::vector<std::vector<ibex::Affine2Vector>> zonotopes(nS_);
    for (int s = 0; s < nS_; s++)
    {
        eval_zonotope_segment(P, s, zonotopes[s]);
    }
    return zonotopes;
}
void SB2::update_point(const std::vector<std::vector<double>> &P, const int i, std::vector<std::vector<std::vector<double>>> &points)
{
    check_control_points(P, n_);
    if ((int)points.size() != nS_)
        throw std::runtime_error("The elements to update do not come from this B-spline");
    const std::pair<int, int> seg = impacted_segments(i);
    for (int s = seg.first; s <= seg.second; s++)
    {
        eval_point_segment(P, s, points[s]);
    }
}
void SB2::update_box(const std::vector<ibex::IntervalVector> &P, const int i, std::vector<std::vector<ibex::IntervalVector>> &boxes)
{
    check_control_points(P, n_);
    if ((int)boxes.size() != nS_)
        throw std::runtime_error("The elements to update do not come from this B-spline");
    const std::pair<int, int> seg = impacted_segments(i);
    for (int s = seg.first; s <= seg.second; s++)
    {
        eval_box_segment(P, s, boxes[s]);
    }
}
void SB2::update_zonotope(const std::vector<ibex::Affine2Vector> &P, const int i, std::vector<std::vector<ibex::Affine2Vector>> &zonotopes)
{
    check_control_points(P, n_);
    if ((int)zonotopes.size() != nS_)
        throw std::runtime_error("The elements to update do not come from this B-spline");
    const std::pair<int, int> seg = impacted_segments(i);
    for (int s = seg.first; s <= seg.second; s++)
    {
        eval_zonotope_segment(P, s, zonotopes[s]);
    }
}

}
