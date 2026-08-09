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
#include <cmath>
#include <list>
#include <stdexcept>
#include <string>
#include <utility>

// When the parameter is affine and the control points are boxes, evaluate the
// whole linear combination in affine arithmetic and take the hull of the
// result, instead of taking the hull of each basis coefficient first. The
// dependencies between the coefficients are then kept and the enclosure is
// narrower, at the price of the affine operations. Defined here so that a
// build may override it with -DSB2L_AFFINE_BEFORE_HULL=0.
#ifndef SB2L_AFFINE_BEFORE_HULL
#define SB2L_AFFINE_BEFORE_HULL 1
#endif

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
    if (p_ < 0)
        throw std::runtime_error("The degree of the curve p must be positive");
    if (p_ >= n_ + 1)
        throw std::runtime_error("The degree of the curve p must be lower than the number of control points n+1");
    if (d_ < 1)
        throw std::runtime_error("The number of subdivisions d must be at least 1");
    if (t_ < -1)
        throw std::runtime_error("The Taylor order t must be positive, or -1 to select it automatically");
    if (f_ == Form::TAYLOR)
    {
        if (t_ == -1) // by default, the taylor order is chosen to minimize the overestimation
        {
            // One order below the degree, above which the remainder vanishes.
            // 11 is the numerical limit: the remainder divides by (t+1)! and
            // the factorial is computed on an int.
            t_ = std::min(std::max(p_ - 1, 0), 11);
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
        else if ((int)W_.size() != n_ + 1)
            throw std::runtime_error("The rational weight vector must hold one weight per control point");
    }
    else if (W_.size() != 0) // weights only mean something for a rational basis
        throw std::runtime_error("Rational weights were given for a non-rational curve type");
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
int SB2::get_p() const { return p_; }
int SB2::get_n() const { return n_; }
int SB2::get_nS() const { return nS_; }
int SB2::get_d() const { return d_; }
void SB2::set_affine_noise_number(const unsigned int n) { ibex::AF_fAFFullI::setAffineNoiseNumber(n); }
unsigned int SB2::get_affine_noise_number() { return ibex::AF_fAFFullI::getAffineNoiseNumber(); }
std::vector<std::vector<std::vector<double>>> SB2::get_reBf() const
{
    if (ps_ == ParameterSet::R)
        return reBf_;
    // Converted on a local copy: a getter never rewrites the basis the
    // evaluation reads.
    std::vector<std::vector<std::vector<double>>> out(nS_, std::vector<std::vector<double>>(n_ + 1, std::vector<double>(d_, 0.0)));
    for (int s = 0; s < nS_; s++)
        for (int i = s; i < s + p_ + 1; i++)
            for (int du = 0; du < d_; du++)
                out[s][i][du] = (ps_ == ParameterSet::IR) ? ieBf_[s][i][du].mid() : aeBf_[s][i][du].mid();
    return out;
}
std::vector<std::vector<std::vector<ibex::Interval>>> SB2::get_ieBf() const
{
    if (ps_ == ParameterSet::R)
        throw std::runtime_error("Can not build an interval basis from a real one");
    if (ps_ == ParameterSet::IR)
        return ieBf_;
    std::vector<std::vector<std::vector<ibex::Interval>>> out(nS_, std::vector<std::vector<ibex::Interval>>(n_ + 1, std::vector<ibex::Interval>(d_, ibex::Interval(0.0))));
    for (int s = 0; s < nS_; s++)
        for (int i = s; i < s + p_ + 1; i++)
            for (int du = 0; du < d_; du++)
                out[s][i][du] = aeBf_[s][i][du].itv();
    return out;
}
std::vector<std::vector<std::vector<ibex::Affine2>>> SB2::get_aeBf() const
{
    if (ps_ == ParameterSet::R)
        throw std::runtime_error("Can not build an affine basis from a real one");
    if (ps_ == ParameterSet::IR)
        throw std::runtime_error("Can not build an affine basis from an interval one");
    return aeBf_;
}
std::vector<std::vector<std::vector<double>>> SB2::get_reBf_diff(const int order) const
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
std::vector<std::vector<std::vector<ibex::Interval>>> SB2::get_ieBf_diff(const int order) const
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
std::vector<std::vector<std::vector<ibex::Affine2>>> SB2::get_aeBf_diff(const int order) const
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
std::vector<std::vector<ibex::Interval>> SB2::get_rdu() const
{
    if (ps_ == ParameterSet::R)
        return rdu_;
    std::vector<std::vector<ibex::Interval>> out(nS_);
    for (int s = 0; s < nS_; s++)
        for (int du = 0; du < d_; du++)
            out[s].push_back(ibex::Interval((ps_ == ParameterSet::IR) ? idu_[s][du].mid() : adu_[s][du].mid()));
    return out;
}
std::vector<std::vector<ibex::Interval>> SB2::get_idu() const
{
    if (ps_ == ParameterSet::R)
        throw std::runtime_error("Can not build an interval decomposition from a real one");
    if (ps_ == ParameterSet::IR)
        return idu_;
    std::vector<std::vector<ibex::Interval>> out(nS_);
    for (int s = 0; s < nS_; s++)
        for (int du = 0; du < d_; du++)
            out[s].push_back(adu_[s][du].itv());
    return out;
}
std::vector<std::vector<ibex::Affine2>> SB2::get_adu() const
{
    if (ps_ == ParameterSet::R)
        throw std::runtime_error("Can not build an affine decomposition from a real one");
    if (ps_ == ParameterSet::IR)
        throw std::runtime_error("Can not build an affine decomposition from an interval one");
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
                U_.push_back(i);
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
                    compute_horner(N_[0][s][i], p_); // a basis function has the degree of the curve
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
            }
            // Once the whole denominator is summed. It is a combination of the
            // basis functions, so it has their degree.
            compute_horner(buffer[s], p_);
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
                    // Differentiating a polynomial lowers its degree, so the
                    // Horner form only needs that many coefficients.
                    if (ct_ == CurveType::UNIFORM_NONRATIONAL || ct_ == CurveType::CLAMPED_NONRATIONAL)
                        compute_horner(N_[t][s][i], p_ - t);
                }
            }
        }
    }
}
void SB2::compute_horner(SymEngine::Expression &expr, const int degree)
{
    const int deg = std::max(degree, 0);
    std::vector<SymEngine::Expression> Coefficients;
    Coefficients.push_back(expr.subs({{u_, SymEngine::Expression(0)}}));
    expr = expr.diff(u_);
    for (int i = 1; i < deg + 1; i++)
    {
        Coefficients.push_back(expr.subs({{u_, SymEngine::Expression(0)}}) / SymEngine::factorial(i));
        expr = expr.diff(u_);
    }
    expr = (Coefficients[Coefficients.size() - 1]);
    for (int i = 2; i < deg + 2; i++)
    {
        expr = Coefficients[Coefficients.size() - i] + u_ * expr;
    }
}
void SB2::compute_iN_()
{
    // Only the p+1 basis functions of the support of a segment are non-zero;
    // the rows are kept at their full width so that a function is still read
    // at its control point index, but the entries outside the support stay
    // empty instead of parsing an expression which is zero. On a curve with
    // many control points, that is the bulk of the construction.
    auto to_ibex = [](std::string expr) { for(size_t pos=0; (pos=expr.find("**",pos))!=std::string::npos; expr.replace(pos,2,"^"), ++pos); return expr; };
    const int nT = (f_ == Form::TAYLOR) ? t_ + 1 : 0;
    iN_.assign(nT + 1, std::vector<std::vector<std::shared_ptr<ibex::Function>>>(
                           nS_, std::vector<std::shared_ptr<ibex::Function>>(n_ + 1)));
    for (int t = 0; t <= nT; t++)
    {
        for (int s = 0; s < nS_; s++)
        {
            for (int i = s; i < s + p_ + 1; i++)
            {
                iN_[t][s][i] = std::make_shared<ibex::Function>(
                    "u", to_ibex(N_[t][s][i].get_basic()->__str__()).c_str());
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
    if (P.size() == 0)
        throw std::runtime_error("The control points carry no coordinate: P must hold one container per dimension");
    for (size_t dim = 0; dim < P.size(); dim++)
    {
        if ((int)P[dim].size() < n + 1)
            throw std::runtime_error("A control point dimension has fewer than n+1 control points");
        if (P[dim].size() != P[0].size())
            throw std::runtime_error("The control point dimensions do not hold the same number of control points");
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
                if (SB2L_AFFINE_BEFORE_HULL) // the basis and the B-spline are combined in affine arithmetic, then the hull is extracted
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

Zonotope zonotope_of(const ibex::Affine2Vector &v)
{
    Zonotope z;
    z.dim = v.size();
    if (z.dim < 1)
        return z;
    z.center.resize(z.dim);
    for (int i = 0; i < z.dim; i++)
        z.center[i] = v[i].val(0);

    // Every coordinate draws on the same noise terms but keeps only its own
    // non-zero ones, and each list is sorted by term number. Walking the lists
    // together, smallest number first, gives one generator per term any
    // coordinate uses, with the coordinates which do not use it left at zero.
    typedef std::list<std::pair<int, double>>::const_iterator It;
    std::vector<const std::list<std::pair<int, double>> *> rays(z.dim);
    std::vector<It> it(z.dim);
    for (int i = 0; i < z.dim; i++)
    {
        rays[i] = &v[i].rays();
        it[i] = rays[i]->begin();
    }
    for (;;)
    {
        int number = -1; // smallest term number still ahead, over all coordinates
        for (int i = 0; i < z.dim; i++)
            if (it[i] != rays[i]->end() && (number < 0 || it[i]->first < number))
                number = it[i]->first;
        if (number < 0)
            break;
        const size_t base = z.generators.size();
        z.generators.resize(base + z.dim, 0.0);
        for (int i = 0; i < z.dim; i++)
            if (it[i] != rays[i]->end() && it[i]->first == number)
                z.generators[base + i] = (it[i]++)->second;
        z.m++;
    }
    // The error term each coordinate accumulated, along its own axis. It is
    // part of the element, so leaving it out would describe less than it.
    for (int i = 0; i < z.dim; i++)
    {
        if (v[i].err() <= 0.0)
            continue;
        const size_t base = z.generators.size();
        z.generators.resize(base + z.dim, 0.0);
        z.generators[base + i] = v[i].err();
        z.m++;
    }
    return z;
}

}
