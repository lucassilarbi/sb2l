/**
 * \date 03/2026
 * \author CEA/DRT/LIST/DIASI/SRI/LCSR
 * \author L.Si Larbi
 *
 * @par Licence
 * Copyright © 2026 CEA
 */

#ifndef INTERVAL_BSPLINE_HPP_
#define INTERVAL_BSPLINE_HPP_

#include <symengine/expression.h>
#include <ibex/ibex.h>
#include <memory>

namespace sb2l {

/**
 * \brief The type of the curve direcly link to the knot vector
 * UNIFORM : Uniform knot vector
 * CLAMPED : The first and lasts knots have a multiplicity of p+1, the resulting B-spline is clamped
 * RATIONAL : The basis is rational. Each control point is link to a dedicated weight
 */
enum class CurveType{UNIFORM_RATIONAL, UNIFORM_NONRATIONAL, CLAMPED_RATIONAL, CLAMPED_NONRATIONAL};
inline std::ostream& operator<<(std::ostream& os, CurveType ct) 
{
    switch (ct) 
    {
        case CurveType::UNIFORM_RATIONAL: return os << "UNIFORM_RATIONAL";
        case CurveType::UNIFORM_NONRATIONAL: return os << "UNIFORM_NONRATIONAL";
        case CurveType::CLAMPED_RATIONAL: return os << "CLAMPED_RATIONAL";
        case CurveType::CLAMPED_NONRATIONAL: return os << "CLAMPED_NONRATIONAL";
        default: return os << "UNKNOWN";
    }
}
/**
 * \brief The form of the equation used to evaluate the B-spline
 * NATURAL : Natural extension of the function using extended usual operators from R
 * TAYLOR : Taylor form using extended usual operators from R, the talyor order need to be set
 */
enum class Form{NATURAL, TAYLOR};
inline std::ostream& operator<<(std::ostream& os, Form f) 
{
    switch (f) 
    {
        case Form::NATURAL: return os << "NATURAL";
        case Form::TAYLOR: return os << "TAYLOR";
        default: return os << "UNKNOWN";
    }
}
/**
 * \brief The parameter set
 * R : Reals
 * IR : Intervals
 * Z : Zonotopes
 */
enum class ParameterSet{R, IR, Z};
inline std::ostream& operator<<(std::ostream& os, ParameterSet ps) 
{
    switch (ps) 
    {
        case ParameterSet::R: return os << "R";
        case ParameterSet::IR: return os << "IR";
        case ParameterSet::Z: return os << "Z";
        default: return os << "UNKNOWN";
    }
}

class SB2
{
    public:

    /**
     * @brief Constructor
     */
    SB2(const int p = 3,
         const int nCP = 5,
         const CurveType ct = CurveType::UNIFORM_NONRATIONAL,
         const Form f = Form::TAYLOR,
         const ParameterSet ps = ParameterSet::IR,
         const int d = 10, 
         const int t=-1,
         const std::vector<SymEngine::Expression> W=std::vector<SymEngine::Expression>({})
        );
    /**
     * @brief get p_
     */
    int get_p();
    /**
     * @brief get n_
     */
    int get_n();
    /**
     * @brief get nS_
     */
    int get_nS();
    /**
     * @brief get d_
     */
    int get_d();
    /**
     * @brief get ieBf_
     */
    std::vector<std::vector<std::vector<double>>> get_reBf();
    std::vector<std::vector<std::vector<ibex::Interval>>> get_ieBf();
    std::vector<std::vector<std::vector<ibex::Affine2>>> get_aeBf();
    /**
     * @brief get idu_
     */
    std::vector<std::vector<ibex::Interval>> get_rdu();
    std::vector<std::vector<ibex::Interval>> get_idu();
    std::vector<std::vector<ibex::Affine2>> get_adu();
    /**
     * @brief evaluate the B-spline. Takes a list of control points as input and return the list corresponding elements.
     * The resulted elements are arranged by segment
     */
    std::vector<std::vector<std::vector<double>>> eval_point(const std::vector<std::vector<double>> &P);
    std::vector<std::vector<ibex::IntervalVector>> eval_box(const std::vector<ibex::IntervalVector> &P);
    std::vector<std::vector<ibex::Affine2Vector>> eval_zonotope(const std::vector<ibex::Affine2Vector> &P);
    
    private:

    int p_; // Curve degree
    int n_; // number of control boxes minus one
    CurveType ct_; // Curve type
    Form f_; // basis equation desired form
    ParameterSet ps_; // chosen set for the B-spline parameter
    int d_; // decomposition: number of parameter intervals between two knots
    int t_; // taylor order (0 correspond to the centered form)
    std::vector<SymEngine::Expression> W_; // rational weights used to compute the rational basis
    int k_; // k+1 knots in the knot vector
    int nS_; // number of segments
    std::vector<SymEngine::Expression> U_; // Knot vector
    const SymEngine::Expression u_ = SymEngine::Expression("u"); // Symbolic B-spline parameter
    std::vector<std::vector<std::vector<SymEngine::Expression>>> N_; // B-spline basis functions arranged by segment
    std::vector<std::vector<std::vector<std::shared_ptr<ibex::Function>>>> iN_; // interval B-spline basis functions arranged by segment
    std::vector<std::vector<ibex::Interval>> rdu_; // real decomposition of the parameter arranged by segment (used interval for numerical guarantee)
    std::vector<std::vector<ibex::Interval>> idu_; // interval decomposition of the parameter arranged by segment
    std::vector<std::vector<ibex::Affine2>> adu_; // affine decomposition of the parameter arranged by segment
    std::vector<std::vector<std::vector<double>>> reBf_; // real evaluation of Basis functions arranged by segment (used interval for numerical guarantee)
    std::vector<std::vector<std::vector<ibex::Interval>>> ieBf_; // interval evaluation of Basis functions arranged by segment
    std::vector<std::vector<std::vector<ibex::Affine2>>> aeBf_; // affine evaluation of Basis functions arranged by segment

    /**
     * @brief compute U_
     */
    void compute_U_();
    /**
     * @brief compute N_
     */
    void compute_N_();
    /**
     * @brief compute the hroner form of a given SymEngine Expression
     */
    void compute_horner(SymEngine::Expression& expr);
    /**
     * @brief compute iN_
     */
    void compute_iN_();
    /**
     * @brief compute Xdu_
     */
    void compute_rdu_();
    void compute_idu_();
    void compute_adu_();
    /**
     * @brief compute ieBf_
     */
    int factorial(int n); // numerically limited to n <= 12
    void compute_reBf();
    void compute_ieBf();
    void compute_aeBf();
};

}

#endif // INTERVAL_BSPLINE_HPP_