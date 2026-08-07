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
#include <utility>
#include <vector>

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
     * @param p degree of the curve, p >= 1 (p >= 0 is accepted but gives a constant basis)
     * @param nCP number of control points, nCP > p
     * @param d number of parameter subdivisions per segment, d >= 1
     * @param t Taylor order, only used when f is Form::TAYLOR. -1 selects it
     *          automatically (p-1, capped at 11 by the factorial). Any other
     *          negative value is rejected; values above 11 are capped.
     * Every argument is validated: an out-of-range one raises std::runtime_error
     * rather than leaving the object in an undefined state.
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
     * @brief number of noise symbols an affine form keeps before the smallest
     * ones are merged into its error term. This is process-wide state of the
     * affine arithmetic, not a property of one B-spline; the default of the
     * underlying library is 10, which silently reduces the order of any
     * zonotope built from more than 10 generators. Raise it before building
     * control zonotopes if their exact shape matters.
     */
    static void set_affine_noise_number(const unsigned int n);
    static unsigned int get_affine_noise_number();
    /**
     * @brief get p_
     */
    int get_p() const;
    /**
     * @brief get n_
     */
    int get_n() const;
    /**
     * @brief get nS_
     */
    int get_nS() const;
    /**
     * @brief get d_
     */
    int get_d() const;
    /**
     * @brief get basis. The basis is returned in the arithmetic asked for,
     * converted from the one it was built in when they differ (an interval
     * basis is the hull of an affine one, a real basis its midpoint); the
     * conversions that lose the enclosure, from a poorer to a richer
     * arithmetic, raise std::runtime_error. These getters do not modify the
     * B-spline.
     */
    std::vector<std::vector<std::vector<double>>> get_reBf() const;
    std::vector<std::vector<std::vector<ibex::Interval>>> get_ieBf() const;
    std::vector<std::vector<std::vector<ibex::Affine2>>> get_aeBf() const;
    /**
     * @brief evaluate the order-th derivative of the basis functions, using the same form as the basis
     */
    std::vector<std::vector<std::vector<double>>> get_reBf_diff(const int order) const;
    std::vector<std::vector<std::vector<ibex::Interval>>> get_ieBf_diff(const int order) const;
    std::vector<std::vector<std::vector<ibex::Affine2>>> get_aeBf_diff(const int order) const;
    /**
     * @brief get the subdivisions of the parameter, in the arithmetic asked
     * for, under the same conversion rules as the basis getters above.
     */
    std::vector<std::vector<ibex::Interval>> get_rdu() const;
    std::vector<std::vector<ibex::Interval>> get_idu() const;
    std::vector<std::vector<ibex::Affine2>> get_adu() const;
    /**
     * @brief evaluate the B-spline. Takes a list of control points as input and return the list corresponding elements.
     * The resulted elements are arranged by segment: out[s][du] is the element
     * of segment s on subdivision du. P holds one container per coordinate,
     * all of the same size, at least n+1 entries each.
     *
     * eval_box and eval_zonotope enclose the curve. eval_point does not: it
     * returns one point per subdivision, obtained from the midpoint of the
     * basis, so it only coincides with the curve when the parameter set is
     * ParameterSet::R. Over an interval or affine parameter it is a
     * representative of the element, with no enclosure property; use
     * eval_box or eval_zonotope when the guarantee is needed.
     */
    std::vector<std::vector<std::vector<double>>> eval_point(const std::vector<std::vector<double>> &P);
    std::vector<std::vector<ibex::IntervalVector>> eval_box(const std::vector<ibex::IntervalVector> &P);
    std::vector<std::vector<ibex::Affine2Vector>> eval_zonotope(const std::vector<ibex::Affine2Vector> &P);
    /**
     * @brief first and last segment (both included) whose elements depend on the control point i.
     * A segment s is only built from the control points i in [s, s+p], hence the control point i
     * only appears in the segments s in [i-p, i]: moving it invalidates at most p+1 segments.
     * An empty range (first > last) is returned when i is not a control point index.
     */
    std::pair<int, int> impacted_segments(const int i) const;
    /**
     * @brief re-evaluate, in place, the only segments impacted by the control point i.
     * The elements must come from the matching eval_* call on the same B-spline, with the same
     * control points except the i-th one. Every other segment is left untouched.
     */
    void update_point(const std::vector<std::vector<double>> &P, const int i, std::vector<std::vector<std::vector<double>>> &points);
    void update_box(const std::vector<ibex::IntervalVector> &P, const int i, std::vector<std::vector<ibex::IntervalVector>> &boxes);
    void update_zonotope(const std::vector<ibex::Affine2Vector> &P, const int i, std::vector<std::vector<ibex::Affine2Vector>> &zonotopes);

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
    static int factorial(int n); // numerically limited to n <= 12
    void compute_reBf();
    void compute_ieBf();
    void compute_aeBf();
    /**
     * @brief number of parameter values evaluated on the segment s (the real-based B-spline
     * evaluates one last parameter value at the end of the last segment)
     */
    int du_count(const int s) const;
    /**
     * @brief evaluate the single segment s, used by both the whole and the incremental evaluations
     */
    void eval_point_segment(const std::vector<std::vector<double>> &P, const int s, std::vector<std::vector<double>> &points);
    void eval_box_segment(const std::vector<ibex::IntervalVector> &P, const int s, std::vector<ibex::IntervalVector> &boxes);
    void eval_zonotope_segment(const std::vector<ibex::Affine2Vector> &P, const int s, std::vector<ibex::Affine2Vector> &zonotopes);
};

}

#endif // INTERVAL_BSPLINE_HPP_