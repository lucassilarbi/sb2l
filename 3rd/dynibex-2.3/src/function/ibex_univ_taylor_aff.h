/* ============================================================================
 * D Y N I B E X - Definition of the univariated taylor component computation using affine forms
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : 2015
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique,
 * 		 THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef UNIV_TAYLOR_AFF_H
#define UNIV_TAYLOR_AFF_H

#include "ibex_Function.h"
#include "ibex_FwdAlgorithm.h"
#include "ibex_BwdAlgorithm.h"
#include "ibex_Affine2Domain.h"

//decomment if you want "verbose style" of the computation
// #ifndef __PRINT_DEBUG__
// #define __PRINT_DEBUG__
// #endif


using namespace std;

namespace ibex{

class Univ_Taylor_Aff : public FwdAlgorithm, public BwdAlgorithm
{
private:
    /**
     * \brief Max order of the Taylor serie (equals the max number of taylor components)
     */
    unsigned int _max_order;

    /**
     * \brief variable on which we want to compute the Taylor coefficient.
     */
    const ExprSymbol * _t;

public:

    /**
     * \brief constructor of Univ_Taylor
     * \param max_order : max order of the Taylor components we want to compute
     * \param t : variable on which we compute the Taylor components
     */
    Univ_Taylor_Aff(unsigned int max_order, const ExprSymbol * t){
	_max_order = max_order;
        _t = t;
    }

    virtual ~Univ_Taylor_Aff() {};

    inline const ExprSymbol * get_symbol() {
        return _t;
    }

    inline unsigned int get_max_order() {
        return _max_order;
    }
    inline void set_max_order(unsigned int max_order) {
        _max_order = max_order;
    }

    /**
     * \brief Calculate the Taylor components of f and store the result in \a taylor_comps.
     */
    void compute_taylor_aff(const Function& f, std::vector<Affine2Domain>& taylor_comps);

    void print_taylor_comps(std::vector<Affine2Domain> *);

//------------------------------------------------------------------------------------------------------------
//inherited methods from FwdAlgorithm

    /**
     * \brief y becomes the ind^th argument of the vector x
    */
    inline void index_fwd(const ExprIndex&, const ExprLabel& x, ExprLabel& y);

    void vector_fwd(const ExprVector&, const ExprLabel** compL, ExprLabel& y) {
      cout << __func__ << endl;
    }

    /**
    * \brief in case of a constant, only the existing components are copied.
    */
    inline void cst_fwd(const ExprConstant&, ExprLabel& y);
    /**
    * \brief if the symbol x is the variable for which we want the taylor component, its Taylor components are (x, 1), if not it is a constant and the Taylor component are (x)
    */
    inline void symbol_fwd(const ExprSymbol&, ExprLabel& y);

    inline void apply_fwd(const ExprApply&, ExprLabel** x, ExprLabel& y);
    inline void chi_fwd(const ExprChi&, const ExprLabel& x1, const ExprLabel& x2, const ExprLabel& x3, ExprLabel& y);

    /**
    * \brief Addition for two arrays of taylor coefficients is the sum componentwise.
    */
    inline void add_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);

    /**
    * \brief taylor components of the multiplication of two arguments
    */
    inline void mul_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);

    /**
    * \brief same as add_fwd
    */
    inline void sub_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void div_fwd(const ExprDiv&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void max_fwd(const ExprMax&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void min_fwd(const ExprMin&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void atan2_fwd(const ExprAtan2&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void minus_fwd(const ExprMinus&, const ExprLabel& x, ExprLabel& y);
    inline void trans_V_fwd(const ExprTrans&, const ExprLabel& x, ExprLabel& y);
    inline void trans_M_fwd(const ExprTrans&, const ExprLabel& x, ExprLabel& y);
    inline void sign_fwd(const ExprSign&, const ExprLabel& x, ExprLabel& y);
    inline void abs_fwd(const ExprAbs&, const ExprLabel& x, ExprLabel& y);
    inline void power_fwd(const ExprPower& p, const ExprLabel& x, ExprLabel& y);
    inline void sqr_fwd(const ExprSqr&, const ExprLabel& x, ExprLabel& y);
    inline void sqrt_fwd(const ExprSqrt&, const ExprLabel& x, ExprLabel& y);
    inline void exp_fwd(const ExprExp&, const ExprLabel& x, ExprLabel& y);
    inline void log_fwd(const ExprLog&, const ExprLabel& x, ExprLabel& y);
    inline void cos_fwd(const ExprCos&, const ExprLabel& x, ExprLabel& y);
    inline void sin_fwd(const ExprSin&, const ExprLabel& x, ExprLabel& y);
    inline void tan_fwd(const ExprTan&, const ExprLabel& x, ExprLabel& y);
    inline void cosh_fwd(const ExprCosh&, const ExprLabel& x, ExprLabel& y);
    inline void sinh_fwd(const ExprSinh&, const ExprLabel& x, ExprLabel& y);
    inline void tanh_fwd(const ExprTanh&, const ExprLabel& x, ExprLabel& y);
    inline void acos_fwd(const ExprAcos&, const ExprLabel& x, ExprLabel& y);
    inline void asin_fwd(const ExprAsin&, const ExprLabel& x, ExprLabel& y);
    inline void atan_fwd(const ExprAtan&, const ExprLabel& x, ExprLabel& y);
    inline void acosh_fwd(const ExprAcosh&, const ExprLabel& x, ExprLabel& y);
    inline void asinh_fwd(const ExprAsinh&, const ExprLabel& x, ExprLabel& y);
    inline void atanh_fwd(const ExprAtanh&, const ExprLabel& x, ExprLabel& y);

    /**
    * \brief addition between two vectors (Taylor components are added componentwise)
    */
    inline void add_V_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void add_M_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_SV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_SM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_VV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_MV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_VM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void mul_MM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void sub_V_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);
    inline void sub_M_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y);


    	       void index_bwd (const ExprIndex&,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear();x.taylor_comps_aff->clear(); }
	       void vector_bwd(const ExprVector& v,  ExprLabel** compL,               const ExprLabel& y){
		 y.taylor_comps_aff->clear();
		 for (int i = 0; i<v.length(); i++)
		   compL[i]->taylor_comps_aff->clear();
	      }
	inline void symbol_bwd(const ExprSymbol& x,                               const ExprLabel& y){y.taylor_comps_aff->clear();}
	inline void cst_bwd   (const ExprConstant&,                             const ExprLabel& y){y.taylor_comps_aff->clear(); }
	       void apply_bwd (const ExprApply& a,  ExprLabel** x,                const ExprLabel& y){
		 y.taylor_comps_aff->clear();
		 for (int i = 0; i < a.nb_args; i++)
		   x[i]->taylor_comps_aff->clear();

	      }
	       void chi_bwd   (const ExprChi&,    ExprLabel& x1, ExprLabel& x2, ExprLabel& x3, const ExprLabel& y){y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();x3.taylor_comps_aff->clear(); }
	inline void add_bwd   (const ExprAdd&,    ExprLabel& x1, ExprLabel& x2, const ExprLabel& y) {y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	inline void mul_bwd   (const ExprMul&,    ExprLabel& x1, ExprLabel& x2, const ExprLabel& y) {y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	inline void sub_bwd   (const ExprSub&,    ExprLabel& x1, ExprLabel& x2, const ExprLabel& y) {y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	inline void div_bwd   (const ExprDiv&,    ExprLabel& x1, ExprLabel& x2, const ExprLabel& y) {y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	       void max_bwd   (const ExprMax&,    ExprLabel& x1,   ExprLabel& x2,   const ExprLabel& y){y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	       void min_bwd   (const ExprMin&,    ExprLabel& x1,   ExprLabel& x2,   const ExprLabel& y){y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	       void atan2_bwd (const ExprAtan2& , ExprLabel& x1,   ExprLabel& x2,   const ExprLabel& y){y.taylor_comps_aff->clear();x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear(); }
	inline void minus_bwd (const ExprMinus& , ExprLabel& x,                 const ExprLabel& y){y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
        inline void trans_V_bwd(const ExprTrans&, ExprLabel& x,             const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
        inline void trans_M_bwd(const ExprTrans&, ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	       void sign_bwd  (const ExprSign&,   ExprLabel& x,                 const ExprLabel& y){y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	       void abs_bwd   (const ExprAbs&,    ExprLabel& x,                 const ExprLabel& y){y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void power_bwd(const ExprPower& e, ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void sqr_bwd   (const ExprSqr& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void sqrt_bwd  (const ExprSqrt& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void exp_bwd   (const ExprExp& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void log_bwd   (const ExprLog& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void cos_bwd   (const ExprCos& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void sin_bwd   (const ExprSin& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void tan_bwd   (const ExprTan& ,   ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void cosh_bwd  (const ExprCosh& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void sinh_bwd  (const ExprSinh& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void tanh_bwd  (const ExprTanh& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void acos_bwd  (const ExprAcos& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void asin_bwd  (const ExprAsin& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void atan_bwd  (const ExprAtan& ,  ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void acosh_bwd (const ExprAcosh& , ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void asinh_bwd (const ExprAsinh& , ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}
	inline void atanh_bwd (const ExprAtanh& , ExprLabel& x,                 const ExprLabel& y) {y.taylor_comps_aff->clear(); x.taylor_comps_aff->clear();}

	inline void add_V_bwd (const ExprAdd&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void add_M_bwd (const ExprAdd&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_SV_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_SM_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_VV_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_MV_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_MM_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void mul_VM_bwd(const ExprMul&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void sub_V_bwd (const ExprSub&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}
	inline void sub_M_bwd (const ExprSub&, ExprLabel& x1, ExprLabel& x2,    const ExprLabel& y) {y.taylor_comps_aff->clear(); x1.taylor_comps_aff->clear();x2.taylor_comps_aff->clear();}


};

inline void Univ_Taylor_Aff::index_fwd(const ExprIndex& ind, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "index = " << ind.index << endl;
#endif

    int index = ind.index;
    unsigned int i;
    if (x.af2->dim.is_scalar())
        y.af2->i() = x.af2->i();
    else
        y.af2->i() = x.af2->v().operator[](index);
    for (i = 0; i < x.taylor_comps_aff->size(); i++) {
        Affine2Domain af2(Dim::scalar());
        if (x.taylor_comps_aff->operator[](i).dim.is_scalar())
            af2.i() = x.taylor_comps_aff->operator[](i).i();
        else
            af2.i() = x.taylor_comps_aff->operator[](i).v()[index];
        y.taylor_comps_aff->push_back(af2);
    }

#ifdef __PRINT_DEBUG__
    cout << "result (" << __func__ << ") = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}


inline void Univ_Taylor_Aff::cst_fwd(const ExprConstant& c, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "cst = " << *c.deco.af2 << endl;
    if (c.deco.taylor_comps_aff) {
        cout << "number of taylor components in c : " << c.deco.taylor_comps_aff->size() << endl;
    }
    cout.flush();
    if (y.taylor_comps_aff)
        cout << "number of taylor components in y : " << y.taylor_comps_aff->size() << endl;
    cout.flush();
#endif

    if (y.taylor_comps_aff->size() == 0) {
        y.taylor_comps_aff->push_back(*c.deco.af2);
        if (c.deco.taylor_comps_aff) {
            for (unsigned int i = 1; i < c.deco.taylor_comps_aff->size(); i++) {
                y.taylor_comps_aff->push_back(c.deco.taylor_comps_aff->operator[](i));
            }
        }
    }

#ifdef __PRINT_DEBUG__
    cout << "result (" << __func__ << ") = " << endl;
    cout.flush();
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
    cout.flush();
#endif
}

inline void Univ_Taylor_Aff::symbol_fwd(const ExprSymbol& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "x = " << *x.deco.af2 << endl;
#endif

    if (y.taylor_comps_aff->size() == 0) {
        y.taylor_comps_aff->push_back(*x.deco.af2);
        if (strcmp(_t->name,x.name)==0) {
            Affine2Domain af2(x.deco.af2->dim);
            if (af2.dim.is_scalar())
                af2.i() = 1;
            else
                af2.v().init(1.0);
            y.taylor_comps_aff->push_back(af2);
        }
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif

}


inline void Univ_Taylor_Aff::apply_fwd(const ExprApply&appl, ExprLabel** x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    for (int j = 0; j < appl.nb_args; j++) {
        for (unsigned int i = 0; i < x[j]->taylor_comps_aff->size(); i++) {
            Affine2Domain af2(x[j]->taylor_comps_aff->operator[](i).dim);
            if (af2.dim.is_scalar())
                af2.i() = x[j]->taylor_comps_aff->operator[](i).i();
            else
                af2.v() = x[j]->taylor_comps_aff->operator[](i).v();
            appl.func.args()[j].deco.taylor_comps_aff->push_back(af2);
        }
    }
    //appl.func.forward<Univ_Taylor>(*this);
    compute_taylor_aff(appl.func, *y.taylor_comps_aff);

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::chi_fwd(const ExprChi&, const ExprLabel& x1, const ExprLabel& x2, const ExprLabel& x3, ExprLabel& y) {
    cerr << __func__ << endl;
}


inline void Univ_Taylor_Aff::add_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    unsigned int i = 0;
    while (i < x1.taylor_comps_aff->size() && i < x2.taylor_comps_aff->size()) {//We have two operands for the addition
        Affine2Domain af2(x2.af2->dim);
        af2.i() = x1.taylor_comps_aff->operator[](i).i() + x2.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
	i++;

    }
    while (i < x1.taylor_comps_aff->size()) {//Only elements of x1 remains, we add them to the list of Taylor components of y
        Affine2Domain  af2(x1.af2->dim);
        af2.i() = 0;
        af2.i() = x1.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
	i++;
    }
    while (i < x2.taylor_comps_aff->size()) {//Only elements of x2 remains, we add them to the list of Taylor components of y
        Affine2Domain af2(x2.af2->dim);
        af2.i()=0;
        af2.i() = x2.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}


inline void Univ_Taylor_Aff::mul_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "x1 :" << endl;
    print_taylor_comps(x1.taylor_comps_aff);
    cout << "\nx2 :" << endl;
    print_taylor_comps(x2.taylor_comps_aff);
    cout << endl;
#endif

    //First taylor component is just the evaluation of the multiplication
    Affine2Domain  y0(x2.af2->dim);
    y0.i() = x1.taylor_comps_aff->operator[](0).i()*x2.taylor_comps_aff->operator[](0).i();
    y.taylor_comps_aff->push_back(y0);

    unsigned int k = 1;
    while (k <= _max_order && k < (x1.taylor_comps_aff->size() + x2.taylor_comps_aff->size() - 1)) {//Second test is to create only  the required number of taylor components
        //computation of the k^th taylor component of y
        Affine2Domain yk(x2.af2->dim);
        yk.i()= 0;
        unsigned int j = 0;
        while(j <= k && j < x1.taylor_comps_aff->size()) {

            if (k-j < x2.taylor_comps_aff->size()) {
                yk.i() += x1.taylor_comps_aff->operator[](j).i()*x2.taylor_comps_aff->operator[](k-j).i();
            } else {
	    }
            j++;
        }
        y.taylor_comps_aff->push_back(yk);
        k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}


inline void Univ_Taylor_Aff::sub_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    unsigned int i = 0;
    while (i < x1.taylor_comps_aff->size() && i < x2.taylor_comps_aff->size()) {
        Affine2Domain af2(x2.af2->dim);
        af2.i() = x1.taylor_comps_aff->operator[](i).i() - x2.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
	i++;
    }
    while (i < x1.taylor_comps_aff->size()) {
        Affine2Domain  af2(x1.af2->dim);
        af2.i() = x1.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }
    while (i < x2.taylor_comps_aff->size()) {
        Affine2Domain af2(x2.af2->dim);
        af2.i()=0;
        af2.i() = - x2.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::div_fwd(const ExprDiv&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    //Computation of the first taylor component y0 that is the evaluation of the division
    Affine2Domain y0(x1.af2->dim);
    y0.i() = x1.taylor_comps_aff->operator[](0).i() / x2.taylor_comps_aff->operator[](0).i();
    y.taylor_comps_aff->push_back(y0);
    unsigned int k = 1;
    Affine2Domain coeff(y0.dim);
    coeff.i() = 1.0 / ( x2.taylor_comps_aff->operator[](0).i() );
    //computation of the yk other Taylor components
    while (k <= _max_order) {
        Affine2Domain  yk(x1.af2->dim);
        yk.i() = coeff.i();
        Affine2Domain accu(x1.af2->dim);
        accu.i() = 0;
        unsigned int j;
        for (j = 0; j < k; j++) {
            if ((k-j) < x2.taylor_comps_aff->size()) {
                accu.i() += y.taylor_comps_aff->operator[](j).i()*x2.taylor_comps_aff->operator[](k-j).i();
            }
        }
        if (k < x1.taylor_comps_aff->size()) {//x1 has a k^th Taylor component
            yk.i() *= (x1.taylor_comps_aff->operator[](k).i() - accu.i());
        } else {
            yk.i() *= 0.0 - accu.i();
        }
        y.taylor_comps_aff->push_back(yk);
	k++;
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::minus_fwd(const ExprMinus&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    unsigned int i = 0;
    while (i < x.taylor_comps_aff->size()) {
        Affine2Domain af2(x.af2->dim);
        af2.i() = - x.taylor_comps_aff->operator[](i).i();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}


inline void Univ_Taylor_Aff::power_fwd(const ExprPower& p, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    //First taylor component is just the evaluation of the multiplication
    Affine2Domain y0(x.af2->dim);
    y0.i() = 1;
    for (int i = 0; i < p.expon; i++) y0.i() *= x.taylor_comps_aff->operator[](0).i();
    y.taylor_comps_aff->push_back(y0);
    Affine2Domain coeff(y0.dim);
    coeff.i() = 1.0 / ( x.taylor_comps_aff->operator[](0).i() );
    unsigned int k = 1;
    //computation of the k^th taylor component of y
    while (k <= _max_order) {
        Affine2Domain yk(x.af2->dim);
        yk.i() = coeff.i();
        Affine2Domain accu1(x.af2->dim);
        accu1.i() = 0;
        unsigned int j;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size() ; j++) {
            accu1.i() += y.taylor_comps_aff->operator[](k-j).i() * j * x.taylor_comps_aff->operator[](j).i();
        }
        accu1.i() *= p.expon;

        for (j = 1; j < k ; j++) {
            if (k-j < x.taylor_comps_aff->size())
                accu1.i() -= x.taylor_comps_aff->operator[](k-j).i() * j * y.taylor_comps_aff->operator[](j).i();
        }
        yk.i() *= accu1.i();
        yk.i() = yk.i() / k;
        y.taylor_comps_aff->push_back(yk);
	k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::sqr_fwd(const ExprSqr&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    //First taylor component is just the evaluation of the multiplication
    Affine2Domain  y0(x.af2->dim);
    y0.i() = sqr(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(y0);
    unsigned int k = 1;
    while (k <= _max_order && k < (2*x.taylor_comps_aff->size() - 1)) {
        //computation of the k^th taylor component of y
        Affine2Domain yk(x.af2->dim);
        yk.i()= 0;
        unsigned int j = 0;
        while(j <= k && j < x.taylor_comps_aff->size()) {
            if (k-j < x.taylor_comps_aff->size()) {
                yk.i() += x.taylor_comps_aff->operator[](j).i()*x.taylor_comps_aff->operator[](k-j).i();
            }
            j++;
        }
        y.taylor_comps_aff->push_back(yk);
        k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::sqrt_fwd(const ExprSqrt&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    Affine2Domain  y0(x.af2->dim);
    y0.i() = sqrt(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(y0);

    unsigned int k = 1;
    Affine2Domain coeff(y0.dim);
    coeff.i() = 1.0 / ( 2 * y0.i() );
    while (k <= _max_order) {
        Affine2Domain  yk(x.af2->dim);
        yk.i() = coeff.i();
        Affine2Domain accu(x.af2->dim);
        accu.i() = 0;
        unsigned int j;
        if (k > 1) {
            for (j = 1; j < k; j++) {
                accu.i() += y.taylor_comps_aff->operator[](j).i()*y.taylor_comps_aff->operator[](k-j).i();
            }
        }
        if (k < x.taylor_comps_aff->size()) {
            yk.i() *= x.taylor_comps_aff->operator[](k).i() - accu.i();
        } else {
            yk.i() *= 0.0 - accu.i();
        }
        y.taylor_comps_aff->push_back(yk);
	k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::exp_fwd(const ExprExp&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    //calculation of y0:
    Affine2Domain  y0(x.af2->dim);
    y0.i() = exp(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(y0);
    unsigned int k = 1;
    while (k <= _max_order) {
        Affine2Domain  yk(x.af2->dim);
        yk.i() = 0;
        unsigned int j;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size() ; j++) {
            yk.i() += y.taylor_comps_aff->operator[](k-j).i() * j * x.taylor_comps_aff->operator[](j).i();
        }
        yk.i() = yk.i() / k;
        y.taylor_comps_aff->push_back(yk);
        k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::log_fwd(const ExprLog&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    //calculation of y0:
    Affine2Domain  y0(x.af2->dim);
    y0.i() = log(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(y0);


    //Calculation of the coefficient 1/x0
    Affine2Domain coeff(y0.dim);
    coeff.i() = 1.0 / ( x.taylor_comps_aff->operator[](0).i() );
    //Calculation of the yk:
    unsigned int k = 1;
    while (k <= _max_order) {
        Affine2Domain yk(x.af2->dim);
        yk.i() = coeff.i();

        //accu get the result of the sum in the expression of the log
        Affine2Domain accu(x.af2->dim);
        accu.i() = 0;
        unsigned int j;
        for (j = 1; j < k ; j++) {
            if (k-j < x.taylor_comps_aff->size())
                accu.i() += x.taylor_comps_aff->operator[](k-j).i() * j * y.taylor_comps_aff->operator[](j).i();
        }
        if (k < x.taylor_comps_aff->size()) { //xk exists
            yk.i() *= (k*x.taylor_comps_aff->operator[](k).i() - accu.i())/j;
        } else {
            yk.i() *= (0.0 - accu.i())/j;
        }
        y.taylor_comps_aff->push_back(yk);
	k++;
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::cos_fwd(const ExprCos&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    std::vector<Affine2Domain> sin_tayl;
    Affine2Domain s0(x.af2->dim);
    s0.i() = sin(x.taylor_comps_aff->operator[](0).i());
    sin_tayl.push_back(s0);
    Affine2Domain c0(x.af2->dim);
    c0.i() = cos(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(c0);
    for (unsigned int k = 1; k <= _max_order; k++) {
        unsigned int j;
        Affine2Domain ck(x.af2->dim);
        ck.i() = 0;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size(); j++) {
            ck.i() -= j * x.taylor_comps_aff->operator[](j).i() * sin_tayl[k-j].i();
        }
        ck.i() = ck.i() / k;
        y.taylor_comps_aff->push_back(ck);
        Affine2Domain sk(x.af2->dim);
        sk.i() = 0;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size(); j++) {
            sk.i() += j * x.taylor_comps_aff->operator[](j).i() * y.taylor_comps_aff->operator[](k-j).i();
        }
        sk.i() = sk.i() / k;
        sin_tayl.push_back(sk);
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}


inline void Univ_Taylor_Aff::sin_fwd(const ExprSin&, const ExprLabel& x, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    std::vector<Affine2Domain> cos_tayl;
    Affine2Domain s0(x.af2->dim);
    s0.i() = sin(x.taylor_comps_aff->operator[](0).i());
    y.taylor_comps_aff->push_back(s0);
    Affine2Domain c0(x.af2->dim);
    c0.i() = cos(x.taylor_comps_aff->operator[](0).i());
    cos_tayl.push_back(c0);
    for (unsigned int k = 1; k <= _max_order; k++) {
        unsigned int j;
        Affine2Domain ck(x.af2->dim);
        ck.i() = 0;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size(); j++) {
            ck.i() -= j * x.taylor_comps_aff->operator[](j).i() * y.taylor_comps_aff->operator[](k-j).i();
        }
        ck.i() = ck.i() / k;
        cos_tayl.push_back(ck);
        Affine2Domain sk(x.af2->dim);
        sk.i() = 0;
        for (j = 1; j <= k && j < x.taylor_comps_aff->size(); j++) {
            sk.i() += j * x.taylor_comps_aff->operator[](j).i() * cos_tayl[k-j].i();
        }
        sk.i() = sk.i() / k;
        y.taylor_comps_aff->push_back(sk);
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::add_V_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
#endif

    unsigned int i = 0;
    while (i < x1.taylor_comps_aff->size() && i < x2.taylor_comps_aff->size()) {
        Affine2Domain af2(x2.af2->dim);
        af2.v() = x1.taylor_comps_aff->operator[](i).v() + x2.taylor_comps_aff->operator[](i).v();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }
    while (i < x1.taylor_comps_aff->size()) {
        Affine2Domain af2(x1.af2->dim);
        af2.v() = x1.taylor_comps_aff->operator[](i).v();
        y.taylor_comps_aff->push_back(af2);
        i++;

    }
    while (i < x2.taylor_comps_aff->size()) {
        Affine2Domain af2(x2.af2->dim);
        af2.v().clear();
        af2.v() = x2.taylor_comps_aff->operator[](i).v();
        y.taylor_comps_aff->push_back(af2);
        i++;
    }

#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::add_M_fwd(const ExprAdd&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cerr << __func__ << endl;
}

inline void Univ_Taylor_Aff::mul_SV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    print_taylor_comps(x2.taylor_comps_aff);
    print_taylor_comps(x1.taylor_comps_aff);
#endif

    //First taylor component is just the evaluation of the multiplication
    Affine2Domain y0(x2.af2->dim);
    y0.v() = x1.taylor_comps_aff->operator[](0).i()*x2.taylor_comps_aff->operator[](0).v();
    y.taylor_comps_aff->push_back(y0);
    unsigned int k = 1;
    while (k <= _max_order && k < (x1.taylor_comps_aff->size() + x2.taylor_comps_aff->size() - 1)) {
        //computation of the k^th taylor component of y
        Affine2Domain yk(x2.af2->dim);
        yk.v().clear();
        unsigned int j = 0;
        while(j <= k && j < x1.taylor_comps_aff->size()) {
	    if (k-j < x2.taylor_comps_aff->size()) {
	        yk.v() += x1.taylor_comps_aff->operator[](j).i()*x2.taylor_comps_aff->operator[](k-j).v();
            }
            j++;
        }
        y.taylor_comps_aff->push_back(yk);
        k++;
    }
#ifdef __PRINT_DEBUG__
    cout << "result = " << endl;
    print_taylor_comps(y.taylor_comps_aff);
    cout << endl;
#endif
}

inline void Univ_Taylor_Aff::max_fwd(const ExprMax&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::min_fwd(const ExprMin&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::atan2_fwd(const ExprAtan2&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cerr << __func__ << endl;
}

inline void Univ_Taylor_Aff::trans_V_fwd(const ExprTrans&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::trans_M_fwd(const ExprTrans&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::sign_fwd(const ExprSign&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::abs_fwd(const ExprAbs&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}

inline void Univ_Taylor_Aff::tan_fwd(const ExprTan&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::cosh_fwd(const ExprCosh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::sinh_fwd(const ExprSinh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::tanh_fwd(const ExprTanh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::acos_fwd(const ExprAcos&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::asin_fwd(const ExprAsin&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::atan_fwd(const ExprAtan&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::acosh_fwd(const ExprAcosh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::asinh_fwd(const ExprAsinh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}
inline void Univ_Taylor_Aff::atanh_fwd(const ExprAtanh&, const ExprLabel& x, ExprLabel& y) {
    cerr << __func__ << endl;
}

inline void Univ_Taylor_Aff::mul_SM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::mul_VV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::mul_MV_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::mul_VM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::mul_MM_fwd(const ExprMul&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::sub_V_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}
inline void Univ_Taylor_Aff::sub_M_fwd(const ExprSub&, const ExprLabel& x1, const ExprLabel& x2, ExprLabel& y) {
    cout << __func__ << endl;
}



}//namespace ibex

#endif // UNIV_TAYLOR_AFF_H
