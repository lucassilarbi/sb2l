/* ============================================================================
 * D Y N I B E X - Definition of the vertices of weighted dags for LTE
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef VERTEX_WDAG_AFF_H
#define VERTEX_WDAG_AFF_H

#include "ibex_rootedtree.h"
#include "ibex_Affine2Matrix.h"
#include "ibex_Function.h"

// #ifndef __PRINT_DEBUG__
// #define __PRINT_DEBUG__
// #endif

using namespace std;

namespace ibex {
class vertex_wDAG_aff
{

private:
    //tree that the vertex_wDAG_aff represents
    rootedTree * _tree;
    //list of edges starting from this. represented by the pair (son, weight)
    list<pair<vertex_wDAG_aff *, int > > * _list_of_edge;


  // Warning: factorial computation for very low value integer!!
    int factorial(int x, int accu){
      if (x <= 0)
	return 0;
      else if (x == 1)
	return accu;
      else
	return factorial(x - 1, x * accu);
    }

    /**
     * \brief compute the factorial of an integer x
     */
    inline int factorial(int x) {
        return factorial(x, 1);
    }

    //order of the tree (number of the nodes in the tree, NOT the number of son in _list_of_edge)
    int _order;

    //number of symmetries for the tree
    int _symetry;
    unsigned int _density;


    Interval * _alpha;//Number of essentially different ways of labelling the nodes of _tree.
    Interval * _psi;

    Affine2Vector * _delta;

    //flag to know when to read _delta from this or compute it because the value changed
    unsigned int _flag;

//     //Vector * component_wise_product(vector<double> & a, vector<double> & b);
//     Affine2Vector power(Affine2Vector a, int n);

    /**
     * \brief compute the sum of the components of x
     */
    static inline double sum(const Vector & x);

    /**
     * \brief compute the sum of the components of x
     */
    static inline double sum(const IntervalVector & x);

    /**
     * \brief compute the sum of the components of x
     */
    static inline Interval sum(const Affine2Vector & x);

    /**
    * \brief compute psi(this). function described in Lambert book p.167 eq. (5.45)
    */
    Interval get_psi(int i, IntervalMatrix * A, IntervalVector * b);

    /**
     * \brief internal function for the computation of gamma. Should not be called, but get_gamma instead.
     */
    static void rec_gamma(const Vector & i, const Vector & j, int d, Vector k, int pos, double & result);

    Interval power(Interval a, int b);

public:
    vertex_wDAG_aff();

    ~vertex_wDAG_aff(){
      cout << __func__ << endl;
      if(_alpha)
	delete _alpha;
      if (_psi)
	delete _psi;
      if (_delta)
	delete _delta;
    }

    vertex_wDAG_aff(rootedTree * tree);

    vertex_wDAG_aff(int order, int symetry, int density, Interval * alpha);

    //Getters
    inline unsigned int get_flag() {
        return _flag;
    }
    inline rootedTree * get_tree() {
        return _tree;
    };
    inline list<pair<vertex_wDAG_aff *, int > > * get_list_of_edge() {
        return _list_of_edge;
    };

    /**
    * \brief Return the order of the DAG (that is the number of nodes of its corresponding tree)
    **/
    int get_order();

    /**
    * \brief return the number of symetries of the rooted tree
    **/
    int get_symetry();


    unsigned int get_density();

    Interval * get_alpha();

    Interval * get_psi(IntervalMatrix * A, IntervalVector * b);

    void add_edge(vertex_wDAG_aff * vertex);

    static double binomial(double n, int k);


//     static Interval binomial(const Interval & n, const int & k);

    static double multinomial(const Vector &, const Vector &);

    Affine2Vector taylor_coeff(Function * f, const Affine2Vector & x, Affine2Vector s, int d, vector<Function *> * list_g_taylor);

    void get_partial_derivative(Function * f, const Affine2Vector & x, Affine2Matrix & S, Vector & i, Vector j, int d, int pos, int d_rest,
				Affine2Vector * _result, vector<Function *> * list_g_taylor);

    static double get_gamma(const Vector & i, const Vector & j, const int&);

    Affine2Vector * get_delta(Function * f, const Affine2Vector & x, unsigned int flag, vector<Function *> * list_g_taylor);

    Affine2Vector * get_lte(Function * f, const Affine2Vector & x);

    void set_list_of_edge(list<pair<vertex_wDAG_aff *, int > > list_of_edge){_list_of_edge = new list<pair<vertex_wDAG_aff *, int > >(list_of_edge);};
};



double vertex_wDAG_aff::sum(const Vector & x) {
//     cout << "x : " << x << endl;
    int i;
    double result = 0.0;
    for (i = 0; i < x.size(); i++)
        result += x[i];
//     cout << "sum : " << result << endl;
    return result;
}

double vertex_wDAG_aff::sum(const IntervalVector & x) {
    int i;
    double result = 0.0;
    for (i = 0; i < x.size(); i++)
        result += x[i].lb();
    return result;
}

Interval vertex_wDAG_aff::sum(const Affine2Vector & x) {

    int i;
    Interval result = Interval(0.0);
    for (i = 0; i < x.size(); i++)
        result += x[i].itv();
    return result;
}



}//namespace ibex
#endif // VERTEX_WDAG_AFF_H
