/* ============================================================================
 * D Y N I B E X - Definition of the weighted dags used for LTE
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef WDAG_AFF_H
#define WDAG_AFF_H

//decomment if you want "verbose style" of the computation
// #ifndef __PRINT_DEBUG__
// #define __PRINT_DEBUG__
// #endif

#include <list>

#include "ibex_vertex_wdag_aff.h"

#include <fstream>

#include <map>

#include <unistd.h>

using namespace std;

namespace ibex{

class vertex_wDAG_aff;

/**
 * \brief DAGs of rooted trees : encapsulates vertex_wDAG_aff to represent a dag of _max_order
 */
class wDAG_aff
{

private:
  list<rootedTree> * _list_of_trees;//All trees of order <= _max_order
  list<list<ibex::vertex_wDAG_aff *> > * _list_of_dags;//list of size _max_order. Each element i is a list of vertex_Dag of order i+1 with i in [0,_max_order-1]
  int _max_order;//Max order of the generated trees.
  Function * _f; //Function for which we want to compute values.
  IntervalMatrix * _A;//Matrix A in the Butcher tableau
  IntervalVector * _b;//Vector b in the Butcher tableau
  bool _filled_in;//true if A and b have been filled in.
  vector<Function *> * _list_g_taylor;
  unsigned int _flag;//flag to know if the computation of delta in the dag is required

  /**
   * \brief for internal contruction of _list_of_trees.
   */
  void add_root(vertex_wDAG_aff * root);
public:
  ~wDAG_aff(){
//     cout << __func__ << endl;cout.flush();
    if (_list_of_trees) delete _list_of_trees;
    if (_list_of_dags) {
      delete _list_of_dags;
    }
    if (_A) delete _A;
    if (_b) delete _b;
  };
  /**
   * \brief Constructor for when _A and _b have not to be set (for example to perform a Taylor).
   */
  wDAG_aff(Function* f, int order);

  /**
   * \brief Complete constructor
   */
  wDAG_aff(Function * f, int order, IntervalMatrix & A, IntervalVector & b);

  /**
   * \brief getter of the list of dags.
   */
  inline list<list<vertex_wDAG_aff *> > * get_list_of_dags(){return _list_of_dags;};

  /**
   * \brief getter of the max order.
   */
  inline int get_order(){return _max_order;};

  /**
   * \brief return the list of dags representing the trees of order=order
   */
  std::list< vertex_wDAG_aff* > * get_dag(int order);


  /**
   * To know if the Butcher table has been filled in (for internal purpose)
   */
  bool is_filled_in(){return _filled_in;};

  void set_filled_in(bool filled_in){_filled_in = filled_in;};


  /**
   * \brief setter for the interval matrix _A.
   */
  inline void set_A(IntervalMatrix A){
    _A = new IntervalMatrix(A);
  }

  /**
   * \brief getter for the interval matrix _A.
   */
  inline IntervalMatrix * get_A(){return _A;}

  /**
   * \brief setter for the interval vector _b.
   */
  inline void set_b(IntervalVector b){
    _b = new IntervalVector(b);
  }

  /**
   * \brief getter for the interval vector _b.
   */
  inline IntervalVector * get_b(){return _b;}

  inline void set_f(Function * f){
    _f = f;
  }

  /**
   * \brief computes the order^th derivatives of f evalued at y_init. Butcher table has to be filled in priorly.
   */
  Affine2Vector get_derivatives(int order, const Affine2Vector & y_init);

  /**
   * get the local truncature error for given order and y_init. Butcher table has to be filled in priorly.
   */
  Affine2Vector get_lte(int order, const Affine2Vector & y_init);

  /**
   * \brief compute the list of rootedTree of order <= argument
   */
  static list<rootedTree> * generate(const int);

//void generate_rec(int, int, int, int, rootedTree, list<rootedTree> *);

  static void generate_dags(int max_order, list<rootedTree> * list_of_tree,
		     list< list<vertex_wDAG_aff *> > * dags, IntervalMatrix * A, IntervalVector * b);

  bool dags_to_file(const std::string &path = "dags.dat");

  bool file_to_dags(int order, list<list<ibex::vertex_wDAG_aff *> > * result, const std::string &path = "dags.dat");
};






} //namespace

#endif // WDAG_AFF_H
