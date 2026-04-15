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

#include "ibex_vertex_wdag_aff.h"

namespace ibex {
//----------------------------------------------------------------------------------------
vertex_wDAG_aff::vertex_wDAG_aff(rootedTree * tree)
{
    _tree = tree;
    _list_of_edge = new list<pair<vertex_wDAG_aff *, int > >();
    _flag = 0;
    _density = 0;
    _order = 0;
    _symetry = 0;
    _alpha = (NULL);
    _delta = (NULL);
    _psi = (NULL);
}

/**
 * \brief add an edge to this if it is not yet in the list of edges
 **/
void vertex_wDAG_aff::add_edge(vertex_wDAG_aff* vertex)
{
    if (_list_of_edge->size() == 0) { //This was a leaf
	pair<vertex_wDAG_aff *, int > temp(vertex, 1);
        _list_of_edge->push_back(temp);
    } else {
        bool found = false;
        list<pair<vertex_wDAG_aff *, int > >::iterator it_edge = _list_of_edge->begin();
        while (!found && it_edge != _list_of_edge->end()) {
            if ((*it_edge).first == vertex) { //vertex is already a child of this
                found = true;
                (*it_edge).second += 1;
            }
            it_edge++;
        }
        if (!found) { //vertex was not a child of this
            pair<vertex_wDAG_aff *, int > temp(vertex, 1);
	    _list_of_edge->push_back(temp);
        }
    }
}

vertex_wDAG_aff::vertex_wDAG_aff(int order, int symetry, int density, Interval * alpha){
    _list_of_edge = new list<pair<vertex_wDAG_aff *, int > >();
    _order = order;
    _density = density;
    _symetry = symetry;
    _flag = 0;
    _alpha = new Interval(*alpha);
    _delta = (NULL);
    _psi = (NULL);
}


//----------------------------------------------------------------------------------------

/**
 * \brief Return the order of the DAG (that is the number of nodes of its corresponding tree)
 **/
int vertex_wDAG_aff::get_order()
{
  if (_order == 0) {
    _order = 1;
    list<pair<vertex_wDAG_aff*, int> >::iterator it;
    for (it = _list_of_edge->begin(); it != _list_of_edge->end(); it++) {
      _order += (*it).second * (*it).first->get_order();
    }
  }
  return _order;
}


/**
 * \brief return the number of symetries of the rooted tree
 **/
int vertex_wDAG_aff::get_symetry()
{
  if (_symetry == 0) {
    _symetry = 1;
    list<pair<vertex_wDAG_aff*, int> >::iterator it;
    for (it = _list_of_edge->begin(); it != _list_of_edge->end(); it++) {
	_symetry *= factorial((*it).second) * std::pow((*it).first->get_symetry(),(*it).second);
    }
  }
  return _symetry;
}

/**
 * \brief return the density of the rooted tree
 **/
unsigned int vertex_wDAG_aff::get_density()
{
  if (_density == 0) {
    _density = this->get_order();
    list<pair<vertex_wDAG_aff*, int> >::iterator it;
    for (it = _list_of_edge->begin(); it != _list_of_edge->end(); it++) {
	_density *= std::pow((*it).first->get_density(), (*it).second);
    }
  }
  return _density;
}

Interval * vertex_wDAG_aff::get_alpha()
{
  if (_alpha==NULL) {
    _alpha = new Interval(Interval(factorial(get_order()))/Interval(1.0*get_symetry()*get_density()));
  }
  return _alpha;
}

Interval vertex_wDAG_aff::power(Interval a, int b){
  if (b == 0)
      return Interval(1.);
  Interval result(a);
  for (int i = 1; i < b; i++){
    result*= a;
  }
  return result;
}

/**
 * \brief computes psi_i(this). function described in Lambert book p.167 eq. (5.45)
 */
Interval vertex_wDAG_aff::get_psi(int i, IntervalMatrix * A, IntervalVector * b) {
    Interval  result(0.0);
    for (int j = 0; j < A->nb_rows(); j++) {
      Interval temp( (*A)[i][j] );
      list<pair<vertex_wDAG_aff*, int> >::iterator it;
      for (it = _list_of_edge->begin(); it != _list_of_edge->end() && (temp.lb() != 0 || temp.ub()!=0); it++) {
	temp *= power((*it).first->get_psi(j, A, b),(*it).second);
      }
      result += temp;
    }
    return result;
}

/**
 * \brief compute psi(this). function described in Lambert book p.167 eq. (5.45)
 */
Interval* vertex_wDAG_aff::get_psi(IntervalMatrix * A, IntervalVector * b)
{
  if (_psi != NULL) {
    return _psi;
  }

  Interval  result(0.0);
  for (int j = 0; j < A->nb_rows(); j++) {
    Interval temp(b->operator[](j));
    list<pair<vertex_wDAG_aff*, int> >::iterator it;
    for (it = _list_of_edge->begin(); it != _list_of_edge->end() && (temp.lb() != 0 || temp.ub()!=0); it++) {
      temp *= power((*it).first->get_psi(j, A, b),(*it).second);
    }
    result+= temp;
  }
  _psi = new Interval(result);
  return _psi;
}

// /**
//  * \brief return the value of a^n
//  **/
// Affine2Vector vertex_wDAG_aff::power(Affine2Vector a, int n) {
//   cout << a << endl;
//   Affine2Vector result(a.size());
//   for (int i = 0; i < a.size(); i++){
//     Affine2 toto(1.0);
//     for(int k = 0; k < n; k++) {
//       toto *= a[i];
//     }
//     result[i] = toto;
//   }
//   return result;
// }

double vertex_wDAG_aff::binomial(double n, int k){
  double res = 1.0;
  for (int i = 0; i < k; i++){
      res*= (n-i)/(i+1);
  }

  return res;
}


// double vertex_wDAG_aff::binomial(double n, int k){
//     cout << "n = " << n << endl;
//     cout << "k = " << k << endl;
//   if (k < 0 /*|| k > n*/){
//       cout << "result LA : 0" << endl;
//     return 0.0;
//   }
//   if (k == 0 || k == n){
//       cout << "result ICI : 1" << endl;
//     return 1.0;
//   }
//   k = std::min<double>(k, n-k);
//   double c = 1.0;
//   for (int i = 0; i < k; i++){
//     c = c*(n - i)/((double)i + 1);
//   }
//   cout << "result : " << c << endl;
//   return c;
// }


// double vertex_wDAG_aff::binomial(int n, int k){
//     double result = 1.0;
//     for (int i = k; i > 0 && result != 0; i--){
//         result *= ((double)n-i + 1)/(double)i;
//     }
//     cout << "result = " << result << endl;
//     return result;
// }

// double vertex_wDAG_aff::binomial(const int & n, const int & k) {
// #ifdef __PRINT_DEBUG__
//   cout << "******************************************************************************************" << endl;
//   cout << __func__ << endl;
//   cout << "n = " << n << "; k = " << k << endl;
// #endif
//
//   if (n == k)
//       return 1;
//   if (k == 0)
//       return 1;
//   if (n < k)
//       return 0;
//   if (n == 0)
//       return 0;
//   if (k == 1)
//       return n;
//
// #ifdef __PRINT_DEBUG__
//   double result = ((n - k + 1)/(1.0*k)) * binomial2(n, k-1);
//   cout << "result = " << result << endl;
//   return result;
// #else
//   return (n - k + 1)/(1.0*k) * binomial(n, k-1);
// #endif
// }

//----------------------------------------------------------------------------------------
// Interval vertex_wDAG_aff::binomial(const Interval & n, const int & k) {
// #ifdef __PRINT_DEBUG__
//     cout << "******************************************************************************************" << endl;
//     cout << __func__ << endl;
//     cout << "n = " << n << "; k = " << k << endl;
// #endif
//
//     return binomial((int)n.lb(), k);
//
//
// //     if (n == k)
// //         return 1;
// //     if (k == 0)
// //         return 1;
// //     if (n.ub() < k)
// //         return 0;
// //     if (n == 0)
// //         return 0;
// //     if (k == 1)
// //         return n;
// // #ifdef __PRINT_DEBUG__
// //     Interval result = ((n - k + 1)/k) * binomial(n, k-1);
// //     cout << "result = " << result << endl;
// //     return result;
// // #else
// //     return (n - k + 1)/k * binomial(n, k-1);
// // #endif
// }

double vertex_wDAG_aff::multinomial(const Vector& i, const Vector& j) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "i = " << i << "\nj = " << j << endl;
#endif

    double result = binomial(i[0], (int)j[0]);
    int l;
    for (l = 1; l < j.size() && (result !=0); l++) {
        result *= binomial(i[l], (int)j[l]);
    }

#ifdef __PRINT_DEBUG__
    cout << "fin de " <<__func__ << endl;
    cout << "result = " << result << endl;
#endif
    return result;
}


// void vertex_wDAG_aff::rec_gamma(const Vector & i, const Vector & j, int d, Vector k, int pos, Affine2 & result) {
// #ifdef __PRINT_DEBUG__
//     cout << "******************************************************************************************" << endl;
//     cout << __func__ << endl;
//     cout << "i = " << i << "\nj = " << j << endl << "\nk = " << k << endl << "d = " << d << endl << "\npos = " << pos << endl << "\nresult = " << result << endl;
// #endif
//     if (pos == i.size()-1) {
//         Affine2Vector i_inter = Affine2Vector(i);
//         Affine2Vector k_inter = Affine2Vector(k);
//         while (k[pos] <= i[pos]) {
//             if (sum(k) != 0) {
// 	      Passer en intervalles!!!
//
// 	      Affine2 temp = std::pow(-1.0, sum(i - k)) * multinomial(i_inter, k) * multinomial((d/ sum(k_inter))*k_inter, j) * pow((sum(k_inter)/(1.0*d)), sum(i));
//
//               Griewank :
// 	      Affine2 temp = pow(-1.0, (int)sum(i - k)) * multinomial(i_inter, k) * multinomial((sum(i)/ sum(k_inter))*k_inter, j) * pow(sum(k_inter)/(1.0*d), sum(i));
//
// 	      Bartha
// 	      Affine2 temp = pow(-1.0, (int)sum(i - k)) * multinomial(IntervalVector(j), k) * multinomial((sum(j)/ sum(k_inter))*k_inter, j) * pow(sum(k_inter)/(1.0*sum(j)), (int)sum(i));
//
//
//                 result += temp;
//             }
//             k[pos] = k[pos] + 1;
//             k_inter[pos] = k_inter[pos] + 1;
//         }
//     } else {
//         while (k[pos] <= i[pos]) {
//             rec_gamma(i, j, d, k, pos+1, result);
//             k[pos] = k[pos] + 1;
//         }
//     }
// }

void vertex_wDAG_aff::rec_gamma(const Vector & i, const Vector & j, int d, Vector k, int pos, double & result) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "i = " << i << "\nj = " << j << endl << "\nk = " << k << endl << "d = " << d << endl << "\npos = " << pos << endl << "\nresult = " << result << endl;
#endif
//     cout << "\tk = " << k << endl;
  if (pos == i.size()-1) {
//       cout << "i = " << i << "\nj = " << j << endl << "\nk = " << k << endl << "d = " << d << endl << "\npos = " << pos << endl << "\nresult = " << result << endl;
//     IntervalVector i_inter = IntervalVector(i);
//     IntervalVector k_inter = IntervalVector(k);
    while (k[pos] <= i[pos]) {
      if (sum(k) != 0) {
//           cout << "k = " << k << endl;
//           cout << "i = " << i << endl;
//           cout << "d = " << d << endl;
//           cout << "==================" << endl;
//           cout << "(1.0/ sum(k)*(1.0*d)*k) = " << (1.0/ sum(k)*(1.0*d)*k) << endl;
//           cout << "\tk = " << k << endl;
	  double temp = std::pow(-1.0, (int)sum(i - k)) * multinomial(i, k) * multinomial((1.0/ sum(k)*(1.0*d)*k), j) * std::pow(sum(k)/(1.0*d), (int)sum(i));
	  result += temp;
      }
      k[pos] = k[pos]+1;

//       k_inter[pos] = k_inter[pos] + 1;
    }
  } else {
    while (k[pos] <= i[pos]) {
      rec_gamma(i, j, d, k, pos+1, result);
      k[pos] = k[pos] + 1;

    }
  }
}

double vertex_wDAG_aff::get_gamma(const Vector & i, const Vector & j, const int& d) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "i = " << i << "\nj = " << j << endl << "d = " << d << endl;
#endif
    double result = 0.0;
    Vector k = Vector::zeros(i.size());
//     cout << "i = " << i << endl;
    rec_gamma(i, j, d, k, 0, result);
#ifdef __PRINT_DEBUG__
    cout << "result (" << __func__ << ") = " << result << endl;
#endif
    return result;
}

Affine2Vector vertex_wDAG_aff::taylor_coeff(Function * f, const Affine2Vector & x, Affine2Vector s, int d, vector<Function *> * list_g_taylor) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "f : " << *f << "\nx = " << x << endl;
    cout.flush();
    cout << "s = " << s << endl;
    cout << "s.itv() = " << s.itv() << endl;
    cout << "d = " << d << endl;
    cout << "x.size() : " << x.size() << endl;
    cout << "s.size() : " << s.size() << endl;
#endif

  Affine2Vector result(f->image_dim());
  result.clear();
  for (int i = 0; i < f->image_dim(); i++) {
    Function * f_temp = (list_g_taylor->operator[](i));
    Affine2Vector t_inst(1+x.size()+s.size());
    t_inst[0] = 0.0;
    int j;
    for (j = 1; j <= x.size(); j++)
      t_inst[j] = x[j-1];
    while (j < t_inst.size()){
      t_inst[j] = s[j-x.size()-1];
      j++;
    }

    //muliplication with variable and affine2vector not allowed (reason why cst_fwd is not called)
    std::vector<Affine2Domain> taylor_comps;
    f_temp->taylor_aff(t_inst, taylor_comps, d, &(f_temp->arg(0)));
    if (d >= taylor_comps.size()) {
	result[i] = 0.0;
    } else {
	result[i] = taylor_comps[d].i();
    }
  }
#ifdef __PRINT_DEBUG__
  cout << "fin de " <<__func__ << endl;
  cout << "result (" << __func__ << ") = " << result << endl;
#endif

  return result;
}

void vertex_wDAG_aff::get_partial_derivative(Function * f, const Affine2Vector & x, Affine2Matrix & S, Vector & i, Vector j, int d, int pos, int d_rest,
					     Affine2Vector * _result, vector<Function *> * list_g_taylor) {
#ifdef __PRINT_DEBUG__
  cout << "******************************************************************************************" << endl;
  cout << __func__ << endl;
  cout << "f : " << *f << endl;
  cout << "x = " << x << endl;
  cout.flush();
  cout << "S = " << endl <<  S << endl;
  cout << "i = " << i << endl;
  cout.flush();
  cout << "j = " << j << endl;
  cout.flush();
  cout << "d = " << d << endl;
  cout << "pos = " << pos << endl;
  cout << "d_rest = " << d_rest << endl;
  cout << "result = " << *_result << endl;
  cout << get_tree()->to_string() << endl;
#endif


  if (pos == j.size()-1) {
    j[pos] = d_rest;

    Interval gamma = get_gamma(i, j, d);
//     cout << "gamma : " << gamma << endl;
    if (gamma.lb() !=0 || gamma.ub() !=0){//if gamma != 0
//         cout << _tree->to_string() << endl;
//         cout << "\tj : " << j << endl;
//         cout << "\tgamma = " << gamma << endl;

      Affine2Vector Sj(S*j);
      Affine2Vector tayl_coef = taylor_coeff(f, x, Sj, (int)sum(i), list_g_taylor);
//       cout << "\ttaylor coeff : " << tayl_coef << endl;
      *_result += gamma * tayl_coef;
    } else {
      cout << "\tj : " << j << endl;
//       cout << "\tgamma = 0 " << endl;
    }
  } else {
    while (j[pos] <= d_rest) {
      get_partial_derivative(f, x, S, i, j, d, pos + 1, d_rest - (int)j[pos], _result, list_g_taylor);
      j[pos]++;
    }
  }

#ifdef __PRINT_DEBUG__
  cout << "fin de " <<__func__ << endl;
  cout << "result (" << __func__ << ") = " << *_result << endl;
#endif
}

/**
 * \brief computes delta(tree); Cf. Computing of B-series by automatic differentiation (Bartha et al.), eq. 2.10
 */
Affine2Vector * vertex_wDAG_aff::get_delta(Function * f, const Affine2Vector & x, unsigned int flag, vector<Function *> * list_g_taylor) {
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "f : " << *f << "\nx = " << x << endl;
    cout.flush();
    cout << "arbre courant " << endl;
    cout << _tree->to_string() << endl;
    cout << "\tdensity : " << get_density() << endl;
    cout << "\tsymetry : " << get_symetry() << endl;
    cout << "\torder : " << get_order() << endl;
    cout << "\talpha : " << *get_alpha() << endl;
#endif
    if (_delta){
      if (flag == _flag) { //we have already computed _delta
#ifdef __PRINT_DEBUG__
        cout << "delta in memory = " << *_delta << endl;flush(cout);
#endif
// 	cout << "delta\n" << endl;
        return _delta;
      }
    } else {
      _delta  = new Affine2Vector(x.size());
    }
    if (_list_of_edge->size() == 0) { //We're on a leaf of the DAG
        *_delta = f->eval_affine2_vector(x);//delta is f(x)
        _flag = flag;//memorize the fact that _delta has been computed

#ifdef __PRINT_DEBUG__
        cout << "delta is a leaf of the dag = " << *_delta << endl;flush(cout);
#endif
        return _delta;
    }

    //_delta is computed using the children of "this".
    _delta->clear();
    Affine2Matrix S(f->image_dim(), _list_of_edge->size());//S will contain the deltas of the children of "this".
    list<pair<vertex_wDAG_aff*, int> >::iterator it;
    Vector i(_list_of_edge->size());
    Vector j = Vector::zeros(_list_of_edge->size());

    int l = 0;
    for (it = _list_of_edge->begin(); it!= _list_of_edge->end(); it++) {
        pair<vertex_wDAG_aff*, int> current_child = (*it);
        Affine2Vector delta_child = *(current_child).first->get_delta(f, x, flag, list_g_taylor);
        S.set_col(l, delta_child);
        i[l] = current_child.second;
        l++;
    }
//     int dim = factorial((int)sum(i));
    int dim = ((int)sum(i));
//         int dim = 10;
//     cout << "10x" << endl;
//     cout << "1x" << endl;
//     if (dim < factorial((int)sum(i)))
//         cout << factorial((int)sum(i)) << endl;
    //TODO lower dim to the relevant directions. Should be done through get_partial_derivative by generating directly the set of "good" directions j.

    get_partial_derivative(f, x, S, i, j, dim, 0, dim, _delta, list_g_taylor);
    _flag = flag;
#ifdef __PRINT_DEBUG__
    cout << "delta from partial derivatives = " << *_delta << endl;flush(cout);
#endif
    return _delta;
}


}//namespace ibex
