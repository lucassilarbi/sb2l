/* ============================================================================
 * D Y N I B E X - Definition of the Lobatto III C Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Author(s) of modification   : Julien Alexandre dit Sandretto and Olivier Mullier
 * Modified    : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_LC3_H
#define IBEX_SOL_LC3_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"

namespace ibex{

class solution_j_lc3 : public solution_j
{
  public:

	//method to define

	IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
	{
	  return picard_tayl(y0,_ode,2);
	}

	//the LTE
	Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode, double h)
	{
	  if (_ode->get_lte_computation() == AUTODIF){
	    if (!_ode->get_dags()->is_filled_in()){
	      ibex::IntervalMatrix A(3, 3, 0.0);
	      A[0][0] = (1.0)/(6.);
	      A[0][1] = (-1.0)/(3.);
	      A[0][2] = (1.0)/(6.);

	      A[1][0] = (1.0)/(6.);
	      A[1][1] = (5.0)/(12.);
	      A[1][2] = (-1.0)/(12.);

	      A[2][0] = (1.0)/(6.);
	      A[2][1] = (2.0)/(3.);
	      A[2][2] = (1.0)/(6.);

	      ibex::Vector b(3);
	      b[0] = (1.0)/(6.);
	      b[1] = (2.0)/(3.);
	      b[2] = (1.0)/(6.);

	      //cout << A << endl << b << endl;
	      _ode->set_Butcher_table(A, b);
	    }
	    Affine2Vector err_aff = _ode->get_dags()->get_lte(5, Affine2Vector(y0, true));
	    err_aff*=(pow(Interval(h),5) / 120.0);
// 	    cout << "y0 : " << y0 << "; LTE = " << err_aff.itv() << endl;cout.flush();
	    return err_aff;
	  } else {
	    Affine2Vector err_aff = _ode->computeLC3derivative(Affine2Vector(y0,true));
	    err_aff*=(pow(Interval(h),5) / 120);
	    return err_aff;
	  }
// 	  cout << "LTE = " << err_aff.itv() << endl << endl;
	}

	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.2)));
	}


	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){
	  //with LC3 and affine form
	  *box_jnh_aff = remainder_lc3(_ode);
	  return 1;
	}


	//constructor
	solution_j_lc3(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{

	}


	//destructor
	~solution_j_lc3(){

	}



private:

	  //lc3 with remainder
	  Affine2Vector remainder_lc3(ivp_ode* _ode)
	  {
	    double h=time_j.diam();
	    double tol = 1e-20;//atol*0.001;
// 	    cout << "box_j1" << *box_j1 << endl;
	    IntervalVector k1 = _ode->compute_derivatives(1,*box_j1);
	    IntervalVector k2(k1);
	    IntervalVector k3(k1);

	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);
	    IntervalVector k3_old(k3);

///we need a first computation to be consistent on k1, why ?? because its sum of coeff = 0
	    k1 = _ode->compute_derivatives(1,*box_jn + h*((Interval(1)/6.0)*k1 + (-Interval(1)/3.0)*k2 + (Interval(1)/6.0)*k3));
	    int i = 1;
	    do
	    {

	      k1_old=k1;
	      k2_old=k2;
	      k3_old=k3;
	      IntervalVector k1_temp(_ode->compute_derivatives(1,*box_jn + h*((Interval(1)/6.0)*k1 + (-Interval(1)/3.0)*k2 + (Interval(1.)/6.0)*k3)));
	      IntervalVector k2_temp(_ode->compute_derivatives(1,*box_jn + h*((Interval(1)/6.0)*k1 + (Interval(5)/12.0)*k2 + (-Interval(1.)/12.0)*k3)));
	      IntervalVector k3_temp(_ode->compute_derivatives(1,*box_jn + h*((Interval(1)/6.0)*k1 + (Interval(2)/3.0)*k2 + (Interval(1.)/6.0)*k3)));

	      k1 &= k1_temp;
	      k2 &= k2_temp;
	      k3 &= k3_temp;

	      i++;

	    }while((k1.rel_distance(k1_old) > tol)||(k2.rel_distance(k2_old) > tol)||(k3.rel_distance(k3_old) > tol));



	    Affine2Vector k1_aff(k1,true);
	    Affine2Vector k2_aff(k2,true);
	    Affine2Vector k3_aff(k3,true);

	    k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((Interval(1)/6.0)*k1_aff + (-Interval(1.0)/3.0)*k2_aff + (Interval(1.0)/6.0)*k3_aff));

	    k2_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((Interval(1)/6.0)*k1_aff + (Interval(5.0)/12.0)*k2_aff + (-Interval(1.0)/12.0)*k3_aff));

	    k3_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((Interval(1)/6.0)*k1_aff + (Interval(2.0)/3.0)*k2_aff + (Interval(1.0)/6.0)*k3_aff));

	    Affine2Vector lc3=*box_jn_aff + h*( (Interval(1)/6.0)*k1_aff + (Interval(2)/3.0)*k2_aff + (Interval(1)/6.0)*k3_aff);


	    return lc3+*box_err_aff;

	  };



};
}

#endif
