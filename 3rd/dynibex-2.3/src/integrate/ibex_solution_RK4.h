/* ============================================================================
 * D Y N I B E X - Definition of the Runge-Kutta order 4 Method
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
#ifndef IBEX_SOL_RK4_H
#define IBEX_SOL_RK4_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{

class solution_j_rk4 : public solution_j
{
  public:
	//method to define

	IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
	{
	  return picard_tayl(y0,_ode,ordre);
	}

	//the LTE
	Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode, double h)
	{
	  if (_ode->get_lte_computation() == AUTODIF){
	    if (!_ode->get_dags()->is_filled_in()){
	      ibex::IntervalMatrix A(4, 4, 0.0);
	      A[1][0] = Interval(0.5);
	      A[2][1] = Interval(0.5);
	      A[3][2] = Interval(1.0);
	      ibex::IntervalVector b(4);
	      b[0] = Interval(1.0)/6;
	      b[1] = Interval(1.0)/3;
	      b[2] = Interval(1.0)/3;
	      b[3] = Interval(1.0)/6;
	      _ode->set_Butcher_table(A, b);
	    }
	    Affine2Vector temp(y0, true);
	    Affine2Vector err_aff = _ode->get_dags()->get_lte(5, temp);
	    err_aff*=(pow(Interval(h),5) / 120);
	    return err_aff;
	  } else {
	    Affine2Vector err_aff = _ode->computeRK4derivative(Affine2Vector(y0,true));
	    err_aff*=(pow(Interval(h),5) / 120);
	    return err_aff;
	  }
	}

	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.95*std::pow(1.0/test,0.2)));
	}

	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){
	  //with RK4 and affine form
	  *box_jnh_aff = remainder_rk4(_ode);
	  return 1;
	};


	//constructor
	solution_j_rk4(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{


	}


	//destructor
	~solution_j_rk4(){

	}



private:

	  //rk4 with remainder
	  Affine2Vector remainder_rk4(ivp_ode* _ode)
	  {
	    double h=time_j.diam();

	    Affine2Vector boxj1(*box_jn_aff);
	    Affine2Vector k1 = _ode->compute_derivatives_aff(1, boxj1);
	    Affine2Vector boxj2(k1);
	    boxj2*=(Interval(0.5)*h);


	    Affine2Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj2);


	    Affine2Vector boxj3(k2);
	    boxj3*=(Interval(0.5)*h);


	    Affine2Vector k3 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj3);


	    Affine2Vector boxj4(k3);
	    boxj4*=(h);


	    Affine2Vector k4 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj4);

	    k2*=(2.0);
	    k3*=(2.0);


	    Affine2Vector int_rk4 = k1+k2+k3+k4;
	    int_rk4*=(Interval(h)/6.0);
	    int_rk4+=*box_jn_aff;

	    return int_rk4+*box_err_aff;//df

	  };



};
}

#endif
