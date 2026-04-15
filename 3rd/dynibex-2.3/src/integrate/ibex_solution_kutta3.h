/* ============================================================================
 * D Y N I B E X - Definition of the KUTTA3 A Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot and Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_KUTTA3_H
#define IBEX_SOL_KUTTA3_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{

class solution_j_kutta3 : public solution_j
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
                ibex::IntervalMatrix A(3, 3, 0.);
                A[1][0] = Interval(0.5);
                A[2][0] = Interval(-1.0);
                A[2][1] = Interval(2.0);
                ibex::IntervalVector b(3);
                b[0] = Interval(1.)/6;
                b[1] = Interval(2.)/3;
                b[2] = Interval(1.)/6;
                _ode->set_Butcher_table(A, b);
            }

            Affine2Vector err_aff = _ode->get_dags()->get_lte(4, Affine2Vector(y0, true));
            err_aff*=(pow(Interval(h),4) / 24.0);
            return err_aff;
        }
	    //computation of the erk33 error
	    Affine2Vector err_aff = _ode->computeKUTTA3derivative_aff(Affine2Vector(y0,true));
	    err_aff*=(pow(Interval(h),4) / 24.0);

	    return err_aff;
	}

	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.95*std::pow(1.0/test,0.25)));
	}


	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){
	  //with erk33 and affine form
	  *box_jnh_aff = remainder_kutta3(_ode);
	  return 1;
	};


	//constructor
	solution_j_kutta3(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{

	}


	//destructor
	~solution_j_kutta3(){

	}



private:

	  //kutta3 with remainder
	  Affine2Vector remainder_kutta3(ivp_ode* _ode)
	  {
	    double h=time_j.diam();
	    Affine2Vector boxj1(*box_jn_aff);

	    Affine2Vector k1 = _ode->compute_derivatives_aff(1, boxj1);

	    Affine2Vector boxj2(k1);
	    boxj2*=(Interval(0.5)*h);
	    Affine2Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj2);

	    Affine2Vector boxj3(-k1+2*k2);
	    boxj3*=h;

	    Affine2Vector k3 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj3);
	    Affine2Vector kutta3 = Interval(1.0)/6.0*k1 + Interval(2.0)/3.0*k2 + Interval(1.0)/6.0*k3;
	    kutta3*=h;
	    kutta3+=*box_jn_aff;



	    return kutta3+*box_err_aff;

	  };



};
}

#endif
