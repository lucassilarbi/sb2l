/* ============================================================================
 * D Y N I B E X - Definition of the Gauss Legendre 4 Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2015
 * Author(s) of modification   : Julien Alexandre dit Sandretto and Olivier Mullier
 * Modified    : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_GL4_H
#define IBEX_SOL_GL4_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{

class solution_j_gauss4 : public solution_j
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
                ibex::IntervalMatrix A(2, 2, 0.);
                A[0][0] = Interval(0.25000000000000000E0, 0.25000000000000000E0);
                A[0][1] = Interval(-0.38675134594812977E-1, -0.38675134594812810E-1);
                A[1][0] = Interval(0.53867513459481264E0, 0.53867513459481309E0);
                A[1][1] = Interval(0.25000000000000000E0, 0.25000000000000000E0);
                ibex::IntervalVector b(2);
                b[0] = Interval(.5);
                b[1] = Interval(.5);
                _ode->set_Butcher_table(A, b);
            }
            
            Affine2Vector err_aff = _ode->get_dags()->get_lte(5, Affine2Vector(y0, true));  
            err_aff*=(pow(Interval(h),5) / 120);
            return err_aff;
        }
	    Affine2Vector err_aff = _ode->computeGL4derivative_aff(Affine2Vector(y0,true));    
	    err_aff*=(pow(Interval(h),5) / 120);
	    return err_aff;	  

	}

	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.2)));
	}


	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){
	  //with GL4 and affine form
	  *box_jnh_aff = remainder_gl4(_ode);
	  return 1;
	};


	//constructor
	solution_j_gauss4(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{

	}


	//destructor
	~solution_j_gauss4(){

	}



private:

	  //gl4 with remainder
	  Affine2Vector remainder_gl4(ivp_ode* _ode)
	  {
	    double h=time_j.diam();
	    double tol = 1e-20;//atol*0.001;

	    IntervalVector k1 = _ode->compute_derivatives(1,*box_j1);
	    IntervalVector k2(k1);


	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);



	    do
	    {

	      k1_old=k1;
	      k2_old=k2;



	      k1 &= _ode->compute_derivatives(1,*box_jn + h*((Interval(0.25000000000000000E0, 0.25000000000000000E0))*k1 + (Interval(-0.38675134594812977E-1, -0.38675134594812810E-1))*k2 ));
	      k2 &= _ode->compute_derivatives(1,*box_jn + h*((Interval(0.53867513459481264E0, 0.53867513459481309E0))*k1 + (Interval(0.25000000000000000E0, 0.25000000000000000E0))*k2 ));


	    }while((k1.rel_distance(k1_old) > tol)||(k2.rel_distance(k2_old) > tol));



	    Affine2Vector k1_aff(k1,true);
	    Affine2Vector k2_aff(k2,true);


	    k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((Interval(0.25000000000000000E0, 0.25000000000000000E0))*k1_aff + (Interval(-0.38675134594812977E-1, -0.38675134594812810E-1))*k2_aff));

	    k2_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((Interval(0.53867513459481264E0, 0.53867513459481309E0))*k1_aff + (Interval(0.25000000000000000E0, 0.25000000000000000E0))*k2_aff));



	    Affine2Vector gl4=*box_jn_aff + h*( (Interval(0.50000000000000000E0, 0.50000000000000000E0))*k1_aff + (Interval(0.50000000000000000E0, 0.50000000000000000E0))*k2_aff);


	    return gl4+*box_err_aff;

	  };



};
}

#endif
