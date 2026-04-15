/* ============================================================================
 * D Y N I B E X - Definition of the Gauss order 6 Method
 * ============================================================================
 * Copyright   : ENSTA PT
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : JADS and The Chap
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique,
 * THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_GAUSS6_H
#define IBEX_SOL_GAUSS6_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{

class solution_j_gauss6 : public solution_j
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
	    Interval sqrt15 = ibex::sqrt(Interval(15.0));
	    Interval five36 = Interval(5.0) / 36.0;

	    Interval two9 = Interval(2.0) / 9.0;
	    Interval four9 = Interval (4.0) / 9.0;
	    Interval five18 = Interval (5.0) / 18.0;
	    Interval half = Interval(0.5);

	    ibex::IntervalMatrix A(3, 3, 0.);
	    A[0][0] = five36;
	    A[0][1] = two9 - (sqrt15 / 15.0);
	    A[0][2] = five36  - (sqrt15 / 30.0);

	    A[1][0] = five36 + (sqrt15 / 24.0);
	    A[1][1] = two9;
	    A[1][2] = five36 - (sqrt15 / 24.0);

	    A[2][0] = five36 + (sqrt15 / 30.0);
	    A[2][1] = two9 + (sqrt15 / 15.0);
	    A[2][2] =  five36;

	    ibex::IntervalVector b(3);
	    b[0] = five18;
	    b[1] = four9;
	    b[2] = five18;

	    _ode->set_Butcher_table(A, b);
            }

            Affine2Vector err_aff = _ode->get_dags()->get_lte(7, Affine2Vector(y0, true));
            err_aff*=(pow(Interval(h),7) / 5040.0);
            return err_aff;
        }

	  Affine2Vector err_aff = _ode->computeGL6derivative_aff(Affine2Vector(y0,true));
	    err_aff*=(pow(Interval(h),7) / 5040.0);

	    return err_aff;
	}

	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,1.0/7.0)));
	}


	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){
	  //with Gauss6 and affine form
	  *box_jnh_aff = remainder_gauss6(_ode);
	  return 1;
	};


	//constructor
	solution_j_gauss6(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{

	}


	//destructor
	~solution_j_gauss6(){
	}



private:

	  //gauss6 with remainder
	  Affine2Vector remainder_gauss6(ivp_ode* _ode)
	  {
	    double h=time_j.diam();
	    double tol = 1e-20;//atol*0.001;

	    IntervalVector k1 = _ode->compute_derivatives(1,*box_j1);
	    IntervalVector k2(k1);
	    IntervalVector k3(k1);

	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);
	    IntervalVector k3_old(k3);

	    Interval sqrt15 = ibex::sqrt(Interval(15.0));
	    Interval five36 = Interval(5.0) / 36.0;

	    Interval two9 = Interval(2.0) / 9.0;
	    Interval four9 = Interval (4.0) / 9.0;
	    Interval five18 = Interval (5.0) / 18.0;
	    Interval half = Interval(0.5);

	    Interval a11= five36;
	    Interval a12 = two9 - (sqrt15 / 15.0);
	    Interval a13 = five36  - (sqrt15 / 30.0);
	    Interval a21= five36 + (sqrt15 / 24.0);
	    Interval a22 = two9;
	    Interval a23= five36 - (sqrt15 / 24.0);
	    Interval a31= five36 + (sqrt15 / 30.0);
	    Interval a32= two9 + (sqrt15 / 15.0);
	    Interval a33=  five36;


	    Interval b1= five18;
	    Interval b2= four9;
	    Interval b3= five18;

	    do
	    {
	      k1_old=k1;
	      k2_old=k2;
	      k3_old=k3;

	      k1 &= _ode->compute_derivatives(1,*box_jn + h*(a11*k1 + a12*k2 + a13*k3));
	      k2 &= _ode->compute_derivatives(1,*box_jn + h*(a21*k1 + a22*k2 + a23*k3));
	      k3 &= _ode->compute_derivatives(1,*box_jn + h*(a31*k1 + a32*k2 + a33*k3));


	    }while ((k1.rel_distance(k1_old) > atol)||(k2.rel_distance(k2_old) > atol)||(k3.rel_distance(k3_old) > atol));

	    Affine2Vector k1_aff(k1,true);
	    Affine2Vector k2_aff(k2,true);
	    Affine2Vector k3_aff(k3,true);

	    k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a11*k1_aff + a12*k2_aff + a13*k3_aff));
	    k2_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a21*k1_aff + a22*k2_aff + a23*k3_aff));
	    k3_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a31*k1_aff + a32*k2_aff + a33*k3_aff));

	    Affine2Vector gauss6=*box_jn_aff + h*(b1*k1_aff + b2*k2_aff + b3*k3_aff);

	    return gauss6+*box_err_aff;

	  };



};
}

#endif
