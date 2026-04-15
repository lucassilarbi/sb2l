/* ============================================================================
 * D Y N I B E X - Definition of the Radau Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot and Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_SOL_RADAU5_H
#define IBEX_SOL_RADAU5_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"

namespace ibex{

  class solution_j_radau5 : public solution_j
  {
  public:

    //method to define

    IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
    {
      // cerr << "RADAU5 avec Picard Taylor 4" << endl;
      //      return picard_tayl(y0,_ode,ordre);
      return picard_tayl(y0,_ode, ordre);
    }

    //the LTE
    Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode, double h)
    {
        if (_ode->get_lte_computation() == AUTODIF){
            if (!_ode->get_dags()->is_filled_in()){
                ibex::IntervalMatrix A(3, 3, 0.0);
                Interval racine = ibex::sqrt(Interval(6., 6.));

            //       A[0][0] = Interval(11.0)/Interval(45.0) - (Interval(7.0) * racine) / Interval(360.0);
            //       A[0][1] = Interval(37.0)/Interval(225.0) - (Interval(169.0) * racine) / Interval(1800.0);
            //       A[0][2] = -Interval(2.0)/Interval(225.0) + racine / Interval(75.0);
            //
            //       A[1][0] = Interval(37.0)/Interval(225.0) + (Interval(169.0) * racine) / Interval(1800.0);
            //       A[1][1] = Interval(11.0)/Interval(45.0) + (Interval(7.0) * racine) / Interval(360.0);
            //       A[1][2] = -Interval(2.0)/Interval(225.0) - racine / Interval(75.0);
            //
            //       A[2][0] = Interval(4.0)/Interval(9.0) - racine / Interval(36.0);
            //       A[2][1] = Interval(4.0)/Interval(9.0) + racine / Interval(36.0);
            //       A[2][2] = Interval(1.0) / Interval(9.0);
            //
            //       ibex::IntervalVector b(3);
            //       b[0] = Interval(4.0)/Interval(9.0) - racine / Interval(36.0);
            //       b[1] = Interval(4.0)/Interval(9.0) + racine / Interval(36.0);
            //       b[2] = Interval(1.0) / Interval(9.0);

            // 	A[0][0] = (11.0)/(45.0) - ((7.0) * racine) / (360.0);
            // 	A[0][1] = (37.0)/(225.0) - ((169.0) * racine) / (1800.0);
            // 	A[0][2] = -(2.0)/(225.0) + racine / (75.0);
            //
            // 	A[1][0] = (37.0)/(225.0) + ((169.0) * racine) / (1800.0);
            // 	A[1][1] = (11.0)/(45.0) + ((7.0) * racine) / (360.0);
            // 	A[1][2] = -(2.0)/(225.0) - racine / (75.0);
            //
            // 	A[2][0] = (4.0)/(9.0) - racine / (36.0);
            // 	A[2][1] = (4.0)/(9.0) + racine / (36.0);
            // 	A[2][2] = (1.0) / (9.0);
            //
            // 	ibex::IntervalVector b(3);
            // 	b[0] = (4.0)/(9.0) - racine / (36.0);
            // 	b[1] = (4.0)/(9.0) + racine / (36.0);
            // 	b[2] = (1.0) / (9.0);

                Interval a11(0.19681547722366038E0, 0.19681547722366044E0);
                Interval a12(-0.65535425850198448E-1, -0.65535425850198337E-1);
                Interval a13(0.23770974348220147E-1, 0.23770974348220158E-1);
                Interval a21(0.39442431473908718E0, 0.39442431473908740E0);
                Interval a22(0.29207341166522843E0, 0.29207341166522854E0);
                Interval a23(-0.41548752125997943E-1, -0.41548752125997922E-1);
                Interval a31(0.37640306270046725E0, 0.37640306270046731E0);
                Interval a32(0.51248582618842153E0, 0.51248582618842164E0);
                Interval a33(0.11111111111111110E0, 0.11111111111111112E0);


                Interval b1(0.37640306270046725E0, 0.37640306270046731E0);
                Interval b2(0.51248582618842153E0, 0.51248582618842164E0);
                Interval b3(0.11111111111111110E0, 0.11111111111111112E0);

                A[0][0] = a11;
                A[0][1] = a12;
                A[0][2] = a13;

                A[1][0] = a21;
                A[1][1] = a22;
                A[1][2] = a23;

                A[2][0] = a31;
                A[2][1] = a32;
                A[2][2] = a33;

                ibex::IntervalVector b(3);
                b[0] = b1;
                b[1] = b2;
                b[2] = b3;


                _ode->set_Butcher_table(A, b);
            }
            Affine2Vector y0_aff = Affine2Vector(y0, true);
            Affine2Vector err_aff = _ode->get_dags()->get_lte(6, y0_aff);
            //      cerr << "y0 = " << y0 << endl;
            //  cerr << "LTE = " << err_aff << endl;
            err_aff*=(pow(Interval(h),6) / 720.0);
            return err_aff;
        } else {
            Affine2Vector err_aff = _ode->computeRADAU5derivative(Affine2Vector(y0,true));
            err_aff*=(pow(Interval(h),6) / 720.0);
            return err_aff;
        }

    }

    //the factor for the next stepsize computation
    double factor_stepsize(double test)
    {
      return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,1.0/6.0)));
    }


    //compute the sharpest jn+1
    int calcul_jnh(ivp_ode* _ode){
      //with RADAU5 and affine form
      *box_jnh_aff = remainder_radau5(_ode);
      return 1;
    };


    //constructor
  solution_j_radau5(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
      {

      }


    //destructor
    ~solution_j_radau5(){
    }



  private:

    //radau5 with remainder
    Affine2Vector remainder_radau5(ivp_ode* _ode)
    {
      double h=time_j.diam();
      double tol = 1e-20;//atol*0.001;

      IntervalVector k1 = _ode->compute_derivatives(1,*box_j1);
      IntervalVector k2(k1);
      IntervalVector k3(k1);

      IntervalVector k1_old(k1);
      IntervalVector k2_old(k2);
      IntervalVector k3_old(k3);

      Interval a11(0.19681547722366038E0, 0.19681547722366044E0);
      Interval a12(-0.65535425850198448E-1, -0.65535425850198337E-1);
      Interval a13(0.23770974348220147E-1, 0.23770974348220158E-1);
      Interval a21(0.39442431473908718E0, 0.39442431473908740E0);
      Interval a22(0.29207341166522843E0, 0.29207341166522854E0);
      Interval a23(-0.41548752125997943E-1, -0.41548752125997922E-1);
      Interval a31(0.37640306270046725E0, 0.37640306270046731E0);
      Interval a32(0.51248582618842153E0, 0.51248582618842164E0);
      Interval a33(0.11111111111111110E0, 0.11111111111111112E0);


      Interval b1(0.37640306270046725E0, 0.37640306270046731E0);
      Interval b2(0.51248582618842153E0, 0.51248582618842164E0);
      Interval b3(0.11111111111111110E0, 0.11111111111111112E0);

      do
	{
	  k1_old=k1;
	  k2_old=k2;
	  k3_old=k3;

	  k1 &= _ode->compute_derivatives(1,*box_jn + h*(a11*k1 + a12*k2 + a13*k3));
	  k2 &= _ode->compute_derivatives(1,*box_jn + h*(a21*k1 + a22*k2 + a23*k3));
	  k3 &= _ode->compute_derivatives(1,*box_jn + h*(a31*k1 + a32*k2 + a33*k3));


	}while ((k1.rel_distance(k1_old) > tol)||(k2.rel_distance(k2_old) > tol)||(k3.rel_distance(k3_old) > tol));

      Affine2Vector k1_aff(k1,true);
      Affine2Vector k2_aff(k2,true);
      Affine2Vector k3_aff(k3,true);

      k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a11*k1_aff + a12*k2_aff + a13*k3_aff));
      k2_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a21*k1_aff + a22*k2_aff + a23*k3_aff));
      k3_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*(a31*k1_aff + a32*k2_aff + a33*k3_aff));

      Affine2Vector radau5=*box_jn_aff + h*(b1*k1_aff + b2*k2_aff + b3*k3_aff);

      return radau5+*box_err_aff;

    };



  };
}

#endif
