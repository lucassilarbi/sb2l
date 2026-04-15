/* ============================================================================
 * D Y N I B E X - Definition of the Solution of one simulation step for DAE
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Oct, 2015
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique,
 * THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_DAE_H
#define IBEX_SOL_DAE_H

#include <iomanip>
#include <stdlib.h>


namespace ibex{

class solution_j_dae
{
  public:
	IntervalVector* box_j0; //encadrement grossier sur t=[tn,tn+h]
	IntervalVector* box_j1; //encadrement fin sur t=[tn,tn+h]
	Interval time_j; //t=[tn,tn+h]

	IntervalVector* box_jn; //encadrement à tn
	IntervalVector* box_jnh; //encadrement à tn+h
	IntervalVector* box_xnh;

	Affine2Vector* box_err_aff;
	//at t=tn
	Affine2Vector* box_jn_aff;
	Affine2Vector* box_xn_aff;
	//at t=tn+h
	Affine2Vector* box_jnh_aff;
	Affine2Vector* box_xnh_aff;


	//accuracy expected
	double atol;

	///control of stepsize
	double factor; //factor for the next stepsize : if > 1 the current step was easy with a small error
		       // if < 1 the current stepsize used a small stepsize or the error was too large

	//few different picard operator
	//picard taylor for dae
	IntervalVector picard_tayl(IntervalVector s0, ivp_dae_h1* _dae, int ordre)
	{
	  double h = time_j.diam();
	  int n=ordre;
	  int fac_i=1;

	  Affine2Vector sn(s0,true);

	  sn.put(0,*box_jn_aff);
	  Affine2Vector int_tayl(*box_jn_aff);

	  for (int i=1;i<n;i++)
	  {
	      fac_i = fac_i*i;
	      Affine2Vector df = _dae->compute_derivatives_aff(i, sn);
	      df*=( 1.0 /fac_i);
	      df*=Interval(0,std::pow(h,i));
	      int_tayl=int_tayl + df;

	  }

	  IntervalVector err = _dae->compute_derivatives_aff(n, Affine2Vector(s0)).itv();

	  fac_i = fac_i*(n);
	  err*=(1.0/fac_i);
	  err*=Interval(0,std::pow(h,n));

	  return int_tayl.itv()+err;
	};



	//******virtual methods to define to build a new solution scheme********///

	//the picard operator
	virtual IntervalVector picard(IntervalVector y0, ivp_dae_h1* _dae, int ordre)
	{
	  return picard_tayl(y0,_dae,ordre);
	}

	//the LTE
	virtual Affine2Vector LTE(IntervalVector y0, ivp_dae_h1* _dae,double h)
	{

	    Affine2Vector err_aff = _dae->computeRADAU3derivative(Affine2Vector(y0,true));
	    err_aff*=(std::pow(h,4) / 24.0);

	    return err_aff;
	}

	//the factor for the next stepsize computation
	virtual double factor_stepsize(double test)
	{
	  return std::min(5.0,std::max(0.4,0.7*std::pow(1.0/test,0.25)));
	}

	//compute the sharpest jn+1
	virtual int calcul_jnh(ivp_dae_h1* _dae){    ///TODO not virtual, just jnh and jnh with param to virtual and not all the fixpoint

	  //with RADAU3 and affine form
	   IntervalVector sh= IntervalVector(*box_j1);

	  IntervalVector sh_old(sh);
	  IntervalVector sh_tilde(sh);
	  int count=0;

	  *box_jnh_aff = remainder_radau3(_dae,sh);

	  do{
	      sh_old = sh;
	      //with radau3 and affine form
	      *box_jnh_aff = remainder_radau3_param(_dae,sh_tilde,sh);

	      IntervalVector sh_temp(sh);
	      sh_temp.put(0,box_jnh_aff->itv());
	      sh_temp &= sh;
	      if (_dae->embedded_ctc != NULL)
	      {
		  _dae->embedded_ctc->contract(sh_temp);
		  if (sh_temp.is_empty())
		    sh_temp = sh;
	      }
	      if (sh_temp.is_empty())
		    sh_temp = sh;

	      sh &= sh_temp;

	      IntervalVector sh_x = sh.subvector(_dae->nbvar,sh.size()-1);

	      *box_xnh_aff = contract_algebraic(_dae,sh);

	      sh &= sh_old;
	      sh_x &= box_xnh_aff->itv();
	      sh.put(_dae->nbvar,sh_x);

	      count++;
	  }while ((sh.rel_distance(sh_old) > 1e-24)&&(count<30));

	  return 1;
	};



	///***algorithm for scheme, general for all scheme***///
	//verification of existence and uniqueness with Lohner test
	int uniqueness(ivp_dae_h1* _dae){


	  //Lohner test for ydot


	  //Consistency test for g
	  IntervalVector y0 = box_jn_aff->itv();
	  IntervalVector x0 = box_xn_aff->itv();
	  std::cout << "Initial values : " << y0 << ", " << x0 << std::endl;
	  IntervalVector s0 = IntervalVector(x0.size()+y0.size(),0);
	  s0.put(0,y0);
	  s0.put(_dae->nbvar,x0);

	  //IntervalVector x1 = _dae->compute_algebraic_hansen(s0);
	  IntervalVector x1 = _dae->compute_algebraic_krawczyk(s0);
	  if (x1.is_empty())
	  {
	    std::cout << "Not consistence of initial values !" << std::endl;
	    return 0;
	  }
	  else
	  {
	      if (x1.is_subset(x0))
	      {
		*box_xn_aff = Affine2Vector(x1,true);
		*box_xnh_aff = Affine2Vector(x1,true);
		std::cout << "Consistence of initial values : " << x1 << std::endl;
		return 1;
	      }
	      else
	      {
		std::cout << "Not existence proof of initial values !" << std::endl;
		return 0;
	      }
	  }
	 return 1;
	};

	//compute j0 roughly

	//picard-krawcyk checks and valides a step
	int picard_krawczyk(IntervalVector y0,IntervalVector x0, ivp_dae_h1* _dae)
	{

	  IntervalVector y_old(y0);
	  IntervalVector x_old(x0);


	  IntervalVector s0 = IntervalVector(x0.size()+y0.size(),0);

	  s0.put(0,y0);
	  s0.put(_dae->nbvar,x0);

	  double infl = std::max(0.1,s0.diam().max()*0.1);
	  s0.inflate(infl);

	  IntervalVector s1 = IntervalVector(s0);

	  //lindelof-krawczyk
	  y0 = picard_tayl(s0,_dae,3);
	  s0.put(0,y0);

	  x0 = _dae->compute_algebraic_krawczyk(s0);


	  s1.put(0,y0);
	  s1.put(_dae->nbvar,x0);

	  if (s1.is_empty())
	    s1 = s0;


	  int iter = 0;
	  while ((!s1.is_strict_subset(s0))&&(iter <= 2*_dae->nbvar))
	  {

	    s0=s1;
	    iter++;

	    double infl = s0.diam().max()*0.1;

	    //lindelof
	    IntervalVector y1 = picard_tayl(s0,_dae,2);

	    if (!y1.is_subset(y0))
	      y1.inflate(infl);

	    //krawczyk
	    IntervalVector x1 = _dae->compute_algebraic_krawczyk(s0);
	    //IntervalVector x1 = _dae->compute_algebraic_hansen(s0);

	    if (!x1.is_subset(x0))
	      x1.inflate(infl);


	    s1.put(0,y1);
	    s1.put(_dae->nbvar,x1);
	  }


	  *box_j0 = s1;
	  double h;
	  if (s1.is_strict_subset(s0))
	  {


	    //contractor
	    do{
	      s0=s1;
	      IntervalVector y1 = picard_tayl(s0,_dae,2);
	      s1.put(0,y1);
	      IntervalVector x1 = _dae->compute_algebraic_krawczyk(s1);
	      //IntervalVector x1 = _dae->compute_algebraic_hansen(s1);
	      s1.put(_dae->nbvar,x1);

	      if (_dae->embedded_ctc != NULL)
	      {
		  _dae->embedded_ctc->contract(s1);
		  if (s1.is_empty())
		    s1=s0;
	      }
	      if (s1.is_empty())
		    s1=s0;


	      s1&=s0;
	    }while (s1.rel_distance(s0) > 1e-24);



	    if (s1.is_empty())
	    {
	      *box_j0 = s0;
	      LOGGER->inc_rejected_picard();
	      h = std::max(dynibex_hmin,time_j.diam()/2.0);
	      time_j = Interval(time_j.lb(),time_j.lb()+h);
	      return 0;
	    }

	    *box_j1 = s1;
	    Affine2Vector err_aff = LTE(*box_j1,_dae, time_j.diam());
	    IntervalVector err = err_aff.itv();

	    double norm_err = infinite_norm(err);

	    //test if truncature error lower than tolerance (or if we have no choice due to the min step reaching)
	    //from Hairer
	    double rtol = atol;

	    IntervalVector y1 = s1.subvector(0,_dae->nbvar-1);
	    IntervalVector x1 = s1.subvector(_dae->nbvar,s1.size()-1);

	    double err_test = atol + infinite_norm(y1)*rtol;

	    double test = norm_err/err_test;

	    //double err_x = x1.diam().max()/x_old.diam().min();
	    //double err_y = y1.diam().max()/y_old.diam().min();

	    if ( (test <= 1.0) || (time_j.diam()<=dynibex_hmin))
	    {
	      LOGGER->inc_accepted_picard();
	      LOGGER->update_lte_max(norm_err);
	      LOGGER->update_step(time_j.diam());
	      //computation of the factor for the next step
	      factor = factor_stepsize(test);
	      *box_err_aff = err_aff;
	      return 1;
	    }
	    else
	    {
	      LOGGER->inc_rejected_picard();
	      //truncature error not accepted => step division
	      h = std::max(dynibex_hmin,time_j.diam()/2.0);
	      time_j = Interval(time_j.lb(),time_j.lb()+h);
	      return 0;
	    }
	  }
	  else
	  {
	    LOGGER->inc_rejected_picard();
	    //picard rejected => step division
	    h = std::max(1e-7,time_j.diam()/2.0);
	    time_j = Interval(time_j.lb(),time_j.lb()+h);
	    return 0;
	  }


	}

	//TODO go back to j0 and j1 method
	//compute a sharp j1
	/*int calcul_j1(ivp_dae_h1* _dae){

	  double h = time_j.diam();
	  //here : picard accepted on box_j0 => integration rule is now contracting
	  IntervalVector yn_p1(*box_j0);
	  IntervalVector yn_p0(*box_j0);

	  do
	  {
	    yn_p0=yn_p1;

	    //evaluation with affine
	    yn_p1 &= picard(yn_p0,_ode,3);
	    if (_ode->embedded_ctc != NULL)
	    {
	      _ode->embedded_ctc->contract(yn_p1);

	    }

	  } while(yn_p1.rel_distance(yn_p0) > 1e-18); //till a fix point

	  *box_j1 = yn_p1;


	  return 1;
	};*/



	double infinite_norm(const IntervalVector _vec)
	{
	  double res = std::max(std::abs(_vec[0].lb()),std::abs(_vec[0].ub()));
	  for (int i=1;i<_vec.size();i++)
	  {
	   res =  std::max(std::abs(_vec[i].lb()),res);
	   res =  std::max(std::abs(_vec[i].ub()),res);
	  }
	  return res;
	};


	void set_atol(double a)
	{
	  atol = a;
	}

	//flush affine form to keep memory
	void flush()
	{
	  *box_err_aff = box_err_aff->itv();
	  *box_jn_aff = *box_jn;
	  *box_jnh_aff = *box_jnh;
	  *box_xn_aff = box_xn_aff->itv();
	  *box_xnh_aff = box_xnh_aff->itv();
	};


	//empty constructor
	solution_j_dae(){};

	solution_j_dae(const Affine2Vector _box_jn, const Affine2Vector _box_xn, double tn, double h, double a)
	{
	  box_jn_aff = new Affine2Vector(_box_jn);
	  box_xn_aff = new Affine2Vector(_box_xn);

	  //important:
	  time_j = Interval(tn-h,tn);
	  box_jnh_aff = new Affine2Vector(_box_jn);
	  box_xnh_aff = new Affine2Vector(_box_xn);

	  box_j0 = new IntervalVector(_box_jn.size()+_box_xn.size());
	  box_j1 = new IntervalVector(*box_j0);

	  box_err_aff = new Affine2Vector(_box_jn.size());

	  box_jn = new IntervalVector(_box_jn.itv());
	  box_jnh  = new IntervalVector(_box_jn);
	  box_xnh  = new IntervalVector(_box_xn);
	  atol = a;
	  factor = 1.0;
	}

	//constructor
	solution_j_dae(const Affine2Vector _box_jn, const Affine2Vector _box_xn, double tn, double h, ivp_dae_h1* _dae,double a, double fac)
	{
	  factor = fac;
  	  atol = a;
	  time_j = Interval(tn,tn+std::min(dynibex_hmax,std::max(dynibex_hmin,h*fac)));

	  box_jn_aff = new Affine2Vector(_box_jn);
	  box_xn_aff = new Affine2Vector(_box_xn);

	  box_jnh_aff = new Affine2Vector(_box_jn.size());
	  box_xnh_aff = new Affine2Vector(_box_xn.size());

	  box_j0 = new IntervalVector(_box_jn.size()+_box_xn.size());
	  box_j1 = new IntervalVector(*box_j0);

	  box_err_aff = new Affine2Vector(_box_jn.size());


	  box_jn = new IntervalVector(_box_jn.itv());
	  box_jnh  = new IntervalVector(_box_jn);
	  box_xnh  = new IntervalVector(_box_xn);

	}

	int compute_oneStep(const Affine2Vector _box_jn,const Affine2Vector _box_xn, ivp_dae_h1* dae)
	{

	  //guard for degenerated algebraic
	  IntervalVector x0 = _box_xn.itv();
	  if (x0.is_flat())
	    x0.inflate(atol);

	  //operator of picard-krawcyk
	  int ok = picard_krawczyk(_box_jn.itv(),x0, dae);
	  int nb=0;
	  while ((ok != 1)&&(nb < 20*dae->nbvar))
	  {
	       IntervalVector s1(*box_j0);
	       IntervalVector y1 = s1.subvector(0,dae->nbvar-1);
	       IntervalVector x1 = s1.subvector(dae->nbvar,s1.size()-1);
	       ok=picard_krawczyk(y1,x1, dae);
	       nb++;
	  }
	  if (ok != 1)
	  {
		std::cout << "Step:"<<nb<<" refused !" << std::endl;
		LOGGER->Log("Solution at t=%f : ",time_j.lb());
		LOGGER->Log_sol(_box_jn.itv());
		LOGGER->Log_end();
		exit(EXIT_FAILURE);

	  }
	  else
	  {
		if (box_j0->is_unbounded())
		{
		  std::cout << "Step:"<<nb<<" refused !" << std::endl;
		  LOGGER->Log("Solution at t=%f : ",time_j.lb());
		  LOGGER->Log_sol(_box_jn.itv());
		  LOGGER->Log_end();
		  exit(EXIT_FAILURE);
		}
	  }

	  //compute jnh
	  calcul_jnh(dae);
	  box_jnh_aff->compact();   ///compactage
	  *box_jnh = box_jnh_aff->itv();

	  box_xnh_aff->compact();   ///compactage
	  *box_xnh = box_xnh_aff->itv();

	  ///test for dae2 with exact solution
	  /*Interval x(time_j.ub());
	  cout << "sol exacte : " << endl;
	  cout << "t=" << time_j.ub() << endl;

	  Interval y0e= sin(x)+5*cos((x*x)/2.0);
	  Interval y1e= cos(x)+5*sin((x*x)/2.0);
	  Interval x0e= -cos(x);
	  Interval x1e= sin(x);

	  cout << "y[0] : " << y0e<< endl;
	  cout << "y[1] : " <<  y1e<< endl;
	  cout << "x[0] : " << x0e << endl;
	  cout << "x[1] : " << x1e << endl;


	  if (!(y0e.is_subset((*box_jnh)[0]))||
	    !(y1e.is_subset((*box_jnh)[1]))||
	    !(x0e.is_subset((*box_xnh)[0])) ||
	    !(x1e.is_subset((*box_xnh)[1])))
	  {
	    std::cout << "non inclu ! " << std::endl;
	  }*/

	  return 1;
	}


	//printer
	void print_soljn(){
	    std::cout << std::setprecision(20)<< "Solution at t=" << time_j.ub() << " : " << *box_jnh << std::endl;
	    std::cout << "affine form : " << *box_jnh_aff << std::endl;
	}


	//destructor
	~solution_j_dae(){
	}

	void destructor()
	{
	  if (box_jn_aff != NULL)
	    delete box_jn_aff;
	  if (box_jn != NULL)
	    delete box_jn;
	  if (box_j0 != NULL)
	    delete box_j0;
	  if (box_j1  != NULL)
	    delete box_j1;
	  if (box_err_aff != NULL)
	    delete box_err_aff;
	  if (box_jnh  != NULL)
	    delete box_jnh;
	  if (box_jnh_aff != NULL)
	    delete box_jnh_aff;
	  if (box_xnh != NULL)
	    delete box_xnh;
	  if (box_xn_aff != NULL)
	    delete box_xn_aff;
	  if (box_xnh_aff != NULL)
	    delete box_xnh_aff;
	}

	Affine2Vector contract_algebraic(ivp_dae_h1* _dae, IntervalVector sh)
	  {
	    IntervalVector sh_old(sh);
	    IntervalVector xnh = sh.subvector(_dae->nbvar,sh.size()-1);

	    CtcFwdBwd c(*_dae->g);
	    CtcFixPoint fix(c);
	    fix.contract(sh);

	    if (sh.is_empty())
	      sh=sh_old;
	    xnh &= _dae->compute_algebraic_krawczyk(sh);
	    //xnh &= _dae->compute_algebraic_hansen(sh);
	    return Affine2Vector(xnh,true);

	  }


	//radau3 with remainder and parameter
	  Affine2Vector remainder_radau3_param(ivp_dae_h1* _dae, IntervalVector sh, IntervalVector sh_x)
	  {
	    double h=time_j.diam();
	    double tol = 1e-20;

	    IntervalVector k1 = _dae->compute_derivatives_aff(1, Affine2Vector(*box_j1,true)).itv();

	    IntervalVector k2(k1);


	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);
	    IntervalVector s_temp(sh);
	    do
	    {
	      k1_old=k1;
	      k2_old=k2;

	      IntervalVector s1(sh);
	      s1.put(0,box_jn_aff->itv() + h*((5.0/12.0)*k1 + (-1.0/12.0)*k2));
	      k1 &= _dae->compute_derivatives_aff(1,Affine2Vector(s1,true)).itv();


	      IntervalVector s2(sh_x);
	      s2.put(0,box_jn_aff->itv() + h*((3.0/4.0)*k1 + (1.0/4.0)*k2));
	      k2 &= _dae->compute_derivatives_aff(1,Affine2Vector(s2,true)).itv();


	    }while ((k1.rel_distance(k1_old) > tol)||(k2.rel_distance(k2_old) > tol));

	    Affine2Vector k1_aff = Affine2Vector(k1,true);
	    Affine2Vector k2_aff = Affine2Vector(k2,true);


	    Affine2Vector s1_aff(sh);
	    s1_aff.put(0,*box_jn_aff + h*((5.0/12.0)*k1_aff  + (-1.0/12.0)*k2_aff));
	    k1_aff = _dae->compute_derivatives_aff(1,s1_aff);


	    Affine2Vector s2_aff(sh_x); //(sh);
	    s2_aff.put(0,*box_jn_aff + h*((3.0/4.0)*k1_aff + (1.0/4.0)*k2_aff));
	    k2_aff = _dae->compute_derivatives_aff(1,s2_aff);


	    //2 times for the x contraction taken into account
	    s1_aff.put(0,*box_jn_aff + h*((5.0/12.0)*k1_aff  + (-1.0/12.0)*k2_aff));
	    k1_aff = _dae->compute_derivatives_aff(1,s1_aff);
	    s2_aff.put(0,*box_jn_aff + h*((3.0/4.0)*k1_aff + (1.0/4.0)*k2_aff));
	    k2_aff = _dae->compute_derivatives_aff(1,s2_aff);


	    Affine2Vector radau3=*box_jn_aff + h*( (3.0/4.0)*k1_aff + (1.0/4.0)*k2_aff);

	    return radau3+*box_err_aff;

	  };

	  //radau3 with remainder
	  Affine2Vector remainder_radau3(ivp_dae_h1* _dae, IntervalVector sh)
	  {
	    double h=time_j.diam();
	    double tol = 1e-20;

	    IntervalVector k1 = _dae->compute_derivatives_aff(1, Affine2Vector(*box_j1,true)).itv();

	    IntervalVector k2(k1);


	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);

	    do
	    {
	      k1_old=k1;
	      k2_old=k2;

	      IntervalVector s1(sh);
	      s1.put(0,box_jn_aff->itv() + h*((5.0/12.0)*k1 + (-1.0/12.0)*k2));
	      k1 &= _dae->compute_derivatives_aff(1,Affine2Vector(s1,true)).itv();

	      IntervalVector s2(sh);
	      s2.put(0,box_jn_aff->itv() + h*((3.0/4.0)*k1 + (1.0/4.0)*k2));
	      k2 &= _dae->compute_derivatives_aff(1,Affine2Vector(s2,true)).itv();


	    }while ((k1.rel_distance(k1_old) > tol)||(k2.rel_distance(k2_old) > tol));

	    Affine2Vector k1_aff = Affine2Vector(k1,true);
	    Affine2Vector k2_aff = Affine2Vector(k2,true);


	    Affine2Vector s1_aff(sh);
	    s1_aff.put(0,*box_jn_aff + h*((5.0/12.0)*k1_aff  + (-1.0/12.0)*k2_aff));
	    k1_aff = _dae->compute_derivatives_aff(1,s1_aff);


	    Affine2Vector s2_aff(sh);
	    s2_aff.put(0,*box_jn_aff + h*((3.0/4.0)*k1_aff + (1.0/4.0)*k2_aff));
	    k2_aff = _dae->compute_derivatives_aff(1,s2_aff);

	    Affine2Vector radau3=*box_jn_aff + h*( (3.0/4.0)*k1_aff + (1.0/4.0)*k2_aff);

	    return radau3+*box_err_aff;

	  };





};
}

#endif
