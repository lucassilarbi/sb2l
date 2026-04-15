/* ============================================================================
 * D Y N I B E X - Definition of the discrete diff structure
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : March 2020
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_DIS_H
#define IBEX_DIS_H

namespace ibex{

class ivp_discrete
{
public:
//ivp discrete y(k+1)=f(y(k))
  ivp_discrete(const Function _ydot, double _t0, const IntervalVector _yinit){
      nbvar = _yinit.size();
      yinit = new IntervalVector(_yinit);
      yinit_aff = new Affine2Vector(_yinit,true);

      ydot = new Function(_ydot);
      t0 = _t0;
  };

  IntervalMatrix eval_jacobian_init()
  {
    return ydot->jacobian(*yinit);

  };

  ~ivp_discrete(){
    delete yinit;
    delete yinit_aff;
    delete ydot;
  };


public:
  unsigned int nbvar;
  Function* ydot;
  IntervalVector* yinit;
  Affine2Vector* yinit_aff;
  double t0;


};


}


#endif
