/* ============================================================================
 * D Y N I B E X - Definition of the generic Solution of one simulation step
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
#ifndef IBEX_SOL_G_H
#define IBEX_SOL_G_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"
#include "ibex_solution_dae.h"


namespace ibex{

class solution_g : public solution_j
{
public:
  solution_g(Interval _time_j, IntervalVector* _box_j1, IntervalVector* _box_jn, IntervalVector* _box_jnh, // Affine2Vector* _box_jn_aff, Affine2Vector* _box_jnh_aff,
	     int nb_var)
  {
	  time_j=_time_j;
	  box_j1 = new IntervalVector(_box_j1->subvector(0,nb_var-1));
	  box_jn = new IntervalVector(*_box_jn);
	  box_jnh = new IntervalVector(*_box_jnh);
	  // box_jn_aff = new Affine2Vector (*_box_jn_aff);
	  // box_jnh_aff = new Affine2Vector (*_box_jnh_aff);
	}
  ~solution_g(){
  }

  void destructor()
  {
    if (box_jn != NULL)
      delete box_jn;
    if (box_jnh != NULL)
      delete box_jnh;
    if (box_j1  != NULL)
      delete box_j1;
  }

};

}

#endif
