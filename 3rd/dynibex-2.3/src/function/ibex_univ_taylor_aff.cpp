/* ============================================================================
 * D Y N I B E X - Definition of the univariated taylor component computation using affine forms
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#include "ibex_univ_taylor_aff.h"

namespace ibex{



void Univ_Taylor_Aff::compute_taylor_aff(const Function& f, std::vector<Affine2Domain>& taylor_comps_aff){
#ifdef __PRINT_DEBUG__
        std::cout << "forward" << std::endl;
        cout.flush();
#endif


	f.forward<Univ_Taylor_Aff>(*this);
	int i;
	for (i = 0; i < f.expr().deco.taylor_comps_aff->size(); i++) {
// 	  if (i > 0)
// 	      cout << "i : " << i << endl;
            taylor_comps_aff.push_back(f.expr().deco.taylor_comps_aff->operator[](i));
        }
        f.backward<Univ_Taylor_Aff>(*this);

}


void Univ_Taylor_Aff::print_taylor_comps(std::vector< Affine2Domain >* tayl)
{
    if (tayl) {
        unsigned int i;
        cout << "(";
        if (tayl->size()>0) {
            cout << tayl->operator[](0);
            for (i = 1; i < tayl->size(); i++) {
                cout << ", " << tayl->operator[](i);
            }
        }
        cout << ")";
    } else {
        cout << "liste vide" << endl;

    }
    cout.flush();
}

}
