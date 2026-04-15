/* ============================================================================
 * D Y N I B E X - tools for geometry (zonotopes)
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto
 * Created     : April 20th, 2020
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_GEO_TOOLS_H
#define IBEX_GEO_TOOLS_H

#include <iomanip>
#include <stdlib.h>

namespace ibex{

static void one_generator(int max, int j, Affine2Vector vx, std::list<std::vector<double> >* lv)
{
	//std::cout << "passage recurrence : " << j << std::endl;
	//std::cout << "nombre courrant de vertices : " << lv.size() << std::endl;
	
	int n=vx.size();
	if (j==0)
	{
		std::vector<double> v(n);
		for (int k=0;k<n;k++)
			v[k]=vx[k].val(j);
		lv->push_back(v);
		one_generator(max,1,vx,lv);
		return;
	}
		
	if (j>max)
	{
		/*std::list<std::vector<double> >::const_iterator it;
		for(it=lv->begin();it!=lv->end();it++) 
		{
			std::cout << *it << std::endl;
		}*/
		return;
	}
	
	int test_alpha=false;
	for (int k=0;k<n;k++)
	{
		if (vx[k].val(j) != 0.0)
			test_alpha=true;
	}
	if (!test_alpha)
	{
		one_generator(max,j+1,vx,lv);
		return;
	}
		
	std::list<std::vector<double> >::iterator it;
	for(it=lv->begin();it!=lv->end();it++) 
	{
		std::vector<double> v(*it);
		for (int k=0;k<n;k++)
		{
			double alphaj=vx[k].val(j);
			v[k] -= alphaj;
			it->at(k) += alphaj;
		}
		lv->push_front(v);
		
	}
	
	one_generator(max,j+1,vx,lv);
	return;


}

static std::list<std::vector<double> > toVertices(Affine2Vector vx)
{
	std::cout << "In toVertices()" << std::endl;
	std::list<std::vector<double> > lv;
	
	//compact
	vx.compact();
	
	//suppression garbage
	vx.empty_garbage();

	
	//affichage pour test
	/*for (int j=0; j< vx.size(); j++)
	{
		Affine2Main<AF_fAFFullI> x=vx[j];
		//std::cout << x << std::endl;
		std::cout << x.val(0);
		for (int i = 1; i <= x.size(); i++) {
			double v = x.val(i);
			if (v!=0)
			{
				std::cout << std::setprecision(15) <<" + " << v << " eps_" << i;
			}
		}
		std::cout << " + " << x.err() << "[-1,1]";
		std::cout << std::endl;
	}*/
	
	
	//calcul indice max
	int max_ind = 0;
	for (int j=0; j< vx.size(); j++)
	{
		Affine2Main<AF_fAFFullI> x=vx[j];
		max_ind=std::max(max_ind, x.size());
	}
	
	//generation vertices
	one_generator(max_ind, 0, vx, &lv);	

	std::cout << "Number of vertices in list : " << lv.size() << std::endl;
	
	return lv;
}
}

#endif
