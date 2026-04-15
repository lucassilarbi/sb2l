/* ============================================================================
 * D Y N I B E X - Definition of the weighted dags used for LTE
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : Jul 18, 2015
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#include "ibex_wdag_aff.h"

using namespace std;

namespace ibex{


/**
 * \brief constructor of a wDAG_aff when _A and _b are to be given (when Runge-kutta like method is applied)
 * \param f function for which computations will be made.
 * \param order order max of the created rooted trees
 * \param A matrix A in the Butcher tableau
 * \param b vector b in the Butcher tableau
 */
wDAG_aff::wDAG_aff(Function* f, int order, IntervalMatrix & A, IntervalVector & b)
{
#ifdef __PRINT_DEBUG__
  cout << "******************************************************************************************" << endl;
  cout << __func__ << endl;
  cout << "A = " << endl << A << endl;
  cout << "b = " << b << endl;
  cout << "order = " << order << endl;
#endif

  _f = f;

  _max_order = order;

  _list_of_trees = generate(order);

  _list_of_dags = new list < list<vertex_wDAG_aff *> >();

  for (int i = 0; i < order; i++) {
      list<vertex_wDAG_aff *> * current_list = new list<vertex_wDAG_aff *>();
      _list_of_dags->push_back(*current_list);
  }

  _A = new IntervalMatrix(A);
  _b = new IntervalVector(b);

  generate_dags(order, _list_of_trees, _list_of_dags, _A, _b);

  _filled_in = true;

  Variable t(1);;
  Variable x_var(_f->nb_var());
  Variable s_var(_f->image_dim());
  _list_g_taylor = new vector<Function * >();
  for (int i = 0; i < _f->image_dim(); i++){
    Function * g_i = new Function(t, x_var, s_var, (f->operator[](i))(x_var + t * s_var));
    _list_g_taylor->push_back(g_i);
  }

  _flag = 0;

#ifdef __PRINT_DEBUG__
    cout << "Nombre de dags = " << _list_of_dags->size() << endl;
    cout << "_A = " << endl << *_A << endl;
    cout << "_b = " << *_b << endl;
#endif
}

/**
 * \brief constructor of a wDAG_aff when _A and _b are not provided
 * \param f function for which computations will be made.
 * \param order order max of the created rooted trees
 */
wDAG_aff::wDAG_aff(Function* f, int order)
{
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "order = " << order << endl;
#endif

    _f = f;
    _max_order = order;
    _list_of_dags = new list < list<vertex_wDAG_aff *> >();
    _A = NULL;
    _b = NULL;

    FILE *fichier = NULL;
    fichier = fopen("dags.dat",  "r");
    bool dag_from_file = fichier != NULL;
    if (dag_from_file){
//         cout << "creation a partir du fichier" << endl;
        fclose(fichier);
        dag_from_file = file_to_dags(order, _list_of_dags, "dags.dat");
        _list_of_trees = NULL;
//         cout << "création terminée" << endl;
    }
    if (!dag_from_file) {
//       cout << "fichier non trouvé" << endl;
        _list_of_trees = generate(order);
        for (int i = 0; i < order; i++) {
            list<vertex_wDAG_aff *> * current_list = new list<vertex_wDAG_aff *>();
            _list_of_dags->push_back(*current_list);
        }
        generate_dags(order, _list_of_trees, _list_of_dags, _A, _b);
        dags_to_file("dags.dat");
  }
  _filled_in = false;

  Variable t(1);
  Variable x_var(_f->nb_var());
  Variable s_var(_f->image_dim());

  // Construction of the list of functions used for the computation of the elementary differentials.
  _list_g_taylor = new vector<Function * >();
    for (int i = 0; i < _f->image_dim(); i++){
      Function * g_i = new Function(t, x_var, s_var, (f->operator[](i))(x_var + t * s_var));
      _list_g_taylor->push_back(g_i);
  }

  _flag = 0;
//   cout << "fin du constructeur" << endl;
#ifdef __PRINT_DEBUG__
    cout << "Nombre de dags = " << _list_of_dags->size() << endl;
#endif
}

/**
 * \brief return the dag with trees of order <= \param order
 * \return the entry point of all rooted trees of order=\param order in the dag.
 */
list<vertex_wDAG_aff *> * wDAG_aff::get_dag(int order)
{
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "order = " << order << endl;
    cout << "Number of dags : " << (*this).get_list_of_dags()->size() << endl;
#endif

  if (order == _max_order){
    return &_list_of_dags->back();
  }

  list<list<vertex_wDAG_aff *> >::iterator it = _list_of_dags->begin();
  for (int i = 0; i < order-1; i++)
    it++;

#ifdef __PRINT_DEBUG__
    cout << "Nombre de points d'entrée du dag trouvé : " << (*it).size() << endl;cout.flush();
#endif
  return &(*it);
}

/**
 * \brief computes value of the derivatives of given order of y_init
 */
Affine2Vector wDAG_aff::get_derivatives(int order, const Affine2Vector & y_init){
#ifdef __PRINT_DEBUG__
  cout << "******************************************************************************************" << endl;
  cout << __func__ << endl;
  cout << "order = " << order << endl;
  cout << "y_init = " << y_init << endl;
  cout << *_f << endl;
  cout << "ordre du dag : " << get_order() << endl;
#endif

  if (order == 1)
    return _f->eval_affine2_vector(y_init);

  Affine2Vector result(_f->image_dim());
  result.clear();
  list<vertex_wDAG_aff *> * dag = get_dag(order);
  list<vertex_wDAG_aff *>::iterator it;
  unsigned int flag = ++_flag;//This flag is to assure the computation and not the use of previous results registered in the dag.
  int i = 1;

  for (it = dag->begin(); it != dag->end(); it++){
#ifdef __PRINT_DEBUG__
    cout << "DEBUT iteration sur les arbres du dag" << endl;
#endif

    Interval alpha = *(*it)->get_alpha();

#ifdef __PRINT_DEBUG__
    cout << "with tree " << __func__ << " : " << (*it)->get_tree()->to_string() << endl;
    cout << "alpha is " << __func__ << " : " << alpha << endl;
    cout << " DAG tree " << i << " " << *(*it)->get_delta(_f, y_init, flag, _list_g_taylor) << endl; flush(cout);
#endif

    if (alpha.lb() != 0.0 || alpha.ub() != 0.0){//if alpha != 0
	result += alpha * *(*it)->get_delta(_f, y_init, flag, _list_g_taylor);
    }

#ifdef __PRINT_DEBUG__
    cout << "FIN iteration sur les arbres du dag" << endl << endl;
#endif

    i++;
  }
#ifdef __PRINT_DEBUG__
  cout << "result of " << __func__ << " : " << result << endl;
#endif
  return result;
}

Affine2Vector wDAG_aff::get_lte(int order, const Affine2Vector & y_init){
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "order = " << order << endl;
    cout << "y_init = " << y_init << endl;
#endif

  Affine2Vector result(_f->image_dim());
  result.clear();
  list<vertex_wDAG_aff *> * dag = get_dag(order);
  list<vertex_wDAG_aff *>::iterator it;
  unsigned int flag = ++_flag;//This flag is to assure the computation and not the use of previous results registered in the dag.
  for (it = dag->begin(); it != dag->end(); it++){
      vertex_wDAG_aff * current_vertex = *it;
      Affine2Vector * current_delta = current_vertex->get_delta(_f, y_init, flag, _list_g_taylor);
      Interval coeff = 1.0 - 1.0*current_vertex->get_density() * (*current_vertex->get_psi(_A, _b));

#ifdef __PRINT_DEBUG__
      cout << "\talpha = " << *current_vertex->get_alpha() << endl;
      cout << "\tdensity = " << current_vertex->get_density() << endl;
      cout << "\tpsi = " << *current_vertex->get_psi(_A, _b) << endl;
      cout << "\t[1 - psi * phi] = " << coeff << endl;
      cout << "\tdelta = " << *current_delta << endl;
      cout << endl << endl;
#endif

    result+= (*current_vertex->get_alpha()) * coeff * (*current_delta);
  }

#ifdef __PRINT_DEBUG__
    cout << "result of " << __func__ << " : " << result << endl;
#endif
  return result;
}


/**
 * \brief find a vertex in the list of trees. Vertex must be present in the list.
 **/
vertex_wDAG_aff * find_vertex(list<vertex_wDAG_aff *> * list_of_vertices, rootedTree * tree) {
#ifdef __PRINT_DEBUG__
  cout << "******************************************************************************************" << endl;
  cout << __func__ << endl;
#endif

    list<vertex_wDAG_aff *>::iterator it_vertices;
    for (it_vertices = list_of_vertices->begin(); it_vertices != list_of_vertices->end(); it_vertices++) {
        vertex_wDAG_aff * current_vertex = *it_vertices;
        if (current_vertex->get_tree()->operator==(*tree))
            return current_vertex;
    }
    cerr << "Error: child has not been found but should be." << endl;
    return 0;
}

/**
 * \brief generation of the dags from a list of trees
 * \param max_order order of the largest rooted tree in list_of_tree
 * \param[in,out] dags list of entry points in the dag according to the order we want.
 */
void wDAG_aff::generate_dags(int max_order, list<rootedTree> * list_of_tree,
                   list< list<vertex_wDAG_aff *> > * dags, IntervalMatrix * A, IntervalVector * b)
{
#ifdef __PRINT_DEBUG__
    cout << "******************************************************************************************" << endl;
    cout << __func__ << endl;
    cout << "max_order : " << max_order << endl;
#endif

    list<rootedTree>::iterator it_tree;
    for (it_tree = list_of_tree->begin(); it_tree != list_of_tree->end(); it_tree++) {
        rootedTree * current_tree = new rootedTree(*it_tree);
        list <list<vertex_wDAG_aff *> >::iterator it_dag;

        //We get the list of entry point where current_tree has to be inserted (by construction, it is not already inside this list)
        it_dag = dags->begin();
        std::advance(it_dag, current_tree->size()-2);

        //creation of the vertex corresponding to current_tree
        vertex_wDAG_aff * current_vertex = new vertex_wDAG_aff(current_tree);

        //we add a new vertex in the list of entry points of order equal to the one of current_tree
        (*it_dag).push_back(current_vertex);
        int i = 0;

        //Recovery of the children of current_tree
        list<rootedTree> children = current_tree->get_children();
        list<rootedTree>::iterator it_children;

        //We add the edges corresponding to the children of current_tree
        //To do it, we need to find each child in its corresponding list of entry points
        for (it_children = children.begin(); it_children != children.end(); it_children++) {
            list <list<vertex_wDAG_aff *> >::iterator it_dag2;
            //We get the list of entry points with order corresponding to the current child
            it_dag2 = dags->begin();

            std::advance(it_dag2, (*it_children).size()-2);
            vertex_wDAG_aff * child = find_vertex(&(*it_dag2), &(*it_children));

            //edge is inserted in the list.
            current_vertex->add_edge(child);
        }
    }
}


/**
 * \brief intern recursive function for the generation of the rooted trees of order max_order.
 * User should call generate(const int max_order) instead which calls this function with the right initialized parameters.
 * \param p current order of the generated tree temp
 * \param s represent the current root in temp where the next vertex has to be added
 * \param L order of the subtree of temp which has s as root
 * \param temp current tree that is constructed
 * \param[in,out] res list the trees generated during the recursive call of this function
 */
void generate_rec(int p, int s, int L, int max_order, rootedTree temp, list<rootedTree> *res)
{
    if (p <= max_order) {
        if (L == 0) {
            temp.add(p-1);
        } else if (temp.get(p-L) < s) {
            temp.add(temp.get(s));
        } else {
            temp.add(L + temp.get(p - L));
        }
        res->push_back(temp);
	generate_rec(p+1, s, L, max_order, temp, res);
        while (temp.get(p)>1) {
            s = temp.get(p);
            temp.set(p, temp.get(s));
	    res->push_back(temp);
            generate_rec(p+1, s, p-s, max_order, temp, res);
        }
    }
}

/**
 * \brief Generation of the set of rooted trees of order max_order.
 * \param max_order maximum number of vertices for the trees.
 */
list<rootedTree> * wDAG_aff::generate(const int max_order)
{
    list<rootedTree> * result = new list<rootedTree>();
    rootedTree ptemp;
    ptemp.add(0);
    generate_rec(1, 0, 0, max_order, ptemp, result);
    result->sort();
    return result;
}

//save the generated dags to a file
bool wDAG_aff::dags_to_file(const std::string &path){

    std::map<vertex_wDAG_aff*, int > map_of_vertices;
    int index = 0;

    std::ofstream FILE;
    FILE.open(path.c_str(), std::ios::out | std::ofstream::binary);

    if (!FILE.is_open()){
        return false;
    }


    //sauvegarde de la taille de la liste de dags
    int nb_list_of_dags = _list_of_dags->size();
//     cout << "nb list of dags : " << nb_list_of_dags << endl;
    FILE.write(reinterpret_cast<const char *>(&nb_list_of_dags), sizeof(nb_list_of_dags));

    //calcul de la taille d'un noeud


    list<list<vertex_wDAG_aff *> >::iterator it_list_of_dags = _list_of_dags->begin();
    while (it_list_of_dags != _list_of_dags->end()){
        int nb_current = it_list_of_dags->size();
        FILE.write(reinterpret_cast<const char *>(&nb_current), sizeof(nb_current));

//         cout << "\tnb current : " << nb_current << endl;

        list<vertex_wDAG_aff *>::iterator it_vertices = it_list_of_dags->begin();
        while (it_vertices != it_list_of_dags->end()){
            map_of_vertices[*it_vertices] = index;
            index++;
//             cout << "tree : " << (*it_vertices)->get_tree()->to_string() << endl;
            //sauvegarde des infos du noeud
//             cout << (*it_vertices)->get_tree()->to_string() << endl;
            int current_order = (*it_vertices)->get_order();
            int current_symetry = (*it_vertices)->get_symetry();
            int current_density = (*it_vertices)->get_density();
            Interval current_alpha = *(*it_vertices)->get_alpha();
//             cout << "\tordre : " << current_order << endl;
//             cout << "\tsymetrie : " << current_symetry << endl;
//             cout << "\tdensité : " << current_density << endl;
//             cout << "\talpha : " << current_alpha << endl;
            FILE.write(reinterpret_cast<const char *>(&( current_order )), sizeof(int));
            FILE.write(reinterpret_cast<const char *>(&( current_symetry )), sizeof(int));
            FILE.write(reinterpret_cast<const char *>(&( current_density )), sizeof(int));
            FILE.write(reinterpret_cast<const char *>(&( current_alpha )), sizeof(Interval));

            //sauvegarde de la taille de la liste des arcs
            int nb_sons = (*it_vertices)->get_list_of_edge()->size();
            FILE.write(reinterpret_cast<const char *>(&nb_sons), sizeof(nb_sons));
            list<pair<vertex_wDAG_aff*, int > >::iterator it_edge = (*it_vertices)->get_list_of_edge()->begin();
            while (it_edge != (*it_vertices)->get_list_of_edge()->end()){
//                 cout << "\t\tfils : " << map_of_vertices[it_edge->first] << endl;
                FILE.write(reinterpret_cast<const char *>(&map_of_vertices[it_edge->first]), sizeof(int));
                FILE.write(reinterpret_cast<const char *>(&it_edge->second), sizeof(int));

                it_edge++;
            }
            it_vertices++;
        }
        it_list_of_dags++;
    }
    return true;
}

bool wDAG_aff::file_to_dags(int order, list<list<ibex::vertex_wDAG_aff *> > * result, const std::string &path){
    ifstream FILE;
    FILE.open(path.c_str(), std::ios::in | std::ifstream::binary);

    if (!FILE.is_open()){
        return false;
    }
    std::map<int, vertex_wDAG_aff* > map_of_vertices;
    int index = 0;

    int nb_list_of_dags = 0;
    FILE.read(reinterpret_cast<char *>(&nb_list_of_dags), sizeof(nb_list_of_dags));
//     cout << "nombre de listes de dags : " << nb_list_of_dags << endl;
    if (nb_list_of_dags < order)
        return false;
    for (int i = 0; i < order; i++){
        list<ibex::vertex_wDAG_aff *> temp;
        //récupération de chaque liste
        int nb_current = 0;
        FILE.read(reinterpret_cast<char *>(&nb_current), sizeof(nb_current));
//         cout << "\tnombre de sommets dans la liste courante : " << nb_current << endl;
        for (int j = 0; j < nb_current; j++){
            int current_order = 0;
            int current_symetry = 0;
            int current_density = 0;
            Interval current_alpha = Interval(0);
            FILE.read(reinterpret_cast<char *>(&current_order), sizeof(current_order));
            FILE.read(reinterpret_cast<char *>(&current_symetry), sizeof(current_symetry));
            FILE.read(reinterpret_cast<char *>(&current_density), sizeof(current_density));
            FILE.read(reinterpret_cast<char *>(&current_alpha), sizeof(current_alpha));
            vertex_wDAG_aff * current = new vertex_wDAG_aff(current_order, current_symetry, current_density, &current_alpha);
//             cout << "\tadresse du sommet courant : " << current << endl;
            int nb_sons = 0;
            FILE.read(reinterpret_cast<char *>(&nb_sons), sizeof(nb_sons));
//             cout << "\t\tnombre de fils : " << nb_sons << endl;
            list<pair<vertex_wDAG_aff *, int > > current_list_of_edge;
            for (int k = 0; k < nb_sons; k++){
                int current_index = 0;
                FILE.read(reinterpret_cast<char *>(&current_index), sizeof(current_index));
                int weight = 0;
                FILE.read(reinterpret_cast<char *>(&weight), sizeof(weight));
                pair<vertex_wDAG_aff *, int > current_son(map_of_vertices[current_index], weight);
//                 cout << "\t\t\tadresse du fils : " << map_of_vertices[current_index] << endl;
                current_list_of_edge.push_back(current_son);
            }
            current->set_list_of_edge(current_list_of_edge);
            temp.push_back(current);
            map_of_vertices[index] = current;
            index++;
        }
        result->push_back(temp);
    }

    return true;
}



} //namespace
