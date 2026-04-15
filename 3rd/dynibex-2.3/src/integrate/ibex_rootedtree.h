/* ============================================================================
 * D Y N I B E X - Definition of the RootedTree computation
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created     : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */


#ifndef ROOTEDTREE_H
#define ROOTEDTREE_H

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

typedef std::vector<int> my_list_t;

namespace ibex {

/**
 * \brief this class represents a rooted tree with an adjacency list.
 * _adjacency_list[i] = j means that the parent of i is j.
 * for the root, its parent is set to 0.
 */
class rootedTree
{
private:
    //Rooted tree represented by an adjacency list.
    my_list_t *_adjacency_list;

public:
    //constructors
    rootedTree();
    rootedTree(const rootedTree& other);

    //destructors
    ~rootedTree();


    //getters
    inline my_list_t get_adjacency_list() {
        return *_adjacency_list;
    };

    /**
     * \brief add an element at the end of the list
     */
    void add(int);

    /**
     * \brief obtain the parent of the i^th element
     */
    inline int get(int i) const {
        return _adjacency_list->at(i);
    };

    inline void clear(){
      _adjacency_list->clear();
    }

    /**
     * \brief For 2 rootedTrees a and b, a < b iff (a.size() < b.size() or with lexicographic order, a < b)
     */
    bool operator<(const rootedTree &) const;

    /**
     * \brief test if 2 trees are equals
     */
    bool operator==(rootedTree &);

    /**
     * \brief get the number of children of the argument
     */
    int number_of_children(int);

    /**
     * \brief get the arg1^th child of vertex arg2
     */
    int get_child(int, int);

    /**
     * \brief return the list of subtrees of this.
     */
    list<rootedTree> get_children();

    /**
     * \brief return the size of the list (e.g. the number of vertices)
     */
    inline int size()const {
        return this->_adjacency_list->size();
    };


    /**
     * \brief put arg2 to the arg1 position in the adjacency list.
     */
    void set(int, int);

    string to_string();

};


} //namespace
#endif // ROOTEDTREE_H
