/* ============================================================================
 * D Y N I B E X - Definition of the RootedTree computation
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Olivier Mullier
 * Created    : Jul 25, 2018
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#include "ibex_rootedtree.h"
#include "ibex_Array.h"

using namespace std;

namespace ibex{

rootedTree::rootedTree()
{
  _adjacency_list = new my_list_t();
}

rootedTree::rootedTree(const rootedTree& other)
{
  _adjacency_list = new my_list_t(*other._adjacency_list);
}

rootedTree::~rootedTree()
{
//     cout << __func__ << endl;
  if (_adjacency_list)
    delete _adjacency_list;
}

void rootedTree::add(int n)
{
  _adjacency_list->push_back(n);
}

int rootedTree::get_child(int i, int n)
{
  int cpt = 0;
  for (int j = 1; j < size(); j++){
    if (get(j)==n){
      cpt++;
      if (cpt == i){
	return j;
      }
    }
  }
  return 0;
}

list< rootedTree > rootedTree::get_children()
{
  list< rootedTree > result;
  if (size() == 2){
    return result;
  }
  int last_child = 0, current_child_index = 0;

  rootedTree current_tree;
  current_tree.add(0);
  current_tree.add(0);
  if (size()==3){
    result.push_back(current_tree);
    return result;
  }
  for (int i = 3; i < size(); i++){
    current_child_index++;
    if (this->get(i)>1){
      current_tree.add(this->get(i)-1-last_child);
      if (i+1 == size()){
	result.push_back(current_tree);
      }
    }
    else {//new subtree
      last_child = current_child_index;
      result.push_back(current_tree);
      current_tree.clear();
//       rootedTree * temp = new rootedTree();
//       current_tree = *temp;
      current_tree.add(0);
      current_tree.add(0);
//       delete temp;
    }
    if (i+1 == size() && current_tree.size()==2)
      result.push_back(current_tree);
  }
  return result;
}


int rootedTree::number_of_children(int n)
{
  if (n == 0) return 1;
  int res = 0;
  for(int i = 1;i < n; i++){
    if (i == n)
      res++;
  }
  return res;
}

bool rootedTree::operator<(const rootedTree& other) const
{
  if (this->size() == other.size()){
    for (int i = 1; i < size(); i++){
      if (this->get(i) > other.get(i)){
	return false;
      } else if (this->get(i) < other.get(i)){
	return true;
      }
    }
    return false;
  }
  return (this->size() < other.size());
}

bool rootedTree::operator==(rootedTree& other)
{
  if (size() == other.size()){
    for (int i = 1; i < size(); i++){
      if (get(i) != other.get(i))
	return false;
    }
    return true;
  }
  return false;
}

void rootedTree::set(int i, int n)
{
  assert(i < size());
  _adjacency_list->operator[](i) = n;
}

string rootedTree::to_string()
{
  string res("");
  if (this->size() > 1){
    res += "(";
    for(int i = 1; i < this->size(); i++){
      string s = std::to_string(this->get(i));
      res += s;
    }
    res += ")";
  } else {
    res += "Empty list";
  }
  return res;
//   std::stringstream sstm;
//   if (this->size() > 1){
//     sstm << "(";
//     for(int i = 1; i < this->size(); i++){
//       sstm << this->get(i);
//     }
//     sstm << ")";
//   } else {
//    sstm << "Empty list";
//   }
//   return sstm.str();
}

} //namespace
