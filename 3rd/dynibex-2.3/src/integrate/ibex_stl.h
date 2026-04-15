/* ============================================================================
 * D Y N I B E X - Definition of the Signal Temporal Logic (STL)
 * ============================================================================
 * Copyright   : ENSTA
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Joris Tillet
 * Created     : Feb 06, 2025
 * Sponsored   : This research benefited from the support of the CIEDS within the
 * STARTS Project.
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_STL_H
#define IBEX_STL_H

namespace ibex{

/** subset: [simu(t)] \subset [box]? true - false if \cap == \emptyset - maybe otherwise (\cap != \emptyset) */
inline BoolInterval subset(simulation& simu, float t, IntervalVector box) {
  if (box.size() != simu.nb_var || t < 0 || t > simu.time_T)
    {
      std::ostringstream msg;
      msg << "error in subset function: box.size() (" << box.size() << ") != simu.nb_var (" << simu.nb_var << ") or t (" << t <<") < 0 or t > simu.time_T (" << simu.time_T << ")";
      throw std::invalid_argument(msg.str());
      }
  IntervalVector x = simu.get(t);
  if (x.is_subset(box))
    return YES;
  if ((x & box).is_empty())
    return NO;
  
  x = simu.get_tight(t);
  if (x.is_subset(box))
    return YES;
  if ((x & box).is_empty())
    return NO;
  return MAYBE;
}

/** globally_subset: same as subset but on the time interval [t] */
inline BoolInterval globally_subset(simulation& simu, Interval t, IntervalVector box) {
  if (box.size() != simu.nb_var || t.lb() < 0 || t.ub() > simu.time_T || simu.list_solution_g.empty())
    {
      std::ostringstream msg;
      msg << "error in globally_subset function: box.size() (" << box.size() << ") != simu.nb_var (" << simu.nb_var << ") or t.lb() (" << t.lb() <<") < 0 or t.ub() (" << t.ub() << ") > simu.time_T (" << simu.time_T << ") or simu.list_solution_g.empty()";
      throw std::invalid_argument(msg.str());
      }

    std::list<solution_g>::iterator iterator_list;
    BoolInterval result = YES;
    for (iterator_list = simu.list_solution_g.begin(); iterator_list != simu.list_solution_g.end(); iterator_list++)
    {
      if (!(iterator_list->time_j & t).is_empty()) {
        IntervalVector y_temp = *(iterator_list->box_j1);
        if ((y_temp & box).is_empty())
          return NO;
        if (result == YES && !y_temp.is_subset(box))
          result = MAYBE;
      }
    }
    return result;
}

/** globally_not_subset */
inline BoolInterval globally_not_subset(simulation& simu, Interval t, IntervalVector box){
  if (box.size() != simu.nb_var || t.lb() < 0 || t.ub() > simu.time_T || simu.list_solution_g.empty())
    {
      std::ostringstream msg;
      msg << "error in globally_not_subset function: box.size() (" << box.size() << ") != simu.nb_var (" << simu.nb_var << ") or t.lb() (" << t.lb() <<") < 0 or t.ub() (" << t.ub() << ") > simu.time_T (" << simu.time_T << ") or simu.list_solution_g.empty()";
      throw std::invalid_argument(msg.str());
      }

    std::list<solution_g>::iterator iterator_list;
    BoolInterval result = YES;
    for (iterator_list = simu.list_solution_g.begin(); iterator_list != simu.list_solution_g.end(); iterator_list++)
    {
      if (!(iterator_list->time_j & t).is_empty()) {
        IntervalVector y_temp = *(iterator_list->box_j1);
        if (y_temp.is_subset(box))
          return NO;
        if (result == YES && !(y_temp & box).is_empty())
          result = MAYBE;
      }
    }
    return result;
}

/** finally_subset: there exists a t' in the time interval [t] such that subset is true */
inline BoolInterval finally_subset(simulation& simu, Interval t, IntervalVector box){
  if (box.size() != simu.nb_var || t.lb() < 0 || t.ub() > simu.time_T || simu.list_solution_g.empty())
    {
      std::ostringstream msg;
      msg << "error in finally_subset function: box.size() (" << box.size() << ") != simu.nb_var (" << simu.nb_var << ") or t.lb() (" << t.lb() <<") < 0 or t.ub() (" << t.ub() << ") > simu.time_T (" << simu.time_T << ") or simu.list_solution_g.empty()";
      throw std::invalid_argument(msg.str());
      }

    std::list<solution_g>::iterator iterator_list;
    BoolInterval result = NO;
    for (iterator_list = simu.list_solution_g.begin(); iterator_list != simu.list_solution_g.end(); iterator_list++)
    {
      if (!(iterator_list->time_j & t).is_empty()) {
        IntervalVector y_temp = *(iterator_list->box_j1);
        if (y_temp.is_subset(box))
          return YES;
        if (result == NO && !(y_temp & box).is_empty())
          result = MAYBE;
      }
    }
    return result;
}

/** finally_not_subset: same as finally_subset but with a t' such that subset is false */
inline BoolInterval finally_not_subset(simulation& simu, Interval t, IntervalVector box) {
  if (box.size() != simu.nb_var || t.lb() < 0 || t.ub() > simu.time_T || simu.list_solution_g.empty())
    {
      std::ostringstream msg;
      msg << "error in finally_not_subset function: box.size() (" << box.size() << ") != simu.nb_var (" << simu.nb_var << ") or t.lb() (" << t.lb() <<") < 0 or t.ub() (" << t.ub() << ") > simu.time_T (" << simu.time_T << ") or simu.list_solution_g.empty()";
      throw std::invalid_argument(msg.str());
      }

    std::list<solution_g>::iterator iterator_list;
    BoolInterval result = NO;
    for (iterator_list = simu.list_solution_g.begin(); iterator_list != simu.list_solution_g.end(); iterator_list++)
    {
      if (!(iterator_list->time_j & t).is_empty()) {
        IntervalVector y_temp = *(iterator_list->box_j1);
        if ((y_temp & box).is_empty())
          return YES;
        if (result == NO && !y_temp.is_subset(box))
          result = MAYBE;
      }
    }
    return result;
}

/** until: same as globally_subset with box1 until a time t' within time interval [t] such that subset is true with box2 on time t' */
inline BoolInterval until(simulation& simu, Interval t, IntervalVector box1, IntervalVector box2){
  if (box1.size() != simu.nb_var || box2.size() != simu.nb_var || t.lb() < 0 || t.ub() > simu.time_T || simu.list_solution_g.empty())
    {
      std::ostringstream msg;
      msg << "error in finally_not_subset function: box.size() (" << box1.size() << " or " << box2.size() << ") != simu.nb_var (" << simu.nb_var << ") or t.lb() (" << t.lb() <<") < 0 or t.ub() (" << t.ub() << ") > simu.time_T (" << simu.time_T << ") or simu.list_solution_g.empty()";
      throw std::invalid_argument(msg.str());
      }

    Interval tprime = simu.has_crossed_when(box2);
    // std::cout << tprime << std::endl; DEBUG
    if (!tprime.is_subset(t))
      return MAYBE; // TODO check NO
    return globally_subset(simu, Interval(0,tprime.ub()), box1);
}

}

#endif
