#include "oslom_net_global.h"

int oslom_net_global::check_intersection(module_collection & Mcoll)
{
  paras.print_flag_subgraph = false;

  std::deque<int> to_check;
  for (auto itM = Mcoll.module_bs.begin(); itM != Mcoll.module_bs.end(); ++itM)
    to_check.push_back(itM->first);

  return check_intersection(to_check, Mcoll);
}

int oslom_net_global::check_intersection(std::deque<int> & to_check, module_collection & Mcoll)
{
  std::set<std::pair<int, int>> pairs_to_check;

  for (auto itM = to_check.begin(); itM != to_check.end(); ++itM)
  {
    if (Mcoll.module_bs.find(*itM) != Mcoll.module_bs.end())
    {
      std::deque<int> & c = Mcoll.modules[*itM];
      std::map<int, int> com_ol;						// it maps the index of the modules into the overlap (overlap=number of overlapping nodes)

      for (int i : c)
      {
        for (auto itj = Mcoll.memberships[i].begin(); itj != Mcoll.memberships[i].end(); ++itj)
        {
          int_histogram(*itj, com_ol);
        }
      }

      for (auto cit = com_ol.begin(); cit != com_ol.end(); ++cit)
      {
        if (cit->first != *itM)
        {
          if (double(cit->second) / std::min(Mcoll.modules[cit->first].size(), c.size()) > paras.check_inter_p)		// they have a few nodes in common
          {
            pairs_to_check.insert(std::make_pair(std::min(*itM, cit->first), std::max(*itM, cit->first)));
          }
        }
      }
    }
  }

  return fusion_intersection(pairs_to_check, Mcoll);
}

int oslom_net_global::fusion_intersection(std::set<std::pair<int, int>> & pairs_to_check, module_collection & Mcoll) {
  std::cout << "pairs to check: " << pairs_to_check.size() << std::endl;

  std::deque<int> new_insertions;

  for (std::set<std::pair<int, int>>::iterator ith = pairs_to_check.begin(); ith != pairs_to_check.end(); ith++) if (ith->first < ith->second)
  {
    if (Mcoll.module_bs.find(ith->first) != Mcoll.module_bs.end())	if (Mcoll.module_bs.find(ith->second) != Mcoll.module_bs.end())
    {
      //		first, you need to check if both the modules in the pair are still in mcoll
      std::deque<int> & a1 = Mcoll.modules[ith->first];
      std::deque<int> & a2 = Mcoll.modules[ith->second];
      int min_s = std::min(a1.size(), a2.size());

      std::deque<int> group_intsec;
      set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_intsec));

      //		if they are, you need to check if they are not almost equal.

      if (double(group_intsec.size()) / min_s >= paras.coverage_inclusion_module_collection)
      {
        int em = ith->first;
        if (a1.size() < a2.size())
          em = ith->second;
        else if (a1.size() == a2.size() && Mcoll.module_bs[ith->first] > Mcoll.module_bs[ith->second])
          em = ith->second;

        Mcoll.erase(em);
      }
      else
        decision_fusion_intersection(ith->first, ith->second, new_insertions, Mcoll, double(group_intsec.size()) / min_s);
    }
  }
  if (!new_insertions.empty())
    return check_intersection(new_insertions, Mcoll);

  return 0;
}

bool oslom_net_global::decision_fusion_intersection(int ai1, int ai2, std::deque<int> & new_insertions, module_collection & Mcoll, double prev_over_percentage)
{
  std::deque<int> & a1 = Mcoll.modules[ai1];
  std::deque<int> & a2 = Mcoll.modules[ai2];

  std::deque<int> group;
  set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group));

  if (int(group.size()) != dim)			//******************  sub_graph_module stuff   ******************
  {
    std::deque<std::deque<int>> link_per_node;
    std::deque<std::deque<std::pair<int, double>>> weights_per_node;
    set_subgraph(group, link_per_node, weights_per_node);
    oslom_net_global sub_graph_module(link_per_node, weights_per_node, group);

    std::deque<std::deque<int>> A;
    A.push_back(a1);
    A.push_back(a2);

    sub_graph_module.translate(A);

    std::deque<int> grc1;
    double bs = sub_graph_module.CUP_check(A[0], grc1);
    std::deque<int> grc2;
    bs = sub_graph_module.CUP_check(A[1], grc2);

    std::deque<int> unions_grcs;
    set_union(grc1.begin(), grc1.end(), grc2.begin(), grc2.end(), back_inserter(unions_grcs));

    if (unions_grcs.size() <= paras.coverage_percentage_fusion_or_submodules*group.size()) {		// in such a case you can take the fusion (if it's good)
      /* actually the right check should be  unions_grcs.size() > paras.coverage_percentage_fusion_or_submodules*group_2.size()
         but this would require more time - it should not make a big difference anyway */

      std::deque<int> group_2;
      bs = CUP_check(group, group_2);

      if (group_2.size() > paras.coverage_percentage_fusion_left * group.size())
      {
        Mcoll.erase(ai1);
        Mcoll.erase(ai2);

        int_matrix _A_;
        std::deque<double> _bss_;
        _A_.push_back(group_2);
        _bss_.push_back(bs);
        check_minimality_all(_A_, _bss_, Mcoll);

        return true;
      }
      else
        return false;
    }

    sub_graph_module.deque_id(grc1);
    sub_graph_module.deque_id(grc2);

    std::deque<int> cg1;
    double bs__1 = CUP_check(grc1, cg1);

    std::deque<int> cg2;
    double bs__2 = CUP_check(grc2, cg2);

    std::deque<int> inters;
    set_intersection(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(), back_inserter(inters));

    std::deque<int> unions;
    set_union(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(), back_inserter(unions));

    if (double(inters.size()) / std::min(cg1.size(), cg2.size()) < prev_over_percentage - 1e-4)
    {
      if (!cg1.empty() && !cg2.empty() && (unions.size() > paras.coverage_percentage_fusion_left * group.size()))
      {
        Mcoll.erase(ai1);
        Mcoll.erase(ai2);

        int newi;

        Mcoll.insert(cg1, bs__1, newi);
        new_insertions.push_back(newi);

        Mcoll.insert(cg2, bs__2, newi);
        new_insertions.push_back(newi);

        //cout<<"pruned module"<< std::endl;

        return true;
      }
    }
  }
  return false;
}
