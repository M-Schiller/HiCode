#pragma once

#ifndef EGOCENTRIC_UNDIR_H
#define EGOCENTRIC_UNDIR_H

#include <utility>
#include <deque>
#include "module_collection.h"
#include "louvain_oslomnet.h"

class egocentric_net : public oslomnet_louvain
{
public:

  egocentric_net() : oslomnet_louvain() {};
  ~egocentric_net() = default;

  int collect_ego_groups_once(int_matrix &);

private:

  int add_this_egomodules(int node, module_collection & Mego);
};

inline int egocentric_net::collect_ego_groups_once(int_matrix & E)
{
  module_collection Mego(dim);
  paras.print_flag_subgraph = false;

  for (unsigned i = 0; i < dim; i++)
  {
    if (vertices[i]->stub_number > 10)
    {
      add_this_egomodules(i, Mego);
    }
  }

  Mego.erase_included();
  Mego.set_partition(E);
  std::ofstream pout("tpp");
  print_id(E, pout);

  return 0;
}

inline int egocentric_net::add_this_egomodules(int node, module_collection & Mego)
{
  std::deque<std::deque<int>> link_per_node;
  std::deque<std::deque<std::pair<int, double>> > weights_per_node;

  std::deque<int> group;
  for (int j = 0; j < vertices[node]->links->size(); j++)
  {
    group.push_back(vertices[node]->links->l[j]);
  }

  std::cout << "........ " << id_of(node) << " " << node << std::endl;

  set_subgraph(group, link_per_node, weights_per_node);
  oslomnet_louvain ego_subgraph;
  ego_subgraph.set_graph(link_per_node, weights_per_node, group);
  int_matrix A;
  ego_subgraph.collect_raw_groups_once(A);

  if (A.size() > 1)
  {
    for (unsigned i = 0; i < A.size(); i++) if (A[i].size() > 1)
    {
      ego_subgraph.deque_id(A[i]);
      A[i].push_back(node);
      std::set<int> sA;
      deque_to_set(A[i], sA);
      set_to_deque(sA, A[i]);

      std::cout << "......... A[i] " << std::endl;
      print_id(A[i], std::cout);

      Mego.insert(A[i], 1.);
    }

    std::cout << "*************************************************" << std::endl;
  }
  return 0;
}
#endif // EGOCENTRIC_UNDIR_H
