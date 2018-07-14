// File: community.h
// -- community detection header file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This file is part of Louvain algorithm.
//
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details
#pragma once

#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "graph_binary.h"

class Community
{
public:
  std::vector<double> m_neighborWeights;
  std::vector<unsigned int> m_neighPositions;
  unsigned int neigh_last;

  Graph m_graph; // network to compute communities for
  int m_Size; // nummber of nodes in the network and size of all vectors
  std::vector<int> m_node2communityMap; // community to which each node belongs
  std::vector<double> in, tot; // used to compute the modularity participation of each community

  // number of pass for one level computation
  // if -1, compute as many passes as needed to increase modularity
  int m_numberPasses;

  // a new pass is computed if the last one has generated an increase
  // greater than min_modularity
  // if 0. even a minor increase is enough to go for one more pass
  double m_minModularity;

  // constructors:
  // reads graph from file using graph constructor
  // type defined the weighted/unweighted status of the graph file
  Community(char *filename, char *filename_w, int type, int nb_pass, double min_modularity);
  // copy graph
  Community(Graph g, int nb_pass, double min_modularity);

  // initiliazes the partition with something else than all nodes alone
  void init_partition(char *filename_part);

  // display the community of each node
  void display();

  // remove the node from its current community with which it has dnodecomm links
  inline void remove(int node, int comm, double dnodecomm);

  // insert the node in comm with which it shares dnodecomm links
  inline void insert(int node, int comm, double dnodecomm);

  // compute the gain of modularity if node where inserted in comm
  // given that node has dnodecomm links to comm.  The formula is:
  // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
  // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
  // where In(comm)    = number of half-links strictly inside comm
  //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
  //       d(node,com) = number of links from node to comm
  //       deg(node)   = node degree
  //       m           = number of links
  inline double modularity_gain(int node, int comm, double dnodecomm, double w_degree);

  // compute the set of neighboring communities of node
  // for each community, gives the number of links from node to comm
  void neigh_comm(unsigned int node);

  // compute the modularity of the current partition
  double modularity();

  // displays the graph of communities as computed by one_level
  void partition2graph();
  // displays the current partition (with communities renumbered from 0 to k-1)
  void display_partition();

  // generates the binary graph of communities as computed by one_level
  Graph partition2graph_binary();

  // compute communities of the graph for one level
  // return true if some nodes have been moved
  bool one_level();
};

inline void Community::remove(int node, int comm, double dnodecomm)
{
  assert(node >= 0 && node < m_Size);

  tot[comm] -= m_graph.getWeightedDegree(node);
  in[comm] -= 2 * dnodecomm + m_graph.getNumberOfSelfLoops(node);
  m_node2communityMap[node] = -1;
}

inline void Community::insert(int node, int comm, double dnodecomm)
{
  assert(node >= 0 && node < m_Size);

  tot[comm] += m_graph.getWeightedDegree(node);
  in[comm] += 2 * dnodecomm + m_graph.getNumberOfSelfLoops(node);
  m_node2communityMap[node] = comm;
}

inline double Community::modularity_gain(int node, int comm, double dnodecomm, double w_degree)
{
  assert(node >= 0 && node < m_Size);

  double totc = tot[comm];
  double degc = w_degree;
  double m2 = m_graph.total_weight;
  double dnc = dnodecomm;

  return (dnc - totc * degc / m2);
}

#endif // COMMUNITY_H
