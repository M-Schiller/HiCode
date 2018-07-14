// File: community.h
// -- community detection source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "community.h"

Community::Community(char * filename, char * filename_w, int type, int nb_pass, double min_modularity)
{
  m_graph = Graph(filename, filename_w, type);
  m_Size = m_graph.nb_nodes;

  m_neighborWeights.resize(m_Size, -1);
  m_neighPositions.resize(m_Size);
  neigh_last = 0;

  m_node2communityMap.resize(m_Size);
  in.resize(m_Size);
  tot.resize(m_Size);

  for (int i = 0; i < m_Size; i++)
  {
    m_node2communityMap[i] = i;
    tot[i] = m_graph.getWeightedDegree(i);
    in[i] = m_graph.getNumberOfSelfLoops(i);
  }

  m_numberPasses = nb_pass;
  m_minModularity = min_modularity;
}

Community::Community(Graph g, int nb_pass, double min_modularity)
{
  m_graph = g;
  m_Size = g.nb_nodes;

  m_neighborWeights.resize(m_Size, -1);
  m_neighPositions.resize(m_Size);
  neigh_last = 0;

  m_node2communityMap.resize(m_Size);
  in.resize(m_Size);
  tot.resize(m_Size);

  for (int i = 0; i < m_Size; i++)
  {
    m_node2communityMap[i] = i;
    in[i] = g.getNumberOfSelfLoops(i);
    tot[i] = g.getWeightedDegree(i);
  }

  m_numberPasses = nb_pass;
  m_minModularity = min_modularity;
}

void Community::init_partition(char * filename) {
  std::ifstream finput;
  finput.open(filename, std::fstream::in);

  // read partition
  while (!finput.eof())
  {
    unsigned int node, comm;
    finput >> node >> comm;

    if (finput)
    {
      int old_comm = m_node2communityMap[node];
      neigh_comm(node);

      remove(node, old_comm, m_neighborWeights[old_comm]);

      unsigned int i;
      for (i = 0; i < neigh_last; ++i)
      {
        unsigned int best_comm = m_neighPositions[i];
        float best_nblinks = m_neighborWeights[m_neighPositions[i]];
        if (best_comm == comm)
        {
          insert(node, best_comm, best_nblinks);
          break;
        }
      }
      if (i == neigh_last)
      {
        insert(node, comm, 0);
      }
    }
  }
  finput.close();
}

void Community::display()
{
  for (int i = 0; i < m_Size; i++)
  {
    std::cerr << " " << i << "/" << m_node2communityMap[i] << "/" << in[i] << "/" << tot[i];
  }
  std::cerr << std::endl;
}

double Community::modularity()
{
  double q = 0.;
  const double m2 = m_graph.total_weight;

  for (int i = 0; i < m_Size; i++)
  {
    if (tot[i] > 0)
    {
      q += in[i] / m2 - (tot[i] / m2)*(tot[i] / m2);
    }
  }
  return q;
}

void Community::neigh_comm(unsigned int node)
{
  for (unsigned int i = 0; i < neigh_last; ++i)
  {
    m_neighborWeights[m_neighPositions[i]] = -1;
  }
  neigh_last = 0;

  std::pair<std::vector<unsigned int>::iterator, std::vector<float>::iterator> p = m_graph.neighbors(node);

  const unsigned int deg = m_graph.getDegree(node);

  m_neighPositions[0] = m_node2communityMap[node];
  m_neighborWeights[m_neighPositions[0]] = 0;
  neigh_last = 1;

  for (unsigned int i = 0; i < deg; i++) {
    unsigned int neigh = *(p.first + i);
    unsigned int neigh_comm = m_node2communityMap[neigh];
    double neigh_w = (m_graph.weights.empty()) ? 1. : *(p.second + i);

    if (neigh != node) {
      if (m_neighborWeights[neigh_comm] == -1) {
        m_neighborWeights[neigh_comm] = 0.;
        m_neighPositions[neigh_last++] = neigh_comm;
      }
      m_neighborWeights[neigh_comm] += neigh_w;
    }
  }
}

void Community::partition2graph() {
  std::vector<int> renumber(m_Size);
  for (int node = 0; node < m_Size; ++node)
  {
    renumber[m_node2communityMap[node]]++;
  }

  int final = 0;
  for (int i = 0; i < m_Size; ++i)
  {
    if (renumber[i] != -1)
    {
      renumber[i] = final++;
    }
  }

  for (int i = 0; i < m_Size; ++i) {
    std::pair<std::vector<unsigned int>::iterator, std::vector<float>::iterator> p = m_graph.neighbors(i);

    int deg = m_graph.getDegree(i);
    for (int j = 0; j < deg; ++j)
    {
      int neigh = *(p.first + j);
      std::cout << renumber[m_node2communityMap[i]] << " " << renumber[m_node2communityMap[neigh]] << std::endl;
    }
  }
}

void Community::display_partition()
{
  std::vector<int> renumber(m_Size, -1);
  for (int node = 0; node < m_Size; node++) {
    renumber[m_node2communityMap[node]]++;
  }

  int final = 0;
  for (int i = 0; i < m_Size; i++)
    if (renumber[i] != -1)
      renumber[i] = final++;

  for (int i = 0; i < m_Size; i++)
    std::cout << i << " " << renumber[m_node2communityMap[i]] << std::endl;
}

Graph Community::partition2graph_binary()
{
  // Renumber communities
  std::vector<int> renumber(m_Size, -1);
  for (int node = 0; node < m_Size; node++)
  {
    renumber[m_node2communityMap[node]]++;
  }

  int final = 0;
  for (int i = 0; i < m_Size; i++)
  {
    if (renumber[i] != -1)
    {
      renumber[i] = final++;
    }
  }

  // Compute communities
  std::vector<std::vector<int>> comm_nodes(final);
  for (int node = 0; node < m_Size; node++)
  {
    comm_nodes[renumber[m_node2communityMap[node]]].push_back(node);
  }

  // Compute weighted graph
  Graph g2;
  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(comm_nodes.size());

  int comm_deg = comm_nodes.size();
  for (int comm = 0; comm < comm_deg; comm++) {
    std::map<int, float> m;
    std::map<int, float>::iterator it;

    int comm_size = comm_nodes[comm].size();
    for (int node = 0; node < comm_size; node++) {
      std::pair<std::vector<unsigned int>::iterator, std::vector<float>::iterator> p = m_graph.neighbors(comm_nodes[comm][node]);
      int deg = m_graph.getDegree(comm_nodes[comm][node]);
      for (int i = 0; i < deg; i++)
      {
        int neigh = *(p.first + i);
        int neigh_comm = renumber[m_node2communityMap[neigh]];
        double neigh_weight = m_graph.weights.empty() ? 1. : *(p.second + i);

        it = m.find(neigh_comm);
        if (it == m.end())
        {
          m.emplace(neigh_comm, neigh_weight);
        }
        else
          it->second += neigh_weight;
      }
    }
    g2.degrees[comm] = (comm == 0) ? m.size() : g2.degrees[comm - 1] + m.size();
    g2.nb_links += m.size();

    for (it = m.begin(); it != m.end(); ++it)
    {
      g2.total_weight += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }

  return g2;
}

bool Community::one_level() {
  bool improvement = false;
  int nb_moves;
  int nb_pass_done = 0;
  double new_mod = modularity();
  double cur_mod = new_mod;

  std::vector<int> random_order(m_Size);
  for (int i = 0; i < m_Size; i++)
  {
    random_order[i] = i;
  }
  for (int i = 0; i < m_Size - 1; i++)
  {
    int rand_pos = rand() % (m_Size - i) + i;
    int tmp = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  // repeat while
  //   there is an improvement of modularity
  //   or there is an improvement of modularity greater than a given epsilon
  //   or a predefined number of pass have been done
  do
  {
    cur_mod = new_mod;
    nb_moves = 0;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp = 0; node_tmp < m_Size; node_tmp++)
    {
      //      int node = node_tmp;
      int node = random_order[node_tmp];
      int node_comm = m_node2communityMap[node];
      double w_degree = m_graph.getWeightedDegree(node);

      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, m_neighborWeights[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm = node_comm;
      double best_nblinks = 0.;
      double best_increase = 0.;
      for (unsigned int i = 0; i < neigh_last; i++)
      {
        const double increase = modularity_gain(node, m_neighPositions[i], m_neighborWeights[m_neighPositions[i]], w_degree);
        if (increase > best_increase) {
          best_comm = m_neighPositions[i];
          best_nblinks = m_neighborWeights[m_neighPositions[i]];
          best_increase = increase;
        }
      }

      // insert node in the nearest community
      insert(node, best_comm, best_nblinks);

      if (best_comm != node_comm)
      {
        nb_moves++;
      }
    }

    double total_tot = 0;
    double total_in = 0;
    for (unsigned int i = 0; i < tot.size(); ++i)
    {
      total_tot += tot[i];
      total_in += in[i];
    }

    new_mod = modularity();
    if (nb_moves > 0)
    {
      improvement = true;
    }
  } while (
    nb_moves > 0
    && new_mod - cur_mod > m_minModularity);

  return improvement;
}
