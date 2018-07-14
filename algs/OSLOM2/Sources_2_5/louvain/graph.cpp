// File: graph.cpp
// -- simple graph handling source file
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

#include "graph.h"
#include <algorithm>

Graph::Graph(char *filename, int type)
{
  std::ifstream finput;
  finput.open(filename, std::fstream::in);

  int nb_links = 0;

  while (!finput.eof())
  {
    unsigned int src, dest;
    double weight = 1.;

    if (type == WEIGHTED)
    {
      finput >> src >> dest >> weight;
    }
    else
    {
      finput >> src >> dest;
    }

    if (finput)
    {
      if (links.size() <= std::max(src, dest) + 1)
      {
        links.resize(std::max(src, dest) + 1);
      }

      links[src].emplace_back(dest, weight);
      if (src != dest)
      {
        links[dest].emplace_back(src, weight);
      }

      nb_links++;
    }
  }

  finput.close();
}

void Graph::renumber(int type)
{
  std::vector<int> linked(links.size(), -1);
  std::vector<int> renum(links.size(), -1);
  int nb = 0;

  for (unsigned int i = 0; i < links.size(); ++i)
  {
    for (unsigned int j = 0; j < links[i].size(); ++j)
    {
      linked[i] = 1;
      linked[links[i][j].first] = 1;
    }
  }

  for (unsigned int i = 0; i < links.size(); i++) {
    if (linked[i] == 1)
      renum[i] = nb++;
  }

  for (unsigned int i = 0; i < links.size(); ++i)
  {
    if (linked[i] == 1)
    {
      for (auto& j : links[i])
      {
        j.first = renum[j.first];
      }
      links[renum[i]] = links[i];
    }
  }
  links.resize(nb);
}

void Graph::clean(int type)
{
  for (auto& link : links)
  {
    std::map<int, float> m;
    std::map<int, float>::iterator it;

    for (auto& i : link)
    {
      it = m.find(i.first);
      if (it == m.end())
      {
        m.insert(i);
      }
      else if (type == WEIGHTED)
      {
        it->second += i.second;
      }
    }

    std::vector<std::pair<int, float>> v;
    for (it = m.begin(); it != m.end(); ++it)
    {
      v.emplace_back(*it);
    }
    link.clear();
    link = v;
  }
}

void Graph::display(int type)
{
  for (unsigned int i = 0; i < links.size(); ++i)
  {
    for (unsigned int j = 0; j < links[i].size(); ++j)
    {
      const int dest = links[i][j].first;
      const float weight = links[i][j].second;
      if (type == WEIGHTED)
      {
        std::cout << i << " " << dest << " " << weight << std::endl;
      }
      else
      {
        std::cout << i << " " << dest << std::endl;
      }
    }
  }
}

void Graph::display_binary(char *filename, char *filename_w, int type)
{
  std::ofstream foutput;
  foutput.open(filename, std::fstream::out | std::fstream::binary);

  unsigned int s = links.size();

  // outputs number of nodes
  foutput.write(reinterpret_cast<char *>(&s), 4);

  // outputs cumulative degree sequence
  long tot = 0;
  for (unsigned int i = 0; i < s; i++)
  {
    tot += static_cast<long>(links[i].size());
    foutput.write(reinterpret_cast<char *>(&tot), 8);
  }

  // outputs links
  for (auto & link : links)
  {
    for (auto& j : link)
    {
      foutput.write(reinterpret_cast<char *>(&j.first), 4);
    }
  }
  foutput.close();

  // outputs weights in a separate file
  if (type == WEIGHTED) {
    std::ofstream foutput_w;
    foutput_w.open(filename_w, std::fstream::out | std::fstream::binary);
    for (auto& link : links)
    {
      for (auto& j : link)
      {
        foutput_w.write(reinterpret_cast<char *>(&j.second), 4);
      }
    }
    foutput_w.close();
  }
}
