// File: graph.cpp
// -- simple graph handling source file
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

#include "graph.h"

Graph::Graph(char *filename, int type)
{
  std::ifstream finput;
  finput.open(filename, std::fstream::in);

  int nb_links = 0;

  while (!finput.eof()) {
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
        links[dest].emplace_back(src, weight);

      nb_links++;
    }
  }

  finput.close();
}

void Graph::renumber(int type) {
  std::vector<int> linked(links.size(), -1);
  std::vector<int> renum(links.size(), -1);
  int nb = 0;

  for (unsigned int i = 0; i < links.size(); i++)
  {
    for (unsigned int j = 0; j < links[i].size(); j++)
    {
      linked[i] = 1;
      linked[links[i][j].first] = 1;
    }
  }

  for (unsigned int i = 0; i < links.size(); i++) {
    if (linked[i] == 1)
      renum[i] = nb++;
  }

  for (unsigned int i = 0; i < links.size(); i++) {
    if (linked[i] == 1) {
      for (unsigned int j = 0; j < links[i].size(); j++) {
        links[i][j].first = renum[links[i][j].first];
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

    for (auto& j : link)
    {
      it = m.find(j.first);
      if (it == m.end())
      {
        m.emplace(j.first, j.second);
      }
      else if (type == WEIGHTED)
      {
        it->second += j.second;
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
      int dest = links[i][j].first;
      float weight = links[i][j].second;
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

void Graph::display_binary(char *filename, char *filename_w, int type) {
  std::ofstream foutput;
  foutput.open(filename, std::fstream::out | std::fstream::binary);

  unsigned int s = links.size();

  // outputs number of nodes
  foutput.write((char *)(&s), 4);

  // outputs cumulative degree sequence
  long tot = 0;
  for (unsigned int i = 0; i < s; i++) {
    tot += (long)links[i].size();
    foutput.write((char *)(&tot), 8);
  }

  // outputs links
  for (auto& link : links) {
    for (auto& j : link)
    {
      int dest = j.first;
      foutput.write((char *)(&dest), 4);
    }
  }
  foutput.close();

  // outputs weights in a separate file
  if (type == WEIGHTED)
  {
    std::ofstream foutput_w;
    foutput_w.open(filename_w, std::fstream::out | std::fstream::binary);
    for (unsigned int i = 0; i < s; i++)
    {
      for (auto& j : links[i])
      {
        float weight = j.second;
        foutput_w.write((char *)(&weight), 4);
      }
    }
    foutput_w.close();
  }
}
