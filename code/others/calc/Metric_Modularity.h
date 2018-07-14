#pragma once

#ifndef METRIC_MODULARITY_H
#define METRIC_MODULARITY_H

#include <iostream>
#include "../../Graph.h"
#include "../../Community.h"

inline double calcModularity_real(Graph g, Community comm)
{
  double Q = 0;
  std::vector<double> e;
  e.clear();
  std::vector<double> a;
  a.clear();
  double M = 0;
  if (!comm.leftover.empty())
  {
    comm.attach(comm.leftover);
  }
  int nc = comm.NC;
  std::cout << "N = " << comm.N << ", NC = " << nc << std::end;
  std::vector<int> cid = comm.getPartitionID();

  for (int i = 0; i <= nc; i++)
  {
    e.push_back(0);
    a.push_back(0);
  }

  if (cid.empty())
  {
    if (!comm.leftover.empty())
      comm.pop();
    return log(0);
  }

  for (auto&[i, j, w] : g.m_edges)
  {
    //std::cout << "i = " << i << ", j = " << j << ", w = " << w << std::endl;
    if (cid[i] == cid[j]
      && cid[i] != -1)
    {
      e[cid[i]] += w;//e_in_in内部边
    }
    if (cid[i] != -1)
    {
      a[cid[i]] += w;//
    }
    M += w;
  }
  for (int i = 1; i <= nc; i++)
  {
    Q += e[i] / M - sqr((a[i]) / (M));
    //std::cout << "e[" << i << "] = " << e[i] << ", a[" << i << "] = " << a[i] << std::endl;
  }
  if (comm.leftover.size() != 0)
    comm.pop();
  /*for(int i=0; i<=cid.size();i++)
    std::cout<<cid[i]<<" ";
  std::cout<<"here"<< std::endl;
  std::cout<<"calcModularity Graph and comm"<< std::endl;*/
  return Q;
}

inline double calcModularity(const std::string& graphFile, const std::string& commFile)
{
  Graph g;
  if (!g.load(graphFile)) return log(0);
  Community comm{};
  comm.load(commFile);
  if (config.find("Metric_CommunitySizeThres") != config.end())
    comm.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
  return calcModularity_real(g, comm);
}

#endif // METRIC_MODULARITY_H
