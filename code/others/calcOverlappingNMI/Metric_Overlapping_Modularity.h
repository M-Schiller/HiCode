#pragma once

#ifndef METRIC_OVERLAPPING_MODULARITY_H
#define METRIC_OVERLAPPING_MODULARITY_H

#include "../../Graph.h"
#include "../../Community.h"

inline double calcModularity(Graph g, Community comm)
{
  double Q = 0;
  std::vector<double> e; e.clear();
  std::vector<double> a; a.clear();
  double M = 0;
  comm.attach(comm.leftover);
  const int nc = comm.NC;
  std::vector<int> cid = comm.getPartitionID();

  for (int i = 0; i <= nc; i++)
  {
    e.push_back(0);
    a.push_back(0);
  }

  if (cid.empty())
  {
    comm.pop();
    return log(0);
  }

  for (auto &[i, j, w] : g.m_edges)
  {
    if (cid[i] == cid[j]
      && cid[i] != -1)
    {
      e[cid[i]] += w;
    }
    if (cid[i] != -1)
    {
      a[cid[i]] += w;
    }
    M += w;
  }
  for (int i = 1; i <= nc; i++)
  {
    Q += e[i] / M - a[i] / M * a[i] / M;
  }
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
  if (!g.load(graphFile))
  {
    return log(0);
  }
  Community comm;
  comm.load(commFile);
  if (config.find("Metric_CommunitySizeThres") != config.end())
  {
    comm.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
  }
  return calcModularity(g, comm);
}

#endif // METRIC_OVERLAPPING_MODULARITY_H
