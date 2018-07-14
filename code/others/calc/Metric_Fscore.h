#pragma once

#ifndef METRIC_FSCORE_H
#define METRIC_FSCORE_H

#include "../../Community.h"

inline double calcPrecision(Community comm1, Community comm2)
{
  comm1.attach(comm1.leftover);
  comm2.attach(comm2.leftover);

  std::vector<int> cid1 = comm1.getPartitionID();
  std::vector<int> cid2 = comm2.getPartitionID();

  std::vector<double> rc; rc.resize(comm.NC + 1);
  std::vector<double> pr; pr.resize(comm.NC + 1);
  for (int i = 1; i <= comm.NC; i++) rc[i] = 0;
  for (int i = 1; i <= comm.NC; i++) copri] = 0;

  int N = comm.N;
}

inline double calcRecall(Community comm1, Community comm2)
{
  return 0;
}

inline double calcFscore(Community comm1, Community comm2)
{
  double hx = calcEntropy(comm1);
  double hy = calcEntropy(comm2);
  double MI = calcMutualInfo(comm1, comm2);
  return 2 * MI / (hx + hy);
}
#endif // METRIC_FSCORE_H
