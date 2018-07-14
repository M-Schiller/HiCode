#pragma once

#ifndef METRIC_NMI_H
#define METRIC_NMI_H

#include <iostream>
#include <map>
#include <vector>

inline double calcEntropy(Community comm)
{
  double E = 0;
  for (int i = 1; i <= comm.NC; i++)
  {
    double p = comm.id[i].size() / (double)(comm.N);
    if (p > 0) E -= p * log(p);
  }
  return E;
}

inline double calcMutualInfo(Community comm1, Community comm2)
{
  double I = 0;
  if (!comm1.leftover.empty())
  {
    //printf("%d %d", comm1.NC, comm1.leftover.size());
    comm1.attach(comm1.leftover);//为什么要把删掉的小社团加回来？
    //printf("%d\n", comm1.NC);
    comm2.attach(comm2.leftover);
  }

  std::vector<int> cid1 = comm1.getPartitionID();
  std::vector<int> cid2 = comm2.getPartitionID();
  std::vector<std::map<int, int>> cnt;
  cnt.resize(comm1.NC + 1);//维护两个commsX,Y中，x∈X,y∈Y，x ∩ y
  double res = log(0);

  if (!cid1.empty()
    && !cid2.empty())
  {
    for (int i = 1; i <= comm1.NC; i++)
    {
      cnt[i].clear();
    }
    int N = comm1.N;
    for (int i = 1; i <= N; i++)
    {
      if (cid1[i] < 0
        || cid2[i] < 0)
      {
        std::cout << "The " << i << "st vertex is not in any community from Comms1 or Comms2" << std::endl;
      }
      if (cid1[i] != -1
        && cid2[i] != -1)
      {
        cnt[cid1[i]][cid2[i]] += 1;
      }
      //如果第i个点存在于comm1和comm2中某个社团，出现再comm1中的第cid1[i]号社团中，也出现在comm2中的cid2[2]号社团
      //那么cid1[i]和cid2[i]的交集就加一
    }

    for (int i = 1; i <= comm1.NC; i++)
    {
      double pi = comm1.id[i].size() / (double)(N);
      for (std::map<int, int>::iterator itr = cnt[i].begin(); itr != cnt[i].end(); ++itr)
      {
        int j = itr->first;;
        double pj = comm2.id[j].size() / (double)(N);
        double pij = itr->second / (double)(N);
        I += pij * log(pij / pi / pj);
      }
    }
    res = I;
  }
  else
  {
    puts("error: \"cid1.size()!=0&&cid2.size()!=0\" not satisfied");
    getchar();
  }
  if (!comm1.leftover.empty())
  {
    comm1.pop();
    comm2.pop();
  }
  return res;
}

inline double calcNMI(Community comm1, Community comm2)
{
  double hx = calcEntropy(comm1);
  double hy = calcEntropy(comm2);
  double MI = calcMutualInfo(comm1, comm2);
  return 2 * MI / (hx + hy);
}

inline double calcNMI(const std::string& commFile1, const std::string& commFile2)
{
  std::cout << commFile1 << ' ' << commFile2 << std::endl;
  Community comm1;
  comm1.load(commFile1);
  Community comm2;
  comm2.load(commFile2);
  if (config.find("Metric_CommunitySizeThres") != config.end())
  {
    comm1.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
    comm2.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
  }

  return  calcNMI(comm1, comm2);
}
#endif // METRIC_NMI_H
