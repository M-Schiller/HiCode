#pragma once

#ifndef FRAMEWORK_REDUCEPP_H
#define FRAMEWORK_REDUCEPP_H

#include <vector>

class Framework_ReducePP : public Framework
{//overlapping reduce++
public:
  virtual Graph calcNextLayerGraph(Graph cur, Community comm)
  {
    auto g = Graph(cur.N);
    double sum = 0;
    std::vector<std::vector<int>> cid = comm.getPartitionIDOver();
    std::vector<double> cs(comm.NC + 1);
    std::vector<double> co(comm.NC + 1);

    for (int i = 1; i <= comm.NC; i++)
    {
      cs[i] = 0;
    }
#ifdef DEBUG
    std::cout << "G_SIZE" << cur.g.size() - 1 << endl;
#endif

    for (auto &[i, j, w] : cur.m_edges)
    {
#ifdef DEBUG
      cout << "here z : " << z << std::endl;
#endif
      sum += w;
      int flag = 0;
      int cid_i_k = 0, cid_j_s = 0;
      if (i <= comm.N && j <= comm.N)
      {
        for (int k = 0; k < cid[i].size(); k++)
        {
          cid_i_k = cid[i][k];
          if (cid_i_k > 0)
          {
            co[cid_i_k] += w;
          }
          for (int s = 0; s < cid[j].size(); s++)
          {
            cid_j_s = cid[j][s];

            if (cid_i_k == cid_j_s
              && cid_i_k > 0)
            {
              flag = 1;
              cs[cid_i_k] += w;
              break;
            }
          }
          //if(flag) break;
        }
      }

      if (!flag)
      {
        g.addEdge(i, j, w);
      }
    }

    for (int edge_it = 0; edge_it < cur.m_edges.size(); edge_it++)
    {
      int i = cur.m_edges[edge_it].i;
      int j = cur.m_edges[edge_it].j;
      int flag = 0;
      int k = 0;
      int s = 0;
      int cid_i_k = 0;
      if (i <= comm.N
        && j <= comm.N)
      {
        for (; k < cid[i].size(); k++)
        {
          for (; s < cid[j].size(); s++)
          {
            cid_i_k = cid[i][k];
            int cid_j_s = cid[j][s];
            if (cid_i_k == cid_j_s
              && cid_i_k > 0)
            {
              flag = 1;
              break;
            }
          }
          if (flag)
          {
            break;
          }
        }
      }

      if (flag)
      {
        int tcid = cid_i_k;
        int tc_size = comm.id[tcid].size();
        int tc_size_sqr = tc_size * tc_size;
        double pi = cs[tcid] / tc_size_sqr;
        double qi = (co[tcid] - cs[tcid]) / (g.N * tc_size - tc_size_sqr);
        if (pi > 0)
        {
          g.addEdge(i, j, cur.m_edges[edge_it].w*qi / pi);
        }
      }
    }
#ifdef DEBUG
    std::cout << "reduce end" << endl;
#endif
    return g;
  }
};
#endif // FRAMEWORK_REDUCEPP_H
