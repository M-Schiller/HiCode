#pragma once

#ifndef FRAMEWORK_REMOVE_H
#define FRAMEWORK_REMOVE_H

#include <vector>

class Framework_Remove : public Framework
{
public:
  //去掉社团内部的边
  Graph calcNextLayerGraph(Graph cur, Community comm) override
  {
    auto g = Graph(cur.N);
    std::cout << "sfsjflksjdfs" << std::endl;
    std::vector<int> cid = comm.getPartitionID();

    //遍历图的每条边
    for (int z = 0; z < cur.m_edges.size(); z++)
    {
      int i = cur.m_edges[z].i;
      int j = cur.m_edges[z].j;
      if (cid[i] == cid[j]) continue;	//i j在同一个社团时 忽略该边
      g.addEdge(i, j, cur.m_edges[z].w);	//i j不在同一个社团 加到结果的图中
    }
    return g;
  }
};
#endif // FRAMEWORK_REMOVE_H
