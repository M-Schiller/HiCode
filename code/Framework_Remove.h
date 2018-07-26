#pragma once

#ifndef FRAMEWORK_REMOVE_H
#define FRAMEWORK_REMOVE_H

#include <vector>

class Framework_Remove : public Framework
{
public:
  //ȥ�������ڲ��ı�
  Graph calcNextLayerGraph(Graph cur, Community comm) override
  {
    auto g = Graph(cur.N);
    std::cout << "sfsjflksjdfs" << std::endl;
    std::vector<int> cid = comm.getPartitionID();

    //����ͼ��ÿ����
    for (int z = 0; z < cur.m_edges.size(); z++)
    {
      int i = cur.m_edges[z].i;
      int j = cur.m_edges[z].j;
      if (cid[i] == cid[j]) continue;	//i j��ͬһ������ʱ ���Ըñ�
      g.addEdge(i, j, cur.m_edges[z].w);	//i j����ͬһ������ �ӵ������ͼ��
    }
    return g;
  }
};
#endif // FRAMEWORK_REMOVE_H
