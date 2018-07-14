#pragma once

#ifndef GRAPH_H
#define GRAPH_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <sstream>
#include <cctype>
#include <fstream>

#define lowbit(x) ((x)&(-(x)))
#define sqr(x) ((x)*(x))

struct Edge
{
  int i;
  int j;
  double w;
  Edge(int i = 0, int j = 0, double w = 0)
    : i(i)
    , j(j)
    , w(w)
  {
  }
};

class Graph
{
public:
  std::vector<Edge> m_edges;
  std::vector<std::vector<std::pair<int, double>>> links;
  std::vector<double> m_nodeWeights;
  double M;
  int N;	//how many nodes
  bool isCalculated;

  Graph(int N = 0)
    : N(N)
  {
    m_edges.clear();
    isCalculated = false;
  }

  bool load(const std::string& filename)
  {
    bool ispair = false;
    if (filename.substr(filename.size() - 4) == "pair")
    {
      ispair = true;
    }

    FILE* fin = fopen(filename.c_str(), "r");
    if (fin == NULL)
    {
      return false;
    }
    fscanf(fin, "%d", &N);
    m_edges.clear();
    int i, j;
    double w;
    std::set<std::pair<int, int>>pairs;
    pairs.clear();

    if (!ispair)
    {
      for (; fscanf(fin, "%d%d%lf", &i, &j, &w) == 3;)
      {
        if (pairs.find(std::make_pair(i, j)) != pairs.end())
        {
          continue;
        }
        //save edge (i,j) and (j,i)
        m_edges.push_back(Edge(i, j, w));
        m_edges.push_back(Edge(j, i, w));
        pairs.emplace(i, j);
        pairs.emplace(j, i);
      }
    }
    else
    {
      int cnt = 0;
      int THRES = 1000000;
      for (; fscanf(fin, "%d%d", &i, &j) == 2;) {
        cnt++; if (cnt % 10000 == 0) printf("\rLoading #%d", cnt);
        if (pairs.find(std::make_pair(i, j)) != pairs.end()) continue;
        if (i > THRES || j > THRES)
        {
          continue;
        }
        m_edges.push_back(Edge(i, j, 1));
        m_edges.push_back(Edge(j, i, 1));
        pairs.emplace(i, j);
        pairs.emplace(j, i);
      }
    }
    fclose(fin);
    isCalculated = false;
    return true;
  }

  void save(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    fprintf(fout, "%d\n", N);
    for (int i = 0; i < m_edges.size(); i++)
    {
      if (m_edges[i].i <= m_edges[i].j)
      {
        fprintf(fout, "%d\t%d\t%f\n", m_edges[i].i, m_edges[i].j, m_edges[i].w);
      }
    }
    fclose(fout);
  }

  void calculateLinks()
  {
    if (isCalculated)
    {
      return;
    }
    links.resize(N + 1);
    m_nodeWeights.resize(N + 1);
    for (int i = 0; i <= N; i++) {
      links[i].clear();
      m_nodeWeights[i] = 0;
    }
    M = 0;
    for (int i = 0; i < m_edges.size(); i++)
    {
      links[m_edges[i].i].push_back(std::make_pair(m_edges[i].j, m_edges[i].w));
      m_nodeWeights[m_edges[i].i] += m_edges[i].w;
      M += m_edges[i].w;
    }
    M /= 2;
    isCalculated = true;
  }

  void saveSingleEdge(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 0; i < m_edges.size(); i++)
      if (m_edges[i].i < m_edges[i].j)
      {
        fprintf(fout, "%d\t%d\t%f\n", m_edges[i].i, m_edges[i].j, m_edges[i].w);
      }
    fprintf(fout, "%d\t%d\t%d\n", N, N, 0);
    fclose(fout);
  }

  void savePairs(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 0; i < m_edges.size(); i++)
      if (m_edges[i].i < m_edges[i].j && getRand() < m_edges[i].w)
        fprintf(fout, "%d\t%d\n", m_edges[i].i, m_edges[i].j);
    //fprintf(fout,"%d\t%d\n",N,N);
    fclose(fout);
  }

  void savePairWeight(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 0; i < m_edges.size(); i++)
    {
      if (m_edges[i].i < m_edges[i].j
        && getRand() < m_edges[i].w)
      {
        fprintf(fout, "%d\t%d\t%f\n", m_edges[i].i, m_edges[i].j, m_edges[i].w);
      }
    }
    //fprintf(fout,"%d\t%d\n",N,N);
    fclose(fout);
  }

  void savePairsBiDir(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 0; i < m_edges.size(); i++)
    {
      if (getRand() < m_edges[i].w)
      {
        fprintf(fout, "%d\t%d\n", m_edges[i].i, m_edges[i].j);
      }
    }
    fprintf(fout, "%d\t%d\n", N, N);
    fclose(fout);
  }

  void addEdge(int i, int j, double w)
  {
#ifdef DEBUG
    std::cout << i << ", " << j << ", " << w << endl;
#endif
    m_edges.push_back(Edge(i, j, w));
    isCalculated = false;
  }

  void sampleByWeights()
  {
    std::vector<Edge> tg = m_edges;
    m_edges.clear();
    std::cout << "GRAPH SIZE : " << tg.size() << std::endl;
    for (int i = 0; i < tg.size(); i++)
    {
      if (getRand(0, 1) <= tg[i].w)
      {
        m_edges.push_back(Edge(tg[i].i, tg[i].j, 1));
      }
    }
  }

  void sub(std::map<int, int> rid)
  {
    N = rid.size();

    for (int z = 0; z < m_edges.size();)
    {
      int i = m_edges[z].i;
      int j = m_edges[z].j;
      if (rid.find(i) == rid.end()
        || rid.find(j) == rid.end())
      {
        std::swap(m_edges[z], m_edges[m_edges.size() - 1]);
        m_edges.pop_back();
      }
      else
      {
        m_edges[z] = Edge(rid[i], rid[j], m_edges[z].w);
        z++;
      }
    }
    M = m_edges.size();
  }

  std::vector<int> randomSample(int cnt)
  {
    std::vector<int> res;
    std::set<int> vi;
    vi.clear();
    calculateLinks();
    for (int cur = 0; res.size() < cnt;)
    {
      if (cur == res.size())
      {
        int k = rand() % N + 1;
        for (; vi.find(k) != vi.end(); k = rand() % N + 1);
        vi.insert(k);
        res.push_back(k);
      }
      else
      {
        int i = res[cur];
        cur++;
        for (int z = 0; z < links[i].size(); z++)
        {
          if (links[i][z].second > 0 && vi.find(links[i][z].first) == vi.end())
          {
            vi.insert(links[i][z].first);
            res.push_back(links[i][z].first);
          }
        }
      }
    }
    res.resize(cnt);
    return res;
  }

  void saveWalkTrap(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 0; i <= N; i++) {
      fprintf(fout, "%d %d 1\n", i, i);
    }
    for (int z = 0; z < m_edges.size(); z++)
    {
      if (m_edges[z].i >= m_edges[z].j)
      {
        continue;
      }
      fprintf(fout, "%d %d %.10f\n", m_edges[z].i, m_edges[z].j, m_edges[z].w);
    }
    fclose(fout);
  }
};

#endif // GRAPH_H
