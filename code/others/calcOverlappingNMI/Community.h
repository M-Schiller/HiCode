#pragma once

#ifndef COMMUNITY_H
#define COMMUNITY_H

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

class Community
{
  std::vector<int> cid;
  void calcid()
  {
    id.clear();
    for (int i = 1; i < cid.size(); i++)
    {
      for (; cid[i] >= id.size();)
      {
        std::vector<int> t;
        t.clear();
        id.push_back(t);
      }
      id[cid[i]].push_back(i);
    }
    NC = id.size() - 1;
  }
public:
  std::vector<std::vector<int>> id;
  std::vector<int> leftover;
  int NC, N;
  bool isPartition;

  bool testPartition()
  {
    std::set<int> t;
    t.clear();
    for (int i = 1; i <= NC; i++)
    {
      for (int z : id[i])
      {
        if (t.find(z) != t.end())
          return false;
        t.insert(z);
      }
    }
    return true;
  }

  Community() = default;

  void load(const std::string& filename)
  {
    if (filename.substr(filename.size() - 3) == "par")
      loadPartition(filename);
    else loadGeneral(filename);
  }

  void save(const std::string& filename)
  {
    if (filename.substr(filename.size() - 3) == "par")
    {
      savePartition(filename);
    }
    else
    {
      saveGeneral(filename);
    }
  }

  void loadPartition(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    cid.clear();
    cid.push_back(0);
    for (int i; fscanf(fin, "%d", &i) == 1;)
      cid.push_back(i);
    fclose(fin);
    N = cid.size() - 1;
    isPartition = true;

    reId();
    leftover.clear();
    calcid();
  }

  void savePartition(const std::string& filename)
  {
    FILE* fout = fopen(filename.c_str(), "w");
    for (int i = 1; i < cid.size(); i++)
      fprintf(fout, "%d\n", cid[i]);
    fclose(fout);
  }

  void loadGeneral(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    std::map<int, int> tcid;
    tcid.clear();
    char s[1000000];
    std::vector<int> cur;
    cur.clear();
    id.clear();
    id.push_back(cur);
    isPartition = true;
    int curid = 0;
    N = 0;
    for (; fgets(s, 1000000, fin);)
    {
      curid++;
      cur.clear();
      for (char* p = strtok(s, "\t \n\r"); p; p = strtok(NULL, "\t \n\r"))
      {
        int i = std::atoi(p);
        cur.push_back(i);
        N = std::max(N, i);
        if (tcid.find(i) != tcid.end()) isPartition = false;
        tcid[i] = curid;
      }
      if (!cur.empty())
      {
        id.push_back(cur);
      }
    }
    fclose(fin);
    if (isPartition)
    {
      cid.clear();
      cid.push_back(0);
      for (int i = 1; i <= N; i++)
        cid.push_back(tcid[i]);
    }
    leftover.clear();
    NC = id.size() - 1;
  }

  void attach(const std::vector<int>& t)
  {
    id.push_back(t);
    NC++;
  }

  void pop()
  {
    id.pop_back();
    NC--;
  }

  void saveGeneral(const std::string& filename)
  {
    //if (!testPartition()) std::cout<<"ERR"<<" "<<filename<< std::endl;
    FILE* fout = fopen(filename.c_str(), "w");
    attach(leftover);
    for (int i = 1; i <= NC; i++)
    {
      if (!id[i].empty())
      {
        for (int j : id[i])
        {
          fprintf(fout, "%d\t", j);
        }
        fprintf(fout, "\n");
      }
    }
    pop();
    //for (int i=0;i<leftover.size();i++) fprintf(fout,"%d\n",leftover[i]);
    fclose(fout);
  }

  void saveSucheta(FILE* fout, const std::vector<int>& rid)
  {
    /*calcid();
    for (int i=1;i<id.size();i++){
      for (int j=0;j<id[i].size();j++)
        fprintf(fout,"%d\t",rid[id[i][j]]);
      if (id[i].size())fprintf(fout,"\n");
    }*/
  }

  void loadModularity(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    cid.clear();
    for (int i, j, k = 0; fscanf(fin, "%d%d", &i, &j) == 2;)
    {
      k += (i == 0);
      if (k == 2) break;
      cid.push_back(j);
    }
    N = cid.size() - 1;
    for (int layer = 1;; layer++)
    {
      std::map<int, int> rid;
      rid.clear();
      for (int i, j; fscanf(fin, "%d%d", &i, &j) == 2;)
      {
        if (i == 0) break;
        rid[i] = j;
      }
      if (rid.empty()) break;
      for (int& id : cid)
      {
        id = rid[id];
      }
    }
    fclose(fin);
    leftover.clear();
    reId();
    calcid();
  }

  void loadInfomap(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    cid.clear();
    cid.push_back(0);
    fscanf(fin, "%*s%*d");
    for (int i; fscanf(fin, "%d", &i) == 1;)
      cid.push_back(i);
    fclose(fin);
    N = cid.size() - 1;
    isPartition = true;

    reId();
    leftover.clear();
    calcid();
  }

  void loadOSLOM(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    std::map<int, int> tcid;
    tcid.clear();
    char s[1000000];
    std::vector<int> cur;
    cur.clear();
    id.clear();
    id.push_back(cur);
    isPartition = true;
    int curid = 0;
    N = 0;
    int linenumber = 0;
    for (; fgets(s, 1000000, fin);)
    {
      linenumber++;
      if (linenumber % 2) continue;
      curid++;
      cur.clear();
      for (char* p = strtok(s, "\t \n\r"); p; p = strtok(NULL, "\t \n\r"))
      {
        int i = std::atoi(p);
        cur.push_back(i);
        N = std::max(N, i);
        if (tcid.find(i) != tcid.end()) isPartition = false;
        tcid[i] = curid;
      }
      if (!cur.empty()) id.push_back(cur);
    }
    fclose(fin);
    if (isPartition)
    {
      cid.clear();
      cid.push_back(0);
      for (int i = 1; i <= N; i++)
        cid.push_back(tcid[i]);
    }
    leftover.clear();
    NC = id.size() - 1;
  }

  void reId()
  {
    std::map<int, int> rid;
    rid.clear();
    for (int i = 1, j = 1; i <= N; i++)
    {
      if (rid.find(cid[i]) == rid.end())
        rid[cid[i]] = j++;
      cid[i] = rid[cid[i]];
    }
  }

  void removeSmallComms(int Thres)
  {
    for (int i = 1; i <= NC; i++)
      if (id[i].size() < Thres)
      {
        for (int z : id[i])
        {
          leftover.push_back(z);
        }
        id[i].clear();
        id[i] = id[NC];
        NC--;
        id.pop_back();
        if (NC + 1 != id.size())
        {
          fprintf(stderr, "Function removeSmallComms in Community.h error");
          exit(-1);
        }
      }
  }
  std::vector<int> getPartitionID()
  {
    std::vector<int> res;
    res.clear();
    //if (!isPartition) return res;
    res.resize(N + 1);
    for (int i = 0; i <= N; i++)
    {
      res[i] = -1;
    }
    for (int i = 1; i <= NC; i++)
    {
      for (int j = 0; j < id[i].size(); j++)
      {
        res[id[i][j]] = i;
      }
    }
    //int cur=NC;
    //for (int i=0;i<=N;i++) if (res[i]==-1) res[i]=++NC;
    return res;
  }

  std::vector<std::vector<int>> getPartitionIDOver()
  {
    std::vector<std::vector<int>> res;
    for (int i = 0; i <= N; i++)
    {
      std::vector<int> tmp;
      res.push_back(tmp);
    }
    //cout<<"aaaaaaaaaa"<<res.size()<<" NC:"<<NC<< std::endl;
    //cout<<id.size()<<" "<<id[1].size()<<" "<<id[2].size()<<" "<<id[1][0]<< std::endl;
    //res.resize(N+1);
    //for (int i=0;i<=N;i++) res[i]=-1;
    /*for (int i=1;i<=NC;i++){
      for (int j=0;j<id[i].size();j++){
        std::cout<<id[i][j]<<" ";
      }
      std::cout<< std::endl;
    }*/
    for (int i = 1; i <= NC; i++)
    {
      for (int j = 0; j < id[i].size(); j++)
      {
        //cout<<"s1 ";
        //cout<<id[i][j]<<" ";
        res[id[i][j]].push_back(i);
      }
      //cout<< std::endl;
    }
    //cout<<"ssssssggggggg"<< std::endl;
    //int cur=NC;
    //for (int i=0;i<=N;i++) if (res[i]==-1) res[i]=++NC;
    return res;
  }

  void loadLinkCommunity(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    std::map<int, int> tcid;
    tcid.clear();
    char s[1000000];
    std::vector<int> cur;
    cur.clear();
    id.clear();
    id.push_back(cur);
    isPartition = true;
    int curid = 0;
    N = 0;
    for (; fgets(s, 1000000, fin);)
    {
      curid++;
      cur.clear();
      std::set<int> tmp;
      tmp.clear();
      for (char* p = strtok(s, "\t \n\r,"); p; p = strtok(NULL, "\t \n\r,"))
      {
        int i = std::atoi(p);
        if (tmp.find(i) != tmp.end()) continue;
        tmp.insert(i);
        cur.push_back(i);
        N = std::max(N, i);
        if (tcid.find(i) != tcid.end()) isPartition = false;
        tcid[i] = curid;
      }
      if (!cur.empty()) id.push_back(cur);
    }
    fclose(fin);
    if (isPartition)
    {
      cid.clear();
      cid.push_back(0);
      for (int i = 1; i <= N; i++)
        cid.push_back(tcid[i]);
    }
    leftover.clear();
    NC = id.size() - 1;
  }

  void sub(std::map<int, int> rid)
  {
    for (int i = 1; i <= NC; i++)
      for (int j = 0; j < id[i].size();)
      {
        if (rid.find(id[i][j]) == rid.end())
        {
          std::swap(id[i][j], id[i][id[i].size() - 1]);
          id[i].pop_back();
        }
        else
        {
          id[i][j] = rid[id[i][j]];
          j++;
        }
      }
    N = rid.size();
  }

  void loadWalkTrap(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    std::map<int, int> tcid;
    tcid.clear();
    char s[1000000];
    std::vector<int> cur;
    cur.clear();
    id.clear();
    id.push_back(cur);
    isPartition = true;
    int curid = 0;
    N = 0;
    bool st = false;
    fgets(s, 1000000, fin);
    if (sscanf(s, "Partition 0 (%d communities)", &N) != 1)
    {
      std::cout << s << std::endl;
      for (;;);
    }
    N--;
    for (; fgets(s, 1000000, fin);)
    {
      if (s[0] == 'M')
      {
        st = true;
        continue;
      }
      if (!st) continue;
      curid++;
      cur.clear();
      int k = 0;
      for (char* p = strtok(s, "{},\t \n\r"); p; p = strtok(NULL, "\t \n\r{},"))
      {
        k++;
        if (k > 3)
        {
          int i = std::atoi(p);
          if (i == 0)
          {
            continue;
          }
          cur.push_back(i);
          N = std::max(N, i);
          if (tcid.find(i) != tcid.end())
          {
            isPartition = false;
          }
          tcid[i] = curid;
        }
      }
      if (!cur.empty()) id.push_back(cur);
    }
    for (int i = 1; i <= N; i++) if (tcid.find(i) == tcid.end())
    {
      curid++;
      cur.clear();
      cur.push_back(i);
      tcid[i] = curid;
      id.push_back(cur);
    }
    fclose(fin);
    if (isPartition)
    {
      cid.clear();
      cid.push_back(0);
      for (int i = 1; i <= N; i++)
        cid.push_back(tcid[i]);
    }
    leftover.clear();
    NC = id.size() - 1;
  }
};
#endif // COMMUNITY_H
