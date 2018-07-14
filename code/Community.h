#pragma once

#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
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
  std::vector<int> cid;  // maps: node (i) -> community
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
  std::vector<int> id_size;
  std::vector<int> leftover;
  // Number of communities
  int NC;
  int N; //max node No.
  bool isPartition;	//is not overlapping

  bool testPartition()
  {
    std::set<int> t;
    t.clear();
    for (int i = 1; i <= NC; i++)
    {
      for (int j = 0; j < id[i].size(); j++)
      {
        if (t.find(id[i][j]) != t.end())
        {
          return false;
        }
        t.insert(id[i][j]);
      }
    }
    return true;
  }

  Community() = default;

  void sync_id_size()
  {
    id_size.resize(id.size());
    for (int i = 0; i < id.size(); i++)
    {
      id_size[i] = id[i].size();
    }
  }

  void load(const std::string& filename)
  {
    const std::string suffix = filename.substr(filename.size() - 3);
    if (suffix == "par")
    {
      loadPartition(filename);
    }
    else
    {
      loadGeneral(filename);
    }
  }

  void save(const std::string& filename)
  {
    const std::string suffix = filename.substr(filename.size() - 3);
    if (suffix == "par")
    {
      savePartition(filename);
    }
    else if (suffix == "gen")
    {
      saveGeneral(filename);
    }
    else
    {
      std::cout << "ERROR: save error: please check the suffix of the community file" << std::endl;
      exit(-1);//exit with error
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
    {
      fprintf(fout, "%d\n", cid[i]);
    }
    fclose(fout);
  }

  void loadGeneral(const std::string& filename)
  {
    FILE* fin = fopen(filename.c_str(), "r");
    std::map<int, int> tcid;	//i belongs to community tcid[i]
    tcid.clear();
    char s[1000000];
    std::vector<int> cur;
    cur.clear();
    id.clear();
    id.push_back(cur);	//id[1] is the first community. id[0] is unused
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

        //check if overlapping
        if (tcid.find(i) != tcid.end()) isPartition = false;
        tcid[i] = curid;
      }
      if (cur.size() > 0)
      {
        id.push_back(cur);
      }
    }
    fclose(fin);

    //not overlapping
    if (isPartition)
    {
      cid.clear();
      cid.push_back(0);
      for (int i = 1; i <= N; i++)
      {
        cid.push_back(tcid[i]);
      }
    }
    leftover.clear();
    NC = id.size() - 1;
  }

  void attach(std::vector<int> t)
  {
    if (t.size() != 0)
    {
      id.push_back(t);
      NC++;
    }
  }

  void pop()
  {
    id.pop_back();
    NC--;
  }

  void saveGeneral(const std::string& filename)
  {
    sync_id_size();
    //if (!testPartition()) std::cout<<"ERR"<<" "<<filename<< std::endl;
    FILE* fout = fopen(filename.c_str(), "w");
    if (config.find("Attach_Leftover") != config.end())
    {
      if (config["Attach_Leftover"] == "TRUE")
      {
        attach(leftover);
      }
    }
    std::vector<int> sorted_index = sort_indexes(id_size);
#ifdef DEBUG
    //assert(sorted_index.size()==id_size.size());
    //assert(comm.id.size()==id_size.size());
    std::cout << endl;
    for (int i = 0; i < id.size(); i++) {
      std::cout << sorted_index[i] << " : " << id_size[sorted_index[i]] << endl;
    }
#endif
    int sorted_i = 0;
    for (int index = 0; index <= NC; index++)
    {
      sorted_i = sorted_index[index];
      if (id[sorted_i].size())
      {
        for (int j = 0; j < id[sorted_i].size(); j++)
        {
          fprintf(fout, "%d\t", id[sorted_i][j]);
        }
        fprintf(fout, "\n");
      }
    }
    fclose(fout);
      }

  void saveSucheta(FILE* fout, std::vector<int> rid)
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
      if (k == 2)
      {
        break;
      }
      cid.push_back(j);
    }
    N = cid.size() - 1;
    for (int layer = 1;; layer++)
    {
      std::map<int, int> rid;
      for (int i, j; fscanf(fin, "%d%d", &i, &j) == 2;)
      {
        if (i == 0) break;
        rid[i] = j;
      }
      if (rid.empty())
      {
        break;
      }
      for (int i = 0; i < cid.size(); i++)
      {
        cid[i] = rid[cid[i]];
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
      if (cur.size() > 0)
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
#ifdef DEBUG
        std::cout << "id[" << i << "].size()=" << id[i].size() << "\tThres=" << Thres << endl;
#endif
        for (int z = 0; z < id[i].size(); z++)
        {
          leftover.push_back(id[i][z]);
        }
        id[i].clear();
        while ((id[NC].size() < Thres) && (NC > i))
        {
#ifdef DEBUG
          std::cout << "id[" << NC << "].size()=" << id[NC].size() << endl;
#endif
          NC--;
          id.pop_back();
        }
        id[i] = id[NC];
        NC--;
        id.pop_back();
        if (NC + 1 != id.size())
        {
#ifdef DEBUG
          fprintf(stderr, "Function removeSmallComms in Community.h error");
#endif
          exit(-1);
        }
      }
  }

  //return res, i belongs to community res[i]
  std::vector<int> getPartitionID()
  {
    std::vector<int> res;
    res.clear();
    //if (!isPartition) return res;
    res.resize(N + 1);
    for (int i = 0; i <= N; i++) res[i] = -1;
    for (int i = 1; i <= NC; i++)
      for (int j = 0; j < id[i].size(); j++)
        res[id[i][j]] = i;
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
        res[id[i][j]].push_back(i);
      }
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
      if (cur.size() > 0) id.push_back(cur);
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
    };
#endif // COMMUNITY_H
