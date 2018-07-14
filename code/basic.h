#pragma once

#ifndef BASIC_H
#define BASIC_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <numeric>
#include <sstream>
#include <cctype>
#include <fstream>

#define lowbit(x) ((x)&(-(x)))
#define sqr(x) ((x)*(x))

template <typename T>
std::vector<int> sort_indexes(const std::vector<T> &v) {
  // initialize original index locations
  std::vector<int> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  if (config.find("ReduceFirst") != config.end())
  {
    if (config["ReduceFirst"] == "MIN")
    {
      sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });
    }
    else if (config["ReduceFirst"] == "MAX")
    {
      sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2]; });
    }
    else
    {
      std::cout << "error parameter with ReduceFirst" << std::endl;
      exit(-1);
    }
  }
  else
  {
    std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2]; });
  }

  return idx;
}

inline bool checkRequiredFile(const std::string& filename)
{
  std::ifstream f(filename, std::fstream::in | std::fstream::trunc);
  if (!f.good())
  {
    std::cerr << "FILE " << filename << " DOES NOT EXIST AND FAILED TO CREATE IT" << std::endl;
  }
  return f.good();
}

inline bool checkFileExist(const std::string& filename)
{
  std::ifstream f(filename);
  return f.good();
}

inline std::string int2str(int i)
{
  return std::to_string(i);
}

inline std::string double2str(double v)
{
  char s[100];
  std::snprintf(s, sizeof s, "%.6f", v);
  return s;
}

inline double str2double(const std::string& s)
{
  return std::stod(s);
}

inline void systemCall(const std::string& cmd)
{
  std::cerr << "Call " << cmd;
  int res = system(cmd.c_str());
  std::cerr << " returned " << res << "." << std::endl;
}

//返回l到r之间的随机数
inline double getRand(double l = 0, double r = 1)
{
  return (rand() % 1000000 / 1000000.) * (r - l) + l;
}

inline std::string escape(const std::string& s)
{
  std::string t;
  for (int i = 0; i < s.size(); i++)
  {
    if (s[i] == '\\'
      || s[i] == '_'
      || s[i] == '#')
    {
      t = t + "\\";
    }
    t = t + s[i];
  }
  return t;
}

inline void mergeMSS(std::map<std::string, std::string> &a, std::map<std::string, std::string> b)
{
  for (auto& itr : b)
  {
    a[itr.first] = itr.second;
  }
}

inline void mergeVSS(std::vector<std::pair<std::string, std::string>>&a, std::vector<std::pair<std::string, std::string>> b)
{
  a.insert(a.end(), b.begin(), b.end());
}

inline double getAverage(std::vector<double> d)
{
  return std::accumulate(d.begin(), d.end(), 0.0) / d.size();
}

inline std::vector<std::vector<std::string>> VSS2VVS(std::vector<std::pair<std::string, std::string>> data)
{
  std::vector<std::vector<std::string>> table;
  for (auto& i : data)
  {
    std::vector<std::string> row;
    row.push_back(i.first);
    row.push_back(i.second);
    table.push_back(row);
  }
  return table;
}

inline std::vector<std::vector<std::string>> transVVS(std::vector<std::vector<std::string>> data)
{
  std::vector<std::vector<std::string>> newData;
  for (int i = 0; i < data[0].size(); i++)
  {
    std::vector<std::string> row;
    for (int j = 0; j < data.size(); j++)
    {
      row.push_back(data[j][i]);
    }
    newData.push_back(row);
  }
  return newData;
}

inline void Tex_Table(std::vector<std::vector<std::string>> table, const std::string& name, FILE* fout)
{
  if (table.empty())
    return;
  fprintf(fout, "\\begin{table}[htb!]\n\\centering\n");
  fprintf(fout, "\\caption{%s}\n", escape(std::move(name)).c_str());
  fprintf(fout, "\\begin{tabular}{|");
  for (int i = 0; i < table[0].size(); i++)
    fprintf(fout, "c|");
  fprintf(fout, "}\n");
  for (int i = 0; i < table.size(); i++)
  {
    fprintf(fout, "\\hline\n");
    for (int j = 0; j < table[i].size(); j++)
    {
      if (j) fprintf(fout, "&");
      fprintf(fout, "%s", escape(table[i][j]).c_str());
    }
    fprintf(fout, "\\\\\n");
  }
  fprintf(fout, "\\hline\n");
  fprintf(fout, "\\end{tabular}\\end{table}\n");
}

inline void Tex_includeGraphics(const std::string& filename, const std::string& name, double width, FILE* fout)
{
  fprintf(fout, "\\begin{figure}[htb!]\n\\centering\n");
  fprintf(fout, "\\caption{%s}\n", escape(std::move(name)).c_str());
  fprintf(fout, "\\includegraphics[width=%.1fin]{%s}\n", width, filename.c_str());
  fprintf(fout, "\\end{figure}\n");
}

inline std::vector<std::string> split(const std::string& str)
{
  std::vector<std::string> res;
  res.clear();
  char s[10000];
  sprintf(s, "%s", str.c_str());
  for (char* p = strtok(s, ";"); p; p = strtok(NULL, ";"))
    res.push_back(p);
  return res;
}

inline std::vector<int> intersect(std::vector<int> Xi, std::vector<int> Yj)
{
  int Xi_size = Xi.size();
  int Yj_size = Yj.size();
  std::vector<int> v(std::max(Xi_size, Yj_size));
  if (Xi_size == 0)
    return Xi;
  if (Yj_size == 0)
    return Yj;

  int Xi_array[Xi_size];
  int Yj_array[Yj_size];

  for (int i = 0; i < Xi_size; i++)
  {
    Xi_array[i] = Xi[i];
  }

  for (int j = 0; j < Yj_size; j++)
  {
    Yj_array[j] = Yj[j];
  }

  int* pXi = Xi_array;
  int* pYj = Yj_array;
  std::sort(pXi, pXi + Xi_size);
  std::sort(pYj, pYj + Yj_size);
  std::vector<int>::iterator it = set_intersection(pXi, pXi + Xi_size, pYj, pYj + Yj_size, v.begin());
  v.resize(it - v.begin());
  return v;
}

inline int intersect_size(std::vector<int> Xi, std::vector<int> Yj)
{
  std::vector<int> v = intersect(std::move(Xi), std::move(Yj));
  return v.size();
}

inline std::vector<int> symmetric_difference(std::vector<int> Xi, std::vector<int> Yj)
{
  int Xi_size = Xi.size();
  int Yj_size = Yj.size();

  if (Yj_size == 0)
  {
    return Xi;
  }

  if (Xi_size == 0)
  {
    return Yj;
  }

  std::vector<int> v(std::max(Xi_size, Yj_size));

  int Xi_array[Xi_size];
  int Yj_array[Yj_size];

  for (int i = 0; i < Xi_size; i++)
  {
    Xi_array[i] = Xi[i];
  }
  for (int j = 0; j < Yj_size; j++)
  {
    Yj_array[j] = Yj[j];
  }

  int* pXi = Xi_array;
  int* pYj = Yj_array;
  std::sort(pXi, pXi + Xi_size);
  std::sort(pYj, pYj + Yj_size);
  const auto it = set_symmetric_difference(pXi, pXi + Xi_size, pYj, pYj + Yj_size, v.begin());
  v.resize(it - v.begin());
  return v;
}

inline std::vector<int> difference(std::vector<int> Xi, std::vector<int> Yj)
{
  int Xi_size = Xi.size();
  int Yj_size = Yj.size();

  if (Xi_size == 0
    || Yj_size == 0)
  {
    return Xi;
  }

  std::vector<int> v(std::max(Xi_size, Yj_size));

  int Xi_array[Xi_size];
  int Yj_array[Yj_size];
  for (int i = 0; i < Xi_size; i++)
  {
    Xi_array[i] = Xi[i];
  }

  for (int j = 0; j < Yj_size; j++)
  {
    Yj_array[j] = Yj[j];
  }

  int* pXi = Xi_array;
  int* pYj = Yj_array;
  std::sort(pXi, pXi + Xi_size);
  std::sort(pYj, pYj + Yj_size);
  const auto it = set_difference(pXi, pXi + Xi_size, pYj, pYj + Yj_size, v.begin());
  v.resize(it - v.begin());
  return v;
}
#endif // BASIC_H
