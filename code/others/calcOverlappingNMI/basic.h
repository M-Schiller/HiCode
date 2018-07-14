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
#include <sstream>
#include <cctype>
#include <fstream>
#include <numeric>

#define lowbit(x) ((x)&(-(x)))
#define sqr(x) ((x)*(x))

inline bool checkRequiredFile(const std::string& filename)
{
  FILE* fin = fopen(filename.c_str(), "r");
  if (fin == NULL)
  {
    systemCall(("make " + filename).c_str());
    fin = fopen(filename.c_str(), "r");
    if (fin == NULL)
    {
      fprintf(stderr, "FILE %s NOT EXIST, AND FAIL TO MAKE\n", filename.c_str());
      return false;
    }
  }
  fclose(fin);
  return true;
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
  std::cerr << "Call: " << cmd;
  const int res = system(cmd.c_str());
  std::cerr << " returned " << res << "." << std::endl;
}

inline double getRand(double l = 0, double r = 1)
{
  return (rand() % 1000000 / 1000000.)*(r - l) + l;
}

inline std::string escape(std::string str)
{
  std::string t;
  for (char c : str)
  {
    if (c == '\\'
      || c == '_'
      || c == '#')
      t += "\\";
    t += c;
  }
  return t;
}

inline void mergeMSS(std::map<std::string, std::string> &a, std::map<std::string, std::string> b)
{
  for (auto& itr : b)
    a[itr.first] = itr.second;
}

inline void mergeVSS(std::vector<std::pair<std::string, std::string>>&a, std::vector<std::pair<std::string, std::string>> b)
{
  a.insert(a.end(), b.begin(), b.end());
}

inline double getAverage(std::vector<double> dVec)
{
  return std::accumulate(dVec.begin(), dVec.end(), 0.0) / dVec.size();
}

inline std::vector<std::vector<std::string>> VSS2VVS(std::vector<std::pair<std::string, std::string>> data)
{
  std::vector<std::vector<std::string>> table;
  for (auto& i : data)
  {
    std::vector<std::string> row;
    row.clear();
    row.push_back(i.first);
    row.push_back(i.second);
    table.push_back(row);
  }
  return table;
}

inline std::vector<std::vector<std::string>> transVVS(std::vector<std::vector<std::string>> data)
{
  std::vector<std::vector<std::string>> newData;
  for (auto i = 0; i < data[0].size(); ++i)
  {
    std::vector<std::string> row;

    for (auto& j : data)
    {
      row.push_back(j[i]);
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
  fprintf(fout, "\\caption{%s}\n", escape(name).c_str());
  fprintf(fout, "\\begin{tabular}{|");

  for (auto i = 0; i < table[0].size(); ++i)
  {
    fprintf(fout, "c|");
  }
  fprintf(fout, "}\n");

  for (auto& i : table)
  {
    fprintf(fout, "\\hline\n");
    for (int j = 0; j < i.size(); j++)
    {
      if (j)
      {
        fprintf(fout, "&");
      }
      fprintf(fout, "%s", escape(i[j]).c_str());
    }
    fprintf(fout, "\\\\\n");
  }
  fprintf(fout, "\\hline\n");
  fprintf(fout, "\\end{tabular}\\end{table}\n");
}

inline void Tex_includeGraphics(const std::string& filename, const std::string& name, double width, FILE* fout)
{
  fprintf(fout, "\\begin{figure}[htb!]\n\\centering\n");
  fprintf(fout, "\\caption{%s}\n", escape(name).c_str());
  fprintf(fout, "\\includegraphics[width=%.1fin]{%s}\n", width, filename.c_str());
  fprintf(fout, "\\end{figure}\n");
}

inline std::vector<std::string> split(const std::string& str)
{
  std::vector<std::string> res;

  char s[10000];
  std::snprintf(s, sizeof s, "%s", str.c_str());
  for (char* p = strtok(s, ";"); p; p = strtok(NULL, ";"))
  {
    res.emplace_back(p);
  }
  return res;
}
#endif // BASIC_H
