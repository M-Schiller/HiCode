#pragma once

#ifndef PARTITION_INCLUDED
#define PARTITION_INCLUDED

#include <deque>
#include <fstream>
#include <set>
#include <map>

int get_partition_from_file(std::string s, std::deque<std::deque<int>> & M)
{
  M.clear();

  std::ifstream inb(s);

  std::string st;
  while (getline(inb, st))
  {
    std::deque<int> v;
    cast_string_to_doubles(st, v);
    sort(v.begin(), v.end());
    if (!v.empty())
      M.push_back(v);
  }
  return 0;
}

int get_partition_from_file_list(std::string s, std::deque<std::deque<int>> & ten)
{
  ten.clear();

  std::ifstream inb(s);

  std::multimap<int, int> id_value;
  std::map <int, int> value_com;
  std::string st;
  while (getline(inb, st)) {
    std::deque<int> v;
    cast_string_to_doubles(st, v);

    for (int i = 1; i < v.size(); i++) {
      value_com.insert(std::make_pair(v[i], value_com.size()));
      id_value.insert(std::make_pair(v[0], v[i]));
    }
  }

  int com = 0;
  for (auto it = value_com.begin(); it != value_com.end(); ++it)
    it->second = com++;

  std::deque <int> first;
  for (int i = 0; i < value_com.size(); i++)
    ten.push_back(first);

  for (auto itm = id_value.begin(); itm != id_value.end(); ++itm)
  {
    ten[value_com[itm->second]].push_back(itm->first);
  }

  return 0;
}

#endif
