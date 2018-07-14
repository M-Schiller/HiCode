#pragma once

#ifndef STANDARD_INCLUDE_INCLUDED
#define STANDARD_INCLUDE_INCLUDED

#include <cmath>
#include <iostream>
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <queue>
#include <list>

typedef std::deque<std::deque<int>> int_matrix;

#include "cast.cpp"
#include "print.cpp"
#include "random.cpp"
#include "combinatorics.cpp"
#include "histograms.cpp"
#include "deque_numeric.cpp"
#include "partition.cpp"
#include "mutual.cpp"
#include "pajek.cpp"

inline void systemCall(const std::string& cmd)
{
  std::cerr << "Call: " << cmd;
  const int res = system(cmd.c_str());
  std::cerr << " returned " << res << "." << std::endl;
}

#endif
