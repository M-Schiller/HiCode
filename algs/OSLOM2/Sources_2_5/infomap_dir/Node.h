#pragma once

#ifndef NODE_H
#define NODE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <stack>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cstring>
class Node {
public:

  Node();
  Node(int nodenr, double tpweight);
  std::vector<int> members;
  std::vector<std::pair<int, double>> inLinks;
  std::vector<std::pair<int, double>> outLinks;
  double selfLink;

  double teleportWeight;
  double danglingSize;
  double exit;
  double size;
  int index;
};

#endif
