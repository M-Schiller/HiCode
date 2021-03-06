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

class Node
{
public:
  ~Node();
  Node();

  Node(int modulenr);

  std::vector<int> members; // If module, lists member nodes in module
  std::vector<std::pair<int, double>> links; // List of identities and link weight of connected nodes/modules

  double exit; // total weight of links to other nodes / modules
  double degree; // total degree of node / module
  int index; // the node / module identity
};

#endif
