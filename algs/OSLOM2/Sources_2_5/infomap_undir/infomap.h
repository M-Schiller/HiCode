#pragma once

#ifndef INFOMAP_H
#define INFOMAP_H

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "MersenneTwister.h"
#include "GreedyBase.h"
#include "Greedy.h"
#include "Node.h"
#define PI 3.14159265

unsigned stou(char *s);

class treeNode {
public:
  std::multimap<double, std::pair<int, std::string>, std::greater<>> members;
  std::multimap<double, treeNode, std::greater<>> nextLevel;
};

template <class T>
std::string to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

void cpyNode(Node *newNode, Node *oldNode) {
  newNode->index = oldNode->index;
  newNode->exit = oldNode->exit;
  newNode->degree = oldNode->degree;

  int Nmembers = oldNode->members.size();
  newNode->members = std::vector<int>(Nmembers);
  for (int i = 0; i < Nmembers; i++)
  {
    newNode->members[i] = oldNode->members[i];
  }
  int Nlinks = oldNode->links.size();
  newNode->links = std::vector<std::pair<int, double>>(Nlinks);
  for (int i = 0; i < Nlinks; i++)
  {
    newNode->links[i].first = oldNode->links[i].first;
    newNode->links[i].second = oldNode->links[i].second;
  }
}
#endif // INFOMAP_H
