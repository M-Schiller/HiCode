#pragma once

#ifndef GREEDY_H
#define GREEDY_H

#include "MersenneTwister.h"
#include "GreedyBase.h"
#include "Node.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <stack>
#include <map>
#include <algorithm>

class Greedy : public GreedyBase
{
public:
  Greedy(MTRand *RR, int nnode, double deg, Node **node);
  virtual ~Greedy();
  virtual void initiate();
  virtual void calibrate();
  virtual void tune();
  virtual void prepare(bool sort);
  virtual void level(Node ***, bool sort);
  virtual void move(bool &moved);
  virtual void determineMove(std::vector<int> &moveTo);

  int Nempty;
  std::vector<int> mod_empty;

  std::vector<double> mod_exit;
  std::vector<double> mod_degree;
  std::vector<int> mod_members;

protected:
  double plogp(double d);
  std::vector<std::pair<int, double>>::iterator link;
  std::vector<int> modWnode;
};

#endif
