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
  Greedy(MTRand *RR, int nnode, Node **node, int nmembers);
  virtual ~Greedy();
  virtual void initiate();
  virtual void calibrate();
  virtual void tune();
  virtual void prepare(bool sort);
  virtual void level(Node ***, bool sort);
  virtual void move(bool &moved);
  virtual void determineMove(std::vector<int> &moveTo);
  virtual void eigenvector();

  std::vector<int> danglings;

  int Nempty;
  std::vector<int> mod_empty;

  std::vector<double> mod_exit;
  std::vector<double> mod_size;
  std::vector<double> mod_danglingSize;
  std::vector<double> mod_teleportWeight;
  std::vector<int> mod_members;

protected:

  std::vector<int> modSnode;
};

#endif
