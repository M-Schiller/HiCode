#pragma once

#ifndef GREEDYBASE_H
#define GREEDYBASE_H

#include "MersenneTwister.h"
#include <cstdio>
#include <vector>

// forward declaration
class Node;
class GreedyBase
{
public:
  GreedyBase() = default;
  virtual ~GreedyBase() = default;
  virtual void initiate() {};
  virtual void tune() {};
  virtual void calibrate() {};
  virtual void prepare(bool sort) {};
  virtual void level(Node ***, bool sort) {};
  virtual void move(bool &moved) {};
  virtual void determineMove(std::vector<int> &moveTo) {};
  virtual void eigenvector() {};
  int Nmod;
  int Nnode;
  int Nmember;
  int Ndanglings;

  double exit;
  double exitFlow;
  double exit_log_exit;
  double size_log_size;
  double nodeSize_log_nodeSize;

  double codeLength;

  Node **node;
  bool bottom;
  double alpha, beta;

protected:

  MTRand * R;
};

#endif
