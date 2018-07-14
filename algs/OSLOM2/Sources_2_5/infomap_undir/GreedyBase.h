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
  int Nmod;
  int Nnode;

  double degree;
  double invDegree;
  double log2;

  double exit;
  double exitDegree;
  double exit_log_exit;
  double degree_log_degree;
  double nodeDegree_log_nodeDegree;

  double codeLength;

  Node **node;

protected:

  MTRand * R;
};

#endif
