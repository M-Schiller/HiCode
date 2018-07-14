#include "Node.h"

Node::Node()
  : selfLink(0.0)
  , exit(0.0)
  , size(0.0)
{
}

Node::Node(int nodenr, double tpweight)
  : selfLink(0.0)
  , teleportWeight(tpweight)
  , exit(0.0)
  , size(0.0)
  , index(nodenr)
{
  members.push_back(nodenr);
}
