#pragma once

#ifndef FRAMEWORK_H
#define FRAMEWORK_H

class Framework
{
public:
  virtual ~Framework() = default;
  std::map<std::string, std::string> config;
  void setConfig(std::map<std::string, std::string> config)
  {
    this->config = config;
  }

  virtual bool checkRequire() {
    return true;
  }

  virtual Graph calcNextLayerGraph(Graph cur, Community comm)
  {
    return cur;
  }

  virtual Graph calcLayerGraph(Graph cur, std::vector<Community> comms)
  {
    for (int i = 0; i < comms.size(); i++)
      cur = calcNextLayerGraph(cur, comms[i]);
    return cur;
  }

  virtual Graph calcLayerGraphAll(Graph cur, std::vector<Community> comms)
  {
    //cout<<"here"<< std::endl;
    //cur=calcLayerGraphAll(cur,comms);
    return cur;
  }
};
#endif // FRAMEWORK_H
